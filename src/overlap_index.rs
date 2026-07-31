use std::collections::HashMap;
use std::io::Read;

use coitrees::{COITree, Interval as CoitInterval, IntervalTree};
use rsomics_common::{Result, RsomicsError};

use crate::bed::{BedReader, BedRecord};

// COITrees returns SIMD metadata by reference and scalar metadata by value.
#[cfg(any(target_feature = "avx2", target_feature = "neon"))]
macro_rules! meta_id {
    ($node:ident) => {
        *$node.metadata
    };
}

#[cfg(not(any(target_feature = "avx2", target_feature = "neon")))]
macro_rules! meta_id {
    ($node:ident) => {
        $node.metadata
    };
}

// i32::MAX overflows the backend's SIMD layout.
const MAX_INDEX_COORDINATE: u64 = (i32::MAX - 1) as u64;

pub(crate) struct IndexedBed {
    records: Vec<IndexedRecord>,
    chrom_ranks: HashMap<String, usize>,
    chromosomes: Vec<ChromIndex>,
}

#[derive(Clone, Copy)]
pub(crate) struct IndexedRecord {
    pub(crate) start: u64,
    pub(crate) end: u64,
    overlap_start: u64,
    overlap_end: u64,
}

impl IndexedRecord {
    pub(crate) fn virtual_bounds(self) -> (u64, u64) {
        (self.overlap_start, self.overlap_end)
    }
}

struct ChromIndex {
    tree: COITree<usize, u32>,
    bounds: Vec<BoundGroup>,
    bound_ids: Vec<usize>,
    zero_ids: Vec<(u64, usize)>,
}

struct BoundGroup {
    first_id: usize,
    last_id: usize,
}

impl IndexedBed {
    pub(crate) fn load(input: impl Read, label: &str) -> Result<Self> {
        let mut reader = BedReader::new(input);
        let mut records = Vec::new();
        let mut chrom_ranks = HashMap::new();
        let mut pending_ids: Vec<Vec<usize>> = Vec::new();

        while let Some(record) = reader.next_coordinates()? {
            let (overlap_start, overlap_end) =
                overlap_bounds(record.chrom(), record.start(), record.end(), label)?;
            let rank = if let Some(&rank) = chrom_ranks.get(record.chrom()) {
                rank
            } else {
                let rank = pending_ids.len();
                chrom_ranks.insert(record.chrom().to_owned(), rank);
                pending_ids.push(Vec::new());
                rank
            };
            let id = records.len();
            records.push(IndexedRecord {
                start: record.start(),
                end: record.end(),
                overlap_start,
                overlap_end,
            });
            pending_ids[rank].push(id);
        }

        let mut chromosomes = Vec::with_capacity(pending_ids.len());
        for ids in &mut pending_ids {
            ids.sort_unstable_by_key(|&id| {
                let record = records[id];
                (record.start, record.end, id)
            });

            let mut nodes = Vec::new();
            let mut bounds = Vec::new();
            let mut bound_ids = Vec::new();
            let mut zero_ids = Vec::new();
            let mut cursor = 0;
            while cursor < ids.len() {
                let record = records[ids[cursor]];
                let mut group_end = cursor + 1;
                while group_end < ids.len() {
                    let candidate = records[ids[group_end]];
                    if (candidate.start, candidate.end) != (record.start, record.end) {
                        break;
                    }
                    group_end += 1;
                }

                if record.start == record.end {
                    zero_ids.extend(ids[cursor..group_end].iter().map(|&id| (record.start, id)));
                } else {
                    let first_id = bound_ids.len();
                    bound_ids.extend_from_slice(&ids[cursor..group_end]);
                    let last_id = bound_ids.len();
                    let group_id = bounds.len();
                    bounds.push(BoundGroup { first_id, last_id });
                    let (first, last) = index_bounds(record.start, record.end, label)?
                        .expect("non-empty parsed BED interval has inclusive index bounds");
                    nodes.push(CoitInterval::new(first, last, group_id));
                }
                cursor = group_end;
            }

            chromosomes.push(ChromIndex {
                tree: COITree::new(&nodes),
                bounds,
                bound_ids,
                zero_ids,
            });
        }

        Ok(Self {
            records,
            chrom_ranks,
            chromosomes,
        })
    }

    pub(crate) fn record(&self, id: usize) -> IndexedRecord {
        self.records[id]
    }

    pub(crate) fn intersection_candidates(
        &self,
        record: &BedRecord,
        label: &str,
        ids: &mut Vec<usize>,
    ) -> Result<()> {
        ids.clear();
        let (start, end) = overlap_bounds(record.chrom(), record.start(), record.end(), label)?;
        self.append_nonzero_ids(record.chrom(), start, end, label, ids)?;
        self.append_zero_ids(record.chrom(), start, end, ids);
        ids.sort_unstable();
        Ok(())
    }

    fn append_nonzero_ids(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        label: &str,
        ids: &mut Vec<usize>,
    ) -> Result<()> {
        let Some(&rank) = self.chrom_ranks.get(chrom) else {
            return Ok(());
        };
        let Some((first, last)) = index_bounds(start, end, label)? else {
            return Ok(());
        };
        let chromosome = &self.chromosomes[rank];
        chromosome.tree.query(first, last, |node| {
            let group = &chromosome.bounds[meta_id!(node)];
            ids.extend_from_slice(&chromosome.bound_ids[group.first_id..group.last_id]);
        });
        Ok(())
    }

    fn append_zero_ids(&self, chrom: &str, low: u64, high: u64, ids: &mut Vec<usize>) {
        let Some(&rank) = self.chrom_ranks.get(chrom) else {
            return;
        };
        let positions = &self.chromosomes[rank].zero_ids;
        let first = positions.partition_point(|&(position, _)| position < low);
        let last = positions.partition_point(|&(position, _)| position <= high);
        ids.extend(positions[first..last].iter().map(|&(_, id)| id));
    }
}

pub(crate) struct CoverageBed {
    chrom_ranks: HashMap<String, usize>,
    coverage: Vec<Vec<(u64, u64)>>,
}

impl CoverageBed {
    pub(crate) fn load(input: impl Read, label: &str) -> Result<Self> {
        let mut reader = BedReader::new(input);
        let mut chrom_ranks = HashMap::new();
        let mut coverage: Vec<Vec<(u64, u64)>> = Vec::new();

        while let Some(record) = reader.next_coordinates()? {
            let span = overlap_bounds(record.chrom(), record.start(), record.end(), label)?;
            let rank = if let Some(&rank) = chrom_ranks.get(record.chrom()) {
                rank
            } else {
                let rank = coverage.len();
                chrom_ranks.insert(record.chrom().to_owned(), rank);
                coverage.push(Vec::new());
                rank
            };
            coverage[rank].push(span);
        }

        for spans in &mut coverage {
            spans.sort_unstable();
            let mut merged: Vec<(u64, u64)> = Vec::with_capacity(spans.len());
            for span in spans.drain(..) {
                if let Some(last) = merged.last_mut()
                    && span.0 <= last.1
                {
                    last.1 = last.1.max(span.1);
                } else {
                    merged.push(span);
                }
            }
            *spans = merged;
        }

        Ok(Self {
            chrom_ranks,
            coverage,
        })
    }

    pub(crate) fn overlaps(
        &self,
        record: &BedRecord,
        label: &str,
        overlaps: &mut Vec<(u64, u64)>,
    ) -> Result<()> {
        overlaps.clear();
        let target = overlap_bounds(record.chrom(), record.start(), record.end(), label)?;
        let Some(&rank) = self.chrom_ranks.get(record.chrom()) else {
            return Ok(());
        };
        let coverage = &self.coverage[rank];
        let first = coverage.partition_point(|&(_, high)| high <= target.0);
        for &(low, high) in &coverage[first..] {
            if low >= target.1 {
                break;
            }
            overlaps.push((target.0.max(low), target.1.min(high)));
        }
        Ok(())
    }
}

fn overlap_bounds(chrom: &str, start: u64, end: u64, label: &str) -> Result<(u64, u64)> {
    if start == 0 && end == 0 {
        return Err(RsomicsError::InvalidInput(format!(
            "{label} zero-length interval {}:0-0 cannot be represented after bedtools-compatible widening",
            chrom
        )));
    }
    if start == end {
        let end = end.checked_add(1).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "{label} zero-length interval {chrom}:{start}-{end} cannot be widened beyond u64"
            ))
        })?;
        Ok((start - 1, end))
    } else {
        Ok((start, end))
    }
}

fn index_bounds(start: u64, end: u64, label: &str) -> Result<Option<(i32, i32)>> {
    if start >= end {
        return Ok(None);
    }
    let last = end - 1;
    for coordinate in [start, last] {
        if coordinate > MAX_INDEX_COORDINATE {
            return Err(RsomicsError::InvalidInput(format!(
                "{label} interval index: coordinate {coordinate} exceeds the interval-index backend maximum {MAX_INDEX_COORDINATE}"
            )));
        }
    }
    Ok(Some((start as i32, last as i32)))
}
