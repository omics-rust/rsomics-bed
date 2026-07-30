use std::collections::HashMap;
use std::io::Read;

use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::{Interval, IntervalIndex, IntervalIndexError, IntervalSet};

use crate::bed::{BedReader, BedRecord};

pub(crate) struct IndexedBed {
    index: IntervalIndex,
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

#[derive(Default)]
struct ChromIndex {
    bounds: Vec<BoundGroup>,
    bound_ids: Vec<usize>,
    zero_ids: Vec<(u64, usize)>,
}

struct BoundGroup {
    start: u64,
    end: u64,
    first_id: usize,
    last_id: usize,
}

impl IndexedBed {
    pub(crate) fn load(input: impl Read, label: &str) -> Result<Self> {
        let mut reader = BedReader::new(input);
        let mut records = Vec::new();
        let mut chrom_ranks = HashMap::new();
        let mut chrom_names = Vec::new();
        let mut pending_ids: Vec<Vec<usize>> = Vec::new();

        while let Some(record) = reader.next_coordinates()? {
            let (overlap_start, overlap_end) =
                overlap_bounds(record.chrom, record.start, record.end, label)?;
            let rank = if let Some(&rank) = chrom_ranks.get(record.chrom) {
                rank
            } else {
                let rank = chrom_names.len();
                chrom_ranks.insert(record.chrom.to_owned(), rank);
                chrom_names.push(record.chrom.to_owned());
                pending_ids.push(Vec::new());
                rank
            };
            let id = records.len();
            records.push(IndexedRecord {
                start: record.start,
                end: record.end,
                overlap_start,
                overlap_end,
            });
            pending_ids[rank].push(id);
        }

        let mut set = IntervalSet::new();
        let mut chromosomes = Vec::with_capacity(chrom_names.len());
        for (rank, ids) in pending_ids.iter_mut().enumerate() {
            ids.sort_unstable_by_key(|&id| {
                let record = records[id];
                (record.start, record.end, id)
            });
            let mut chromosome = ChromIndex::default();
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
                    chromosome
                        .zero_ids
                        .extend(ids[cursor..group_end].iter().map(|&id| (record.start, id)));
                } else {
                    let first_id = chromosome.bound_ids.len();
                    chromosome
                        .bound_ids
                        .extend_from_slice(&ids[cursor..group_end]);
                    let last_id = chromosome.bound_ids.len();
                    chromosome.bounds.push(BoundGroup {
                        start: record.start,
                        end: record.end,
                        first_id,
                        last_id,
                    });
                    set.push(
                        Interval::new(chrom_names[rank].clone(), record.start, record.end)
                            .expect("BED parser rejects inverted intervals"),
                    );
                }
                cursor = group_end;
            }
            chromosomes.push(chromosome);
        }

        let index = IntervalIndex::try_build(&set).map_err(|error| index_error(label, error))?;
        Ok(Self {
            index,
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
        let (start, end) = overlap_bounds(&record.chrom, record.start, record.end, label)?;
        self.append_nonzero_ids(&record.chrom, start, end, label, ids)?;
        self.append_zero_ids(&record.chrom, start, end, ids);
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
        let chromosome = &self.chromosomes[rank];
        self.index
            .try_for_each_overlap(chrom, start, end, |hit| {
                let group = chromosome
                    .bounds
                    .binary_search_by(|candidate| {
                        (candidate.start, candidate.end).cmp(&(hit.start, hit.end))
                    })
                    .ok()
                    .map(|index| &chromosome.bounds[index])
                    .expect("index coordinates originate from the bounds map");
                ids.extend_from_slice(&chromosome.bound_ids[group.first_id..group.last_id]);
            })
            .map_err(|error| index_error(label, error))
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
            let span = overlap_bounds(record.chrom, record.start, record.end, label)?;
            let rank = if let Some(&rank) = chrom_ranks.get(record.chrom) {
                rank
            } else {
                let rank = coverage.len();
                chrom_ranks.insert(record.chrom.to_owned(), rank);
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
        let target = overlap_bounds(&record.chrom, record.start, record.end, label)?;
        let Some(&rank) = self.chrom_ranks.get(&record.chrom) else {
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

fn index_error(label: &str, error: IntervalIndexError) -> RsomicsError {
    RsomicsError::InvalidInput(format!("{label} interval index: {error}"))
}
