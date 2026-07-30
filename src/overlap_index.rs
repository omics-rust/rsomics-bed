use std::collections::HashMap;
use std::io::Read;

use rsomics_common::{Result, RsomicsError};
use rsomics_intervals::{Interval, IntervalIndex, IntervalSet};

use crate::bed::{BedReader, BedRecord};

const BACKEND_END_EXCLUSIVE_LIMIT: u64 = i32::MAX as u64;

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
}

impl IndexedRecord {
    pub(crate) fn virtual_bounds(self) -> (u64, u64) {
        if self.start == self.end {
            // IndexedBed::load rejects both coordinate-domain edges first.
            (self.start - 1, self.end + 1)
        } else {
            (self.start, self.end)
        }
    }
}

#[derive(Default)]
struct ChromIndex {
    bounds: Vec<BoundGroup>,
    bound_ids: Vec<usize>,
    zero_ids: Vec<(u64, usize)>,
    coverage: Vec<(u64, u64)>,
}

struct BoundGroup {
    start: u64,
    end: u64,
    first_id: usize,
    last_id: usize,
}

impl IndexedBed {
    pub(crate) fn load(input: impl Read, label: &str) -> Result<Self> {
        Self::load_inner(input, label, false)
    }

    pub(crate) fn load_for_subtract(input: impl Read, label: &str) -> Result<Self> {
        Self::load_inner(input, label, true)
    }

    fn load_inner(input: impl Read, label: &str, build_coverage: bool) -> Result<Self> {
        let mut reader = BedReader::new(input);
        let mut records = Vec::new();
        let mut chrom_ranks = HashMap::new();
        let mut chrom_names = Vec::new();
        let mut pending_ids: Vec<Vec<usize>> = Vec::new();

        while let Some(record) = reader.next_coordinates()? {
            ensure_indexable(record.chrom, record.start, record.end, label)?;
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
            if build_coverage {
                let mut coverage: Vec<_> =
                    ids.iter().map(|&id| records[id].virtual_bounds()).collect();
                coverage.sort_unstable();
                for span in coverage {
                    if let Some(last) = chromosome.coverage.last_mut()
                        && span.0 <= last.1
                    {
                        last.1 = last.1.max(span.1);
                    } else {
                        chromosome.coverage.push(span);
                    }
                }
            }
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

        let index = IntervalIndex::build(&set);
        Ok(Self {
            index,
            records,
            chrom_ranks,
            chromosomes,
        })
    }

    pub(crate) fn ensure_query(&self, record: &BedRecord, label: &str) -> Result<()> {
        ensure_indexable(&record.chrom, record.start, record.end, label)
    }

    pub(crate) fn record(&self, id: usize) -> IndexedRecord {
        self.records[id]
    }

    pub(crate) fn intersection_candidates(&self, record: &BedRecord, ids: &mut Vec<usize>) {
        ids.clear();
        if record.start == record.end {
            let position = record.start;
            self.append_nonzero_ids(&record.chrom, position - 1, position + 1, ids);
            self.append_zero_ids(&record.chrom, position - 1, position + 1, ids);
        } else {
            self.append_nonzero_ids(&record.chrom, record.start, record.end, ids);
            self.append_zero_ids(&record.chrom, record.start, record.end, ids);
        }
        ids.sort_unstable();
    }

    pub(crate) fn coverage_overlaps(&self, record: &BedRecord, covers: &mut Vec<(u64, u64)>) {
        covers.clear();
        let Some(&rank) = self.chrom_ranks.get(&record.chrom) else {
            return;
        };
        let target = if record.start == record.end {
            (record.start - 1, record.end + 1)
        } else {
            (record.start, record.end)
        };
        let coverage = &self.chromosomes[rank].coverage;
        let first = coverage.partition_point(|&(_, high)| high <= target.0);
        for &(low, high) in &coverage[first..] {
            if low >= target.1 {
                break;
            }
            covers.push((target.0.max(low), target.1.min(high)));
        }
    }

    fn append_nonzero_ids(&self, chrom: &str, start: u64, end: u64, ids: &mut Vec<usize>) {
        let Some(&rank) = self.chrom_ranks.get(chrom) else {
            return;
        };
        let chromosome = &self.chromosomes[rank];
        self.index.for_each_overlap(chrom, start, end, |hit| {
            let group = chromosome
                .bounds
                .binary_search_by(|candidate| {
                    (candidate.start, candidate.end).cmp(&(hit.start, hit.end))
                })
                .ok()
                .map(|index| &chromosome.bounds[index])
                .expect("index coordinates originate from the bounds map");
            ids.extend_from_slice(&chromosome.bound_ids[group.first_id..group.last_id]);
        });
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

fn ensure_indexable(chrom: &str, start: u64, end: u64, label: &str) -> Result<()> {
    if start >= BACKEND_END_EXCLUSIVE_LIMIT || end > BACKEND_END_EXCLUSIVE_LIMIT {
        return Err(RsomicsError::InvalidInput(format!(
            "{label} interval {}:{}-{} exceeds the published interval-index backend maximum inclusive coordinate {}",
            chrom,
            start,
            end,
            BACKEND_END_EXCLUSIVE_LIMIT - 1
        )));
    }
    if start == 0 && end == 0 {
        return Err(RsomicsError::InvalidInput(format!(
            "{label} zero-length interval {}:0-0 cannot be represented after bedtools-compatible widening",
            chrom
        )));
    }
    Ok(())
}
