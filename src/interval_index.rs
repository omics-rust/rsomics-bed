use std::collections::HashMap;

use coitrees::{COITree, Interval as CoitInterval, IntervalTree};
use rsomics_common::{Result, RsomicsError};

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

const MAX_INDEX_COORDINATE: u64 = (i32::MAX - 1) as u64;

pub(crate) struct IntervalIndexBuilder {
    records: Vec<IndexedInterval>,
    chrom_ranks: HashMap<String, usize>,
    pending_ids: Vec<Vec<usize>>,
}

impl IntervalIndexBuilder {
    pub(crate) fn new() -> Self {
        Self {
            records: Vec::new(),
            chrom_ranks: HashMap::new(),
            pending_ids: Vec::new(),
        }
    }

    pub(crate) fn push(&mut self, chrom: &str, start: u64, end: u64, label: &str) -> Result<usize> {
        overlap_bounds(chrom, start, end, label)?;
        let rank = if let Some(&rank) = self.chrom_ranks.get(chrom) {
            rank
        } else {
            let rank = self.pending_ids.len();
            self.chrom_ranks.insert(chrom.to_owned(), rank);
            self.pending_ids.push(Vec::new());
            rank
        };
        let id = self.records.len();
        self.records.push(IndexedInterval { start, end });
        self.pending_ids[rank].push(id);
        Ok(id)
    }

    pub(crate) fn finish(self, label: &str) -> Result<IntervalIndex> {
        self.finish_with_start_order(label).map(|(index, _)| index)
    }

    pub(crate) fn finish_with_start_order(
        self,
        label: &str,
    ) -> Result<(IntervalIndex, Vec<Vec<usize>>)> {
        let Self {
            records,
            chrom_ranks,
            mut pending_ids,
        } = self;
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

        Ok((
            IntervalIndex {
                records,
                chrom_ranks,
                chromosomes,
            },
            pending_ids,
        ))
    }
}

pub(crate) struct IntervalIndex {
    records: Vec<IndexedInterval>,
    chrom_ranks: HashMap<String, usize>,
    chromosomes: Vec<ChromIndex>,
}

#[derive(Clone, Copy)]
pub(crate) struct IndexedInterval {
    pub(crate) start: u64,
    pub(crate) end: u64,
}

impl IndexedInterval {
    pub(crate) fn virtual_bounds(self) -> (u64, u64) {
        if self.start == self.end {
            (self.start - 1, self.end + 1)
        } else {
            (self.start, self.end)
        }
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

impl IntervalIndex {
    pub(crate) fn chromosome_rank(&self, chrom: &str) -> Option<usize> {
        self.chrom_ranks.get(chrom).copied()
    }

    pub(crate) fn record(&self, id: usize) -> IndexedInterval {
        self.records[id]
    }

    pub(crate) fn intersection_candidates(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        label: &str,
        ids: &mut Vec<usize>,
    ) -> Result<()> {
        let (start, end) = overlap_bounds(chrom, start, end, label)?;
        self.range_candidates(chrom, start, end, label, ids)
    }

    pub(crate) fn range_candidates(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        label: &str,
        ids: &mut Vec<usize>,
    ) -> Result<()> {
        ids.clear();
        if start >= end {
            return Ok(());
        }
        self.append_nonzero_ids(chrom, start, end, label, ids)?;
        self.append_zero_ids(chrom, start, end, ids);
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

pub(crate) fn overlap_bounds(chrom: &str, start: u64, end: u64, label: &str) -> Result<(u64, u64)> {
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

#[cfg(test)]
mod tests {
    use super::*;

    fn build(records: &[(&str, u64, u64)]) -> Result<IntervalIndex> {
        let mut builder = IntervalIndexBuilder::new();
        for &(chrom, start, end) in records {
            builder.push(chrom, start, end, "B")?;
        }
        builder.finish("B")
    }

    #[test]
    fn duplicate_and_zero_length_candidates_keep_record_ids() {
        let index = build(&[
            ("chr1", 10, 20),
            ("chr1", 30, 30),
            ("chr1", 10, 20),
            ("chr1", 40, 50),
        ])
        .unwrap();
        let mut ids = Vec::new();
        index
            .intersection_candidates("chr1", 9, 31, "A", &mut ids)
            .unwrap();
        assert_eq!(ids, [0, 1, 2]);
        assert_eq!(index.record(1).virtual_bounds(), (29, 31));
    }

    #[test]
    fn absent_chromosomes_return_no_candidates() {
        let index = build(&[("chr1", 10, 20)]).unwrap();
        let mut ids = vec![99];
        index
            .intersection_candidates("chr2", 10, 20, "A", &mut ids)
            .unwrap();
        assert!(ids.is_empty());
    }

    #[test]
    fn empty_query_ranges_return_no_zero_length_candidates() {
        let index = build(&[("chr1", 10, 10)]).unwrap();
        let mut ids = vec![99];
        index
            .range_candidates("chr1", 10, 10, "A", &mut ids)
            .unwrap();
        assert!(ids.is_empty());
    }

    #[test]
    fn coordinate_backend_boundary_is_checked() {
        let maximum = (i32::MAX - 1) as u64;
        let index = build(&[("chr1", maximum, maximum + 1)]).unwrap();
        assert_eq!(index.record(0).start, maximum);

        let error = build(&[("chr1", maximum + 1, maximum + 2)])
            .err()
            .expect("coordinate outside the backend must fail");
        assert!(error.to_string().contains("B interval index"), "{error}");
        assert!(error.to_string().contains("2147483647"), "{error}");
    }

    #[test]
    fn zero_length_origin_is_rejected() {
        let error = build(&[("chr1", 0, 0)])
            .err()
            .expect("origin zero-length record must fail");
        assert!(
            error.to_string().contains("B zero-length interval"),
            "{error}"
        );
    }
}
