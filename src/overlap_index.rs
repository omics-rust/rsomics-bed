use std::collections::HashMap;
use std::io::Read;

use rsomics_common::Result;

use crate::bed::{BedReader, BedRecord};
use crate::interval_index::{IndexedInterval, IntervalIndex, IntervalIndexBuilder, overlap_bounds};

pub(crate) struct IndexedBed {
    index: IntervalIndex,
}

impl IndexedBed {
    pub(crate) fn load(input: impl Read, label: &str) -> Result<Self> {
        let mut reader = BedReader::new(input);
        let mut builder = IntervalIndexBuilder::new();

        while let Some(record) = reader.next_coordinates()? {
            builder.push(record.chrom(), record.start(), record.end(), label)?;
        }
        Ok(Self {
            index: builder.finish(label)?,
        })
    }

    pub(crate) fn record(&self, id: usize) -> IndexedInterval {
        self.index.record(id)
    }

    pub(crate) fn intersection_candidates(
        &self,
        record: &BedRecord,
        label: &str,
        ids: &mut Vec<usize>,
    ) -> Result<()> {
        self.index
            .intersection_candidates(record.chrom(), record.start(), record.end(), label, ids)
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
