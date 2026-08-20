use std::io::Read;
use std::ops::Range;

use rsomics_common::{Context, Result};

use crate::bed::{BedSliceReader, Strand, invalid, required_field, required_strand};
use crate::interval_index::{IndexedInterval, IntervalIndex, IntervalIndexBuilder};

pub(crate) struct RelationBed {
    raw: Vec<u8>,
    raw_ranges: Vec<Range<usize>>,
    line_numbers: Vec<usize>,
    index: IntervalIndex,
    field_count: usize,
    chromosomes: Vec<ChromOrder>,
}

struct ChromOrder {
    start_ids: Vec<usize>,
    end_ids: Vec<usize>,
}

pub(crate) struct RelationRecord<'a> {
    raw: &'a [u8],
    line_number: usize,
}

impl<'a> RelationRecord<'a> {
    pub(crate) fn raw(&self) -> &'a [u8] {
        self.raw
    }

    pub(crate) fn name(&self, label: &str) -> Result<&'a [u8]> {
        required_field(self.raw, self.line_number, 3, label, "name")
    }

    pub(crate) fn strand(&self, label: &str) -> Result<Strand> {
        required_strand(self.raw, self.line_number, label)
    }

    #[cfg(test)]
    fn write_raw(&self, output: &mut dyn std::io::Write) -> Result<()> {
        output
            .write_all(self.raw)
            .rs_context("writing BED record")?;
        output.write_all(b"\n").rs_context("writing BED record")
    }
}

impl RelationBed {
    pub(crate) fn load(mut input: impl Read, label: &str) -> Result<Self> {
        let mut raw = Vec::new();
        input
            .read_to_end(&mut raw)
            .rs_context(format!("reading {label} BED input"))?;
        let mut reader = BedSliceReader::new(&raw);
        let mut builder = IntervalIndexBuilder::new();
        let mut raw_ranges = Vec::new();
        let mut line_numbers = Vec::new();
        let mut field_width = None;

        while let Some(record) = reader.next_record()? {
            let count = record.field_count();
            if let Some((expected, first_line)) = field_width {
                if count != expected {
                    return Err(invalid(format!(
                        "{label} BED line {} has {count} fields; line {first_line} has {expected}",
                        record.line_number()
                    )));
                }
            } else {
                field_width = Some((count, record.line_number()));
            }

            builder.push(record.chrom(), record.start(), record.end(), label)?;
            raw_ranges.push(record.raw_range());
            line_numbers.push(record.line_number());
        }
        let (index, start_orders) = builder.finish_with_start_order(label)?;
        let chromosomes = start_orders
            .into_iter()
            .map(|start_ids| {
                let mut end_ids = start_ids.clone();
                end_ids.sort_unstable_by_key(|&id| {
                    let interval = index.record(id);
                    let (start, end) = interval.virtual_bounds();
                    (end, start, id)
                });
                ChromOrder { start_ids, end_ids }
            })
            .collect();
        Ok(Self {
            raw,
            raw_ranges,
            line_numbers,
            index,
            field_count: field_width.map_or(3, |(count, _)| count),
            chromosomes,
        })
    }

    pub(crate) fn record(&self, id: usize) -> RelationRecord<'_> {
        RelationRecord {
            raw: &self.raw[self.raw_ranges[id].clone()],
            line_number: self.line_numbers[id],
        }
    }

    pub(crate) fn interval(&self, id: usize) -> IndexedInterval {
        self.index.record(id)
    }

    pub(crate) fn field_count(&self) -> usize {
        self.field_count
    }

    pub(crate) fn len(&self) -> usize {
        self.raw_ranges.len()
    }

    pub(crate) fn checked_strands(&self, label: &str) -> Result<Vec<Strand>> {
        (0..self.len())
            .map(|id| self.record(id).strand(label))
            .collect()
    }

    pub(crate) fn range_candidates(
        &self,
        chrom: &str,
        start: u64,
        end: u64,
        label: &str,
        ids: &mut Vec<usize>,
    ) -> Result<()> {
        self.index.range_candidates(chrom, start, end, label, ids)
    }

    pub(crate) fn left_candidates(
        &self,
        chrom: &str,
        boundary: u64,
    ) -> impl Iterator<Item = usize> + '_ {
        let ids = self.chromosome(chrom).map_or(&[][..], |order| {
            let last = order
                .end_ids
                .partition_point(|&id| self.index.record(id).virtual_bounds().1 <= boundary);
            &order.end_ids[..last]
        });
        ids.iter().rev().copied()
    }

    pub(crate) fn right_candidates(
        &self,
        chrom: &str,
        boundary: u64,
    ) -> impl Iterator<Item = usize> + '_ {
        let ids = self.chromosome(chrom).map_or(&[][..], |order| {
            let first = order
                .start_ids
                .partition_point(|&id| self.index.record(id).virtual_bounds().0 < boundary);
            &order.start_ids[first..]
        });
        ids.iter().copied()
    }

    fn chromosome(&self, chrom: &str) -> Option<&ChromOrder> {
        self.index
            .chromosome_rank(chrom)
            .map(|rank| &self.chromosomes[rank])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn full_records_and_duplicate_ids_are_preserved() {
        let relation = RelationBed::load(
            &b"chr1\t30\t40\tfirst\n\
              chr1\t10\t20\tsecond\n\
              chr1\t30\t40\tlast\n"[..],
            "B",
        )
        .unwrap();
        let mut output = Vec::new();
        relation.record(2).write_raw(&mut output).unwrap();
        assert_eq!(output, b"chr1\t30\t40\tlast\n");

        let mut ids = Vec::new();
        relation
            .range_candidates("chr1", 0, 100, "A", &mut ids)
            .unwrap();
        assert_eq!(ids, [0, 1, 2]);
    }

    #[test]
    fn variable_b_widths_report_both_physical_lines() {
        let error = RelationBed::load(&b"# header\nchr1\t10\t20\nchr1\t30\t40\tname\n"[..], "B")
            .err()
            .expect("variable B widths must fail");
        assert!(
            error
                .to_string()
                .contains("B BED line 3 has 4 fields; line 2 has 3"),
            "{error}"
        );
    }

    #[test]
    fn empty_b_and_absent_chromosomes_are_empty() {
        let relation = RelationBed::load(&b""[..], "B").unwrap();
        let mut ids = vec![99];
        relation
            .range_candidates("chr1", 10, 20, "A", &mut ids)
            .unwrap();
        assert!(ids.is_empty());
    }

    #[test]
    fn range_queries_share_zero_length_virtual_bounds() {
        let relation = RelationBed::load(
            &b"chr1\t9\t9\tleft\nchr1\t10\t20\tspan\nchr1\t11\t11\tright\n"[..],
            "B",
        )
        .unwrap();
        let mut ids = Vec::new();
        relation
            .range_candidates("chr1", 9, 11, "A", &mut ids)
            .unwrap();
        assert_eq!(ids, [0, 1, 2]);
    }

    #[test]
    fn directional_candidates_advance_outward_from_virtual_boundaries() {
        let relation = RelationBed::load(
            &b"chr1\t0\t5\tleft-first\n\
              chr1\t2\t5\tleft-second\n\
              chr1\t9\t9\tleft-point\n\
              chr1\t25\t30\tright-first\n\
              chr1\t25\t35\tright-second\n"[..],
            "B",
        )
        .unwrap();
        assert_eq!(relation.field_count(), 4);
        assert_eq!(relation.interval(2).virtual_bounds(), (8, 10));
        assert_eq!(
            relation.left_candidates("chr1", 10).collect::<Vec<_>>(),
            [2, 1, 0]
        );
        assert_eq!(
            relation.right_candidates("chr1", 20).collect::<Vec<_>>(),
            [3, 4]
        );
        assert_eq!(relation.left_candidates("chr2", 10).count(), 0);
        assert_eq!(relation.right_candidates("chr2", 20).count(), 0);
    }

    #[test]
    fn start_order_remains_monotonic_across_zero_length_widening() {
        let relation = RelationBed::load(
            &b"chr1\t9\t11\tspan-left\n\
              chr1\t10\t10\tpoint\n\
              chr1\t10\t12\tspan-right\n"[..],
            "B",
        )
        .unwrap();
        assert_eq!(
            relation.right_candidates("chr1", 9).collect::<Vec<_>>(),
            [0, 1, 2]
        );
        assert_eq!(
            relation.right_candidates("chr1", 10).collect::<Vec<_>>(),
            [2]
        );
    }
}
