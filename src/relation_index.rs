use std::io::Read;

use rsomics_common::Result;

use crate::bed::{BedReader, BedRecord, invalid};
use crate::interval_index::{IntervalIndex, IntervalIndexBuilder};

pub(crate) struct RelationBed {
    records: Vec<BedRecord>,
    index: IntervalIndex,
}

impl RelationBed {
    pub(crate) fn load(input: impl Read, label: &str) -> Result<Self> {
        let mut reader = BedReader::new(input);
        let mut builder = IntervalIndexBuilder::new();
        let mut records = Vec::new();
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
            records.push(record);
        }

        let index = builder.finish(label)?;
        Ok(Self { records, index })
    }

    pub(crate) fn record(&self, id: usize) -> &BedRecord {
        &self.records[id]
    }

    pub(crate) fn len(&self) -> usize {
        self.records.len()
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
}
