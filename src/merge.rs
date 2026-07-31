//! Streaming merge of sorted BED intervals.

use std::collections::HashSet;
use std::io::{BufWriter, Read, Write};

use rsomics_common::{Context, Result, RsomicsError};

use crate::bed::{BedReader, BedRecord, virtual_bounds};

/// Merge overlapping or touching intervals from sorted, chromosome-grouped BED.
///
/// Output is BED3. Zero-length records use bedtools-compatible virtual
/// widening when that widening remains inside the nonnegative `u64` domain.
///
/// # Errors
///
/// Returns an error for malformed, unsorted, or ungrouped BED; input/output
/// failures; or zero-length widening below zero or beyond `u64`.
pub fn merge(input: impl Read, output: impl Write) -> Result<()> {
    let mut reader = BedReader::new(input);
    let mut output = BufWriter::new(output);
    let mut closed = HashSet::new();
    let mut current: Option<Cluster> = None;
    let mut last_start = 0_u64;

    while let Some(record) = reader.next_record()? {
        match current.as_mut() {
            Some(cluster) if cluster.chrom == record.chrom() => {
                if record.start() < last_start {
                    return Err(RsomicsError::InvalidInput(format!(
                        "merge input is not sorted: {}:{} follows start {}",
                        record.chrom(),
                        record.start(),
                        last_start
                    )));
                }
                last_start = record.start();
                if cluster.absorb(continuation_bounds(&record)?) {
                    continue;
                }
                cluster.write(&mut output)?;
                *cluster = Cluster::from_record(&record, virtual_bounds(&record, "merge")?);
            }
            Some(cluster) => {
                cluster.write(&mut output)?;
                closed.insert(cluster.chrom.clone());
                if closed.contains(record.chrom()) {
                    return Err(RsomicsError::InvalidInput(format!(
                        "merge input is not grouped by chromosome: {:?} reappears",
                        record.chrom()
                    )));
                }
                last_start = record.start();
                *cluster = Cluster::from_record(&record, virtual_bounds(&record, "merge")?);
            }
            None => {
                last_start = record.start();
                current = Some(Cluster::from_record(
                    &record,
                    virtual_bounds(&record, "merge")?,
                ));
            }
        }
    }
    if let Some(cluster) = current {
        cluster.write(&mut output)?;
    }
    output.flush().rs_context("flushing merged BED output")
}

fn continuation_bounds(record: &BedRecord) -> Result<(u64, u64)> {
    if record.start() != record.end() {
        return Ok((record.start(), record.end()));
    }
    let high = record.end().checked_add(1).ok_or_else(|| {
        RsomicsError::InvalidInput(format!(
            "merge interval {}:{}-{} widens beyond u64",
            record.chrom(),
            record.start(),
            record.end()
        ))
    })?;
    Ok((record.start().saturating_sub(1), high))
}

struct Cluster {
    chrom: String,
    low: u64,
    high: u64,
    single_zero: Option<(u64, u64)>,
}

impl Cluster {
    fn from_record(record: &BedRecord, (low, high): (u64, u64)) -> Self {
        Self {
            chrom: record.chrom().to_owned(),
            low,
            high,
            single_zero: (record.start() == record.end()).then_some((record.start(), record.end())),
        }
    }

    fn absorb(&mut self, (low, high): (u64, u64)) -> bool {
        if low > self.high {
            return false;
        }
        self.high = self.high.max(high);
        self.single_zero = None;
        true
    }

    fn write(&self, output: &mut dyn Write) -> Result<()> {
        match self.single_zero {
            Some((start, end)) => writeln!(output, "{}\t{start}\t{end}", self.chrom),
            None => writeln!(output, "{}\t{}\t{}", self.chrom, self.low, self.high),
        }
        .rs_context("writing merged BED record")
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn touching_and_zero_length_records_merge() {
        let input = b"chr1\t10\t20\nchr1\t20\t20\nchr1\t30\t40\n";
        let mut output = Vec::new();
        merge(&input[..], &mut output).unwrap();
        assert_eq!(output, b"chr1\t10\t21\nchr1\t30\t40\n");
    }

    #[test]
    fn later_zero_length_record_does_not_move_cluster_start_backwards() {
        let input = b"chr1\t2\t13\nchr1\t2\t2\n";
        let mut output = Vec::new();
        merge(&input[..], &mut output).unwrap();
        assert_eq!(output, b"chr1\t2\t13\n");
    }

    #[test]
    fn first_zero_length_record_defines_the_widened_cluster_start() {
        let input = b"chr1\t2\t2\nchr1\t2\t13\n";
        let mut output = Vec::new();
        merge(&input[..], &mut output).unwrap();
        assert_eq!(output, b"chr1\t1\t13\n");
    }

    #[test]
    fn later_origin_zero_length_record_does_not_require_negative_output() {
        let input = b"chr1\t0\t5\nchr1\t0\t0\n";
        let mut output = Vec::new();
        merge(&input[..], &mut output).unwrap();
        assert_eq!(output, b"chr1\t0\t5\n");
    }

    #[test]
    fn unsorted_cluster_fails_loud() {
        let error = merge(
            &b"chr1\t10\t100\nchr1\t20\t30\nchr1\t15\t16\n"[..],
            Vec::new(),
        )
        .unwrap_err();
        assert!(error.to_string().contains("not sorted"), "{error}");
    }

    #[test]
    fn maximum_zero_length_interval_fails_before_widening() {
        let input = format!("chr1\t{}\t{}\n", u64::MAX, u64::MAX);
        let error = merge(input.as_bytes(), Vec::new()).unwrap_err();
        assert!(error.to_string().contains("widens beyond u64"), "{error}");
    }

    #[test]
    fn origin_zero_length_interval_fails_instead_of_emitting_negative_bed() {
        let error = merge(&b"chr1\t0\t0\nchr1\t0\t5\n"[..], Vec::new()).unwrap_err();
        assert!(error.to_string().contains("widens below"), "{error}");
    }
}
