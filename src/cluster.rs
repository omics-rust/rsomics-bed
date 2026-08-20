//! Streaming BED clustering with BEDTools-compatible distance and strand order.

use std::collections::HashSet;
use std::io::{BufWriter, Read, Write};

use rsomics_common::{Context, Result};

use crate::bed::{BedReader, BedRecord, Strand, invalid};

/// Options for [`cluster`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ClusterOptions {
    /// Maximum nonnegative gap between records in one cluster.
    pub distance: u64,
    /// Cluster forward and reverse strands independently.
    pub same_strand: bool,
}

/// Append a one-based cluster ID to each sorted BED record.
///
/// Unstranded output preserves input order. Same-strand output follows
/// BEDTools by emitting forward then reverse records within each chromosome.
///
/// # Errors
///
/// Returns an error for malformed or unsorted BED, invalid required strand
/// fields, distance overflow, or input/output failure.
pub fn cluster(input: impl Read, output: impl Write, options: ClusterOptions) -> Result<()> {
    let mut reader = BedReader::new(input);
    let mut output = BufWriter::new(output);
    if options.same_strand {
        cluster_by_strand(&mut reader, &mut output, options.distance)?;
    } else {
        cluster_stream(&mut reader, &mut output, options.distance)?;
    }
    output.flush().rs_context("flushing cluster BED output")
}

fn cluster_stream(
    reader: &mut BedReader<impl std::io::BufRead>,
    output: &mut dyn Write,
    distance: u64,
) -> Result<()> {
    let mut order = InputOrder::default();
    let mut active: Option<ActiveCluster> = None;
    let mut cluster_id = 0_u64;

    while let Some(record) = reader.next_record()? {
        let new_chromosome = order.accept(&record)?;
        let joins = !new_chromosome
            && active
                .as_ref()
                .is_some_and(|cluster| record.start() <= cluster.reach);
        if joins {
            let active = active
                .as_mut()
                .expect("join decision requires an active cluster");
            active.extend(&record, distance)?;
        } else {
            cluster_id += 1;
            active = Some(ActiveCluster::new(&record, distance)?);
        }
        record.write_column(output, cluster_id)?;
    }
    Ok(())
}

fn cluster_by_strand(
    reader: &mut BedReader<impl std::io::BufRead>,
    output: &mut dyn Write,
    distance: u64,
) -> Result<()> {
    let mut order = InputOrder::default();
    let mut records = Vec::new();
    let mut cluster_id = 0_u64;

    while let Some(record) = reader.next_record()? {
        if order.accept(&record)? && !records.is_empty() {
            write_stranded_chromosome(&records, output, distance, &mut cluster_id)?;
            records.clear();
        }
        records.push(record);
    }
    if !records.is_empty() {
        write_stranded_chromosome(&records, output, distance, &mut cluster_id)?;
    }
    Ok(())
}

fn write_stranded_chromosome(
    records: &[BedRecord],
    output: &mut dyn Write,
    distance: u64,
    cluster_id: &mut u64,
) -> Result<()> {
    let strands = records
        .iter()
        .map(|record| record.strand("cluster"))
        .collect::<Result<Vec<_>>>()?;

    for selected in [Strand::Forward, Strand::Reverse] {
        let mut active: Option<ActiveCluster> = None;
        for (record, strand) in records.iter().zip(&strands) {
            if *strand != selected {
                continue;
            }
            if active
                .as_ref()
                .is_some_and(|cluster| record.start() <= cluster.reach)
            {
                active
                    .as_mut()
                    .expect("join decision requires an active cluster")
                    .extend(record, distance)?;
            } else {
                *cluster_id += 1;
                active = Some(ActiveCluster::new(record, distance)?);
            }
            record.write_column(output, *cluster_id)?;
        }
    }
    Ok(())
}

struct ActiveCluster {
    end: u64,
    reach: u64,
}

impl ActiveCluster {
    fn new(record: &BedRecord, distance: u64) -> Result<Self> {
        Ok(Self {
            end: record.end(),
            reach: checked_reach(record, record.end(), distance)?,
        })
    }

    fn extend(&mut self, record: &BedRecord, distance: u64) -> Result<()> {
        if record.end() > self.end {
            self.end = record.end();
            self.reach = checked_reach(record, self.end, distance)?;
        }
        Ok(())
    }
}

fn checked_reach(record: &BedRecord, end: u64, distance: u64) -> Result<u64> {
    end.checked_add(distance).ok_or_else(|| {
        invalid(format!(
            "cluster distance overflows u64 at {}:{}-{} with distance {distance}",
            record.chrom(),
            record.start(),
            record.end()
        ))
    })
}

#[derive(Default)]
struct InputOrder {
    chromosome: Option<String>,
    last_start: u64,
    closed: HashSet<String>,
}

impl InputOrder {
    fn accept(&mut self, record: &BedRecord) -> Result<bool> {
        match self.chromosome.as_deref() {
            Some(chromosome) if chromosome == record.chrom() => {
                if record.start() < self.last_start {
                    return Err(invalid(format!(
                        "cluster input is not sorted: {}:{} follows start {}",
                        record.chrom(),
                        record.start(),
                        self.last_start
                    )));
                }
                self.last_start = record.start();
                Ok(false)
            }
            Some(chromosome) => {
                self.closed.insert(chromosome.to_owned());
                if self.closed.contains(record.chrom()) {
                    return Err(invalid(format!(
                        "cluster input chromosome {:?} reappears",
                        record.chrom()
                    )));
                }
                self.chromosome = Some(record.chrom().to_owned());
                self.last_start = record.start();
                Ok(true)
            }
            None => {
                self.chromosome = Some(record.chrom().to_owned());
                self.last_start = record.start();
                Ok(true)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(input: &[u8], distance: u64, same_strand: bool) -> Result<Vec<u8>> {
        let mut output = Vec::new();
        cluster(
            input,
            &mut output,
            ClusterOptions {
                distance,
                same_strand,
            },
        )?;
        Ok(output)
    }

    #[test]
    fn distance_joins_bookends_and_exact_gaps() {
        let output = run(
            b"chr1\t1\t10\nchr1\t10\t10\nchr1\t15\t20\nchr1\t26\t30\n",
            5,
            false,
        )
        .unwrap();
        assert_eq!(
            output,
            b"chr1\t1\t10\t1\nchr1\t10\t10\t1\nchr1\t15\t20\t1\nchr1\t26\t30\t2\n"
        );
    }

    #[test]
    fn same_strand_matches_chromosome_local_oracle_order() {
        let output = run(
            b"chr1\t1\t10\ta\t0\t-\n\
              chr1\t5\t15\tb\t0\t+\n\
              chr1\t8\t12\tc\t0\t-\n",
            0,
            true,
        )
        .unwrap();
        assert_eq!(
            output,
            b"chr1\t5\t15\tb\t0\t+\t1\n\
              chr1\t1\t10\ta\t0\t-\t2\n\
              chr1\t8\t12\tc\t0\t-\t2\n"
        );
    }

    #[test]
    fn start_order_and_chromosome_reappearance_fail_loud() {
        let unsorted = run(b"chr1\t10\t20\nchr1\t5\t6\n", 0, false).unwrap_err();
        assert!(unsorted.to_string().contains("not sorted"), "{unsorted}");

        let reappears = run(b"chr1\t1\t2\nchr2\t1\t2\nchr1\t3\t4\n", 0, false).unwrap_err();
        assert!(reappears.to_string().contains("reappears"), "{reappears}");
    }

    #[test]
    fn same_strand_requires_valid_bed6() {
        let missing = run(b"chr1\t1\t2\n", 0, true).unwrap_err();
        assert!(missing.to_string().contains("missing strand"), "{missing}");
        let invalid = run(b"chr1\t1\t2\ta\t0\t.\n", 0, true).unwrap_err();
        assert!(invalid.to_string().contains("invalid strand"), "{invalid}");
    }

    #[test]
    fn distance_reach_overflow_fails() {
        let input = format!("chr1\t0\t{}\n", u64::MAX);
        let error = run(input.as_bytes(), 1, false).unwrap_err();
        assert!(
            error.to_string().contains("distance overflows u64"),
            "{error}"
        );
    }
}
