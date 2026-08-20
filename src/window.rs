//! Indexed BED neighborhood queries with typed reporting and strand policies.

use std::io::{BufWriter, Read, Write};

use rsomics_common::{Context, Result};

use crate::bed::{BedReader, Strand, invalid, virtual_bounds};
use crate::relation_index::RelationBed;

pub use crate::StrandFilter;

/// Output reduction applied to each A record.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum WindowReport {
    /// Emit one joined row per matching B record.
    #[default]
    Pairs,
    /// Emit A once when at least one B record matches.
    Any,
    /// Append the number of matching B records to A.
    Count,
    /// Emit A only when no B record matches.
    None,
}

/// Options for [`window`].
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct WindowOptions {
    /// Bases added to the reference-left side of A.
    pub left: u64,
    /// Bases added to the reference-right side of A.
    pub right: u64,
    /// Interpret left and right relative to A's strand.
    pub strand_relative: bool,
    /// Strand relationship required between A and B.
    pub strand: StrandFilter,
    /// Output reduction.
    pub report: WindowReport,
}

impl Default for WindowOptions {
    fn default() -> Self {
        Self {
            left: 1000,
            right: 1000,
            strand_relative: false,
            strand: StrandFilter::Any,
            report: WindowReport::Pairs,
        }
    }
}

/// Report B records within the configured neighborhood of each A record.
///
/// # Errors
///
/// Returns an error for malformed BED, inconsistent B field widths, missing or
/// invalid required strands, coordinate overflow, backend limits, or I/O
/// failure.
pub fn window(
    a_input: impl Read,
    b_input: impl Read,
    output: impl Write,
    options: WindowOptions,
) -> Result<()> {
    let b = RelationBed::load(b_input, "window B")?;
    let b_strands = if options.strand == StrandFilter::Any {
        None
    } else {
        Some(b.checked_strands("window B")?)
    };
    let mut a = BedReader::new(a_input);
    let mut output = BufWriter::new(output);
    let mut candidate_ids = Vec::new();
    let mut eligible_ids = Vec::new();

    while let Some(record) = a.next_record()? {
        let a_strand = if options.strand_relative || options.strand != StrandFilter::Any {
            Some(record.strand("window A")?)
        } else {
            None
        };
        let (left, right) = if options.strand_relative && a_strand == Some(Strand::Reverse) {
            (options.right, options.left)
        } else {
            (options.left, options.right)
        };
        let (start, end) = virtual_bounds(&record, "window")?;
        let start = start.saturating_sub(left);
        let end = end.checked_add(right).ok_or_else(|| {
            invalid(format!(
                "window overflows u64 at {}:{}-{} with right distance {right}",
                record.chrom(),
                record.start(),
                record.end()
            ))
        })?;
        b.range_candidates(
            record.chrom(),
            start,
            end,
            "window A interval index",
            &mut candidate_ids,
        )?;

        eligible_ids.clear();
        for &id in &candidate_ids {
            let eligible = match options.strand {
                StrandFilter::Any => true,
                StrandFilter::Same => {
                    b_strands.as_ref().expect("strand filter is loaded")[id]
                        == a_strand.expect("strand filter requires A strand")
                }
                StrandFilter::Opposite => {
                    b_strands.as_ref().expect("strand filter is loaded")[id]
                        != a_strand.expect("strand filter requires A strand")
                }
            };
            if eligible {
                eligible_ids.push(id);
            }
        }

        match options.report {
            WindowReport::Pairs => {
                for &id in &eligible_ids {
                    record.write_joined_raw(&mut output, b.record(id).raw())?;
                }
            }
            WindowReport::Any if !eligible_ids.is_empty() => record.write_raw(&mut output)?,
            WindowReport::Count => record.write_column(&mut output, eligible_ids.len())?,
            WindowReport::None if eligible_ids.is_empty() => record.write_raw(&mut output)?,
            WindowReport::Any | WindowReport::None => {}
        }
    }
    output.flush().rs_context("flushing window BED output")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn run(a: &[u8], b: &[u8], options: WindowOptions) -> rsomics_common::Result<Vec<u8>> {
        let mut output = Vec::new();
        window(a, b, &mut output, options)?;
        Ok(output)
    }

    #[test]
    fn report_modes_reduce_file_order_candidates() {
        let a = b"chr1\t10\t20\tA\nchr2\t10\t20\tN\n";
        let b = b"chr1\t15\t25\tfirst\nchr1\t5\t15\tsecond\n";
        let base = WindowOptions {
            left: 0,
            right: 0,
            ..WindowOptions::default()
        };

        assert_eq!(
            run(a, b, base).unwrap(),
            b"chr1\t10\t20\tA\tchr1\t15\t25\tfirst\n\
              chr1\t10\t20\tA\tchr1\t5\t15\tsecond\n"
        );
        assert_eq!(
            run(
                a,
                b,
                WindowOptions {
                    report: WindowReport::Count,
                    ..base
                }
            )
            .unwrap(),
            b"chr1\t10\t20\tA\t2\nchr2\t10\t20\tN\t0\n"
        );
        assert_eq!(
            run(
                a,
                b,
                WindowOptions {
                    report: WindowReport::Any,
                    ..base
                }
            )
            .unwrap(),
            b"chr1\t10\t20\tA\n"
        );
        assert_eq!(
            run(
                a,
                b,
                WindowOptions {
                    report: WindowReport::None,
                    ..base
                }
            )
            .unwrap(),
            b"chr2\t10\t20\tN\n"
        );
    }

    #[test]
    fn strand_relative_swaps_asymmetric_sides() {
        let a = b"chr1\t50\t60\tA\t0\t-\n";
        let b = b"chr1\t65\t70\tright\t0\t-\nchr1\t40\t45\tleft\t0\t-\n";
        let output = run(
            a,
            b,
            WindowOptions {
                left: 10,
                right: 2,
                strand_relative: true,
                strand: StrandFilter::Same,
                report: WindowReport::Pairs,
            },
        )
        .unwrap();
        assert_eq!(
            output,
            b"chr1\t50\t60\tA\t0\t-\tchr1\t65\t70\tright\t0\t-\n"
        );
    }

    #[test]
    fn selected_strand_modes_require_bed6() {
        let missing_a = run(
            b"chr1\t10\t20\n",
            b"chr1\t10\t20\tB\t0\t+\n",
            WindowOptions {
                strand: StrandFilter::Same,
                ..WindowOptions::default()
            },
        )
        .unwrap_err();
        assert!(missing_a.to_string().contains("window A"), "{missing_a}");

        let invalid_b = run(
            b"chr1\t10\t20\tA\t0\t+\n",
            b"chr2\t10\t20\tB\t0\t.\n",
            WindowOptions {
                strand: StrandFilter::Same,
                ..WindowOptions::default()
            },
        )
        .unwrap_err();
        assert!(invalid_b.to_string().contains("window B"), "{invalid_b}");
    }

    #[test]
    fn expansion_is_zero_bounded_and_overflow_checked() {
        let output = run(
            b"chr1\t1\t2\tA\n",
            b"chr1\t0\t1\tB\n",
            WindowOptions {
                left: 10,
                right: 0,
                ..WindowOptions::default()
            },
        )
        .unwrap();
        assert_eq!(output, b"chr1\t1\t2\tA\tchr1\t0\t1\tB\n");

        let input = format!("chr1\t{}\t{}\tA\n", u64::MAX - 1, u64::MAX);
        let error = run(
            input.as_bytes(),
            b"",
            WindowOptions {
                left: 0,
                right: 1,
                ..WindowOptions::default()
            },
        )
        .unwrap_err();
        assert!(
            error.to_string().contains("window overflows u64"),
            "{error}"
        );
    }
}
