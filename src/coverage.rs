//! Per-record coverage summaries and depth distributions over BED intervals.

use std::collections::BTreeMap;
use std::io::{BufWriter, Read, Write};

use rsomics_common::{Context, Result};

use crate::StrandFilter;
use crate::bed::{BedReader, BedRecord, Strand, invalid, virtual_bounds};
use crate::depth::DepthBuffer;
use crate::interval_index::{IntervalIndex, IntervalIndexBuilder};

/// Coverage columns appended to each original A record.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum CoverageMode {
    /// Hit count, covered bases, interval length and covered fraction.
    #[default]
    Summary,
    /// Number of qualifying overlaps.
    Counts,
    /// Mean depth across the entire A interval.
    Mean,
    /// Observed depth bins followed by aggregate bins over all A intervals.
    Histogram,
    /// One-based relative position and depth for every base of A.
    Depth,
}

/// Selection and reporting options for [`coverage`].
#[derive(Debug, Clone, Copy, Default)]
pub struct CoverageOptions {
    /// Columns to report.
    pub mode: CoverageMode,
    /// Strand relationship required between A and B.
    pub strand: StrandFilter,
    /// Use BED12 blocks instead of contiguous record spans.
    pub split: bool,
    /// Minimum fraction of A overlapped by each B candidate, in (0, 1].
    pub minimum_a: Option<f32>,
    /// Minimum fraction of B overlapped by A, in (0, 1].
    pub minimum_b: Option<f32>,
    /// Apply the A fraction to both A and B.
    pub reciprocal: bool,
    /// Accept either fraction requirement instead of requiring both.
    pub either: bool,
}

impl CoverageOptions {
    fn validate(self) -> Result<()> {
        for value in [self.minimum_a, self.minimum_b].into_iter().flatten() {
            if !value.is_finite() || value <= 0.0 || value > 1.0 {
                return Err(invalid(
                    "coverage overlap fractions must be finite and in (0, 1]",
                ));
            }
        }
        if self.reciprocal && self.minimum_b.is_some() {
            return Err(invalid(
                "coverage reciprocal selection cannot specify a B fraction",
            ));
        }
        if self.either && self.minimum_b.is_some() && self.minimum_a.is_none() {
            return Err(invalid(
                "coverage either selection with a B fraction also requires an A fraction",
            ));
        }
        if self.split
            && (self.minimum_a.is_some()
                || self.minimum_b.is_some()
                || self.reciprocal
                || self.either)
        {
            return Err(invalid(
                "coverage split with fraction selection is not supported",
            ));
        }
        Ok(())
    }

    #[allow(clippy::cast_precision_loss)]
    fn accepts(self, overlap: u64, a_len: u64, b_len: u64) -> bool {
        let a = self.minimum_a.unwrap_or(1e-9);
        let b = if self.reciprocal {
            a
        } else {
            self.minimum_b.unwrap_or(1e-9)
        };
        let sufficient_a = overlap as f32 / a_len as f32 >= a;
        let sufficient_b = overlap as f32 / b_len as f32 >= b;
        if self.either {
            sufficient_a || sufficient_b
        } else {
            sufficient_a && sufficient_b
        }
    }
}

struct CoverageIndex {
    index: IntervalIndex,
    strands: Vec<Strand>,
    blocks: Vec<Vec<(u64, u64)>>,
}

impl CoverageIndex {
    fn load<R: Read>(
        inputs: impl IntoIterator<Item = R>,
        options: CoverageOptions,
    ) -> Result<Self> {
        let mut builder = IntervalIndexBuilder::new();
        let mut strands = Vec::new();
        let mut blocks = Vec::new();
        for (number, input) in inputs.into_iter().enumerate() {
            let label = format!("coverage B{}", number + 1);
            let mut reader = BedReader::new(input);
            while let Some(record) = reader.next_record()? {
                builder.push(record.chrom(), record.start(), record.end(), &label)?;
                if options.strand != StrandFilter::Any {
                    strands.push(record.strand(&label)?);
                }
                if options.split {
                    blocks.push(record.blocks(&label)?);
                }
            }
        }
        Ok(Self {
            index: builder.finish("coverage B")?,
            strands,
            blocks,
        })
    }
}

/// Pool B inputs and append coverage measurements to each A record in input order.
///
/// BED12 split overlap counts refer to overlapping block pairs. Fractions and
/// means use binary32 division followed by seven decimal places. Unknown or
/// missing strands are errors when a strand filter is selected.
///
/// # Errors
///
/// Returns errors for malformed selected BED fields, unsupported split/fraction
/// selection, invalid fractions, coordinate/index limits, or I/O failure.
pub fn coverage<R: Read>(
    a_input: impl Read,
    b_inputs: impl IntoIterator<Item = R>,
    output: impl Write,
    options: CoverageOptions,
) -> Result<()> {
    options.validate()?;
    let b = CoverageIndex::load(b_inputs, options)?;
    let mut reader = BedReader::new(a_input);
    let mut output = BufWriter::new(output);
    let mut candidates = Vec::new();
    let mut depths = DepthBuffer::default();
    let mut aggregate = BTreeMap::<usize, u128>::new();
    let mut total_length = 0_u128;

    while let Some(record) = reader.next_record()? {
        let (start, end) = virtual_bounds(&record, "coverage A")?;
        let length = end - start;
        b.index.intersection_candidates(
            record.chrom(),
            record.start(),
            record.end(),
            "coverage A",
            &mut candidates,
        )?;
        let strand = if options.strand == StrandFilter::Any {
            None
        } else {
            Some(record.strand("coverage A")?)
        };
        let blocks = if options.split {
            record.blocks("coverage A")?
        } else {
            Vec::new()
        };
        depths.clear_events();
        let mut count = 0_usize;
        for &id in &candidates {
            let eligible = match options.strand {
                StrandFilter::Any => true,
                StrandFilter::Same => strand == Some(b.strands[id]),
                StrandFilter::Opposite => strand != Some(b.strands[id]),
            };
            if !eligible {
                continue;
            }
            if options.split {
                for &(a_start, a_end) in &blocks {
                    for &(b_start, b_end) in &b.blocks[id] {
                        let low = a_start.max(b_start);
                        let high = a_end.min(b_end);
                        if low < high {
                            count += 1;
                            if options.mode != CoverageMode::Counts {
                                depths.add(low, high);
                            }
                        }
                    }
                }
            } else {
                let (b_start, b_end) = b.index.record(id).virtual_bounds();
                let low = start.max(b_start);
                let high = end.min(b_end);
                if low < high && options.accepts(high - low, length, b_end - b_start) {
                    count += 1;
                    if options.mode != CoverageMode::Counts {
                        depths.add(low, high);
                    }
                }
            }
        }
        if options.mode == CoverageMode::Counts {
            record.write_column(&mut output, count)?;
            continue;
        }
        depths.reduce(start, end);
        report(
            &record,
            &depths,
            count,
            start,
            length,
            &mut output,
            options.mode,
        )?;
        if options.mode == CoverageMode::Histogram {
            total_length += u128::from(length);
            for (&depth, &bases) in depths.histogram() {
                *aggregate.entry(depth).or_default() += u128::from(bases);
            }
        }
    }
    for (depth, bases) in aggregate {
        writeln!(
            output,
            "all\t{depth}\t{bases}\t{total_length}\t{:.7}",
            ratio(bases, total_length)
        )
        .rs_context("writing coverage histogram totals")?;
    }
    output.flush().rs_context("flushing coverage output")
}

#[allow(clippy::cast_precision_loss)]
fn ratio(numerator: u128, denominator: u128) -> f32 {
    numerator as f32 / denominator as f32
}

fn report(
    record: &BedRecord,
    depths: &DepthBuffer,
    count: usize,
    start: u64,
    length: u64,
    output: &mut dyn Write,
    mode: CoverageMode,
) -> Result<()> {
    match mode {
        CoverageMode::Summary => {
            let covered = depths.covered_bases(length);
            record.write_column(
                output,
                format_args!(
                    "{count}\t{covered}\t{length}\t{:.7}",
                    ratio(u128::from(covered), u128::from(length))
                ),
            )?;
        }
        CoverageMode::Mean => {
            let sum: u128 = depths
                .histogram()
                .iter()
                .map(|(&depth, &bases)| depth as u128 * u128::from(bases))
                .sum();
            record.write_column(
                output,
                format_args!("{:.7}", ratio(sum, u128::from(length))),
            )?;
        }
        CoverageMode::Histogram => {
            for (&depth, &bases) in depths.histogram() {
                record.write_column(
                    output,
                    format_args!(
                        "{depth}\t{bases}\t{length}\t{:.7}",
                        ratio(u128::from(bases), u128::from(length))
                    ),
                )?;
            }
        }
        CoverageMode::Depth => {
            for run in depths.runs() {
                for position in run.start..run.end {
                    record.write_column(
                        output,
                        format_args!("{}\t{}", position - start + 1, run.depth),
                    )?;
                }
            }
        }
        CoverageMode::Counts => unreachable!("count output bypasses depth reduction"),
    }
    Ok(())
}
