//! Per-file overlap counts and covered fractions in BEDTools annotation order.

use std::fmt::{self, Display};
use std::io::{BufWriter, Read, Write};

use rsomics_common::{Context, Result};

use crate::StrandFilter;
use crate::bed::{BedReader, BedRecord, Strand, invalid, virtual_bounds};
use crate::depth::DepthBuffer;
use crate::interval_index::{IntervalIndex, IntervalIndexBuilder};

/// Columns emitted for each annotation file.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum AnnotateMode {
    /// Fraction of the query interval covered by annotation records.
    #[default]
    Fraction,
    /// Number of overlapping annotation records, including duplicates.
    Counts,
    /// Count followed by covered fraction.
    Both,
}

/// Reporting, strand selection and optional annotation-file labels.
#[derive(Debug, Clone, Default)]
pub struct AnnotateOptions {
    /// Columns emitted per annotation file.
    pub mode: AnnotateMode,
    /// Strand relationship required between query and annotation records.
    pub strand: StrandFilter,
    /// One label per annotation file, or empty to omit the generated header.
    pub names: Vec<String>,
}

struct Query {
    record: BedRecord,
    start: u64,
    end: u64,
    bin: u64,
    strand: Option<Strand>,
}

struct AnnotationIndex {
    index: IntervalIndex,
    strands: Vec<Strand>,
}

impl AnnotationIndex {
    fn load(input: impl Read, label: &str, strand: StrandFilter) -> Result<Self> {
        let mut reader = BedReader::new(input);
        let mut builder = IntervalIndexBuilder::new();
        let mut strands = Vec::new();
        let mut width = None;
        while let Some(record) = reader.next_record()? {
            check_width(&record, &mut width, label)?;
            builder.push(record.chrom(), record.start(), record.end(), label)?;
            if strand != StrandFilter::Any {
                strands.push(record.strand(label)?);
            }
        }
        Ok(Self {
            index: builder.finish(label)?,
            strands,
        })
    }
}

#[derive(Clone, Copy, Default)]
struct Metric {
    count: usize,
    covered: u64,
}

struct Columns<'a> {
    metrics: &'a [Metric],
    length: u64,
    mode: AnnotateMode,
}

impl Display for Columns<'_> {
    #[allow(clippy::cast_precision_loss)]
    fn fmt(&self, output: &mut fmt::Formatter<'_>) -> fmt::Result {
        for (index, metric) in self.metrics.iter().enumerate() {
            if index != 0 {
                output.write_str("\t")?;
            }
            match self.mode {
                AnnotateMode::Counts => write!(output, "{}", metric.count)?,
                AnnotateMode::Fraction => {
                    write!(output, "{:.6}", metric.covered as f64 / self.length as f64)?
                }
                AnnotateMode::Both => write!(
                    output,
                    "{}\t{:.6}",
                    metric.count,
                    metric.covered as f64 / self.length as f64
                )?,
            }
        }
        Ok(())
    }
}

/// Append independent annotation-file results in chromosome/bin/input-tie order.
///
/// Fractions use binary64 division and six decimal places. Query points use a
/// checked two-base virtual span; annotation points count as overlaps but add
/// no covered bases. BED12 is interpreted by its outer span. Names produce a
/// generated header even when the query input is empty.
///
/// # Errors
///
/// Returns errors for missing annotation inputs, mismatched names, malformed or
/// inconsistent BED rows, invalid selected strands, coordinate/index limits,
/// or I/O failure.
pub fn annotate<R: Read>(
    input: impl Read,
    annotations: impl IntoIterator<Item = R>,
    output: impl Write,
    options: AnnotateOptions,
) -> Result<()> {
    let annotations: Vec<R> = annotations.into_iter().collect();
    let file_count = annotations.len();
    if file_count == 0 {
        return Err(invalid("annotate requires at least one annotation file"));
    }
    if !options.names.is_empty() && options.names.len() != file_count {
        return Err(invalid(
            "annotate names must match the annotation-file count",
        ));
    }
    let mut reader = BedReader::new(input);
    let mut queries = Vec::new();
    let mut width = None;
    while let Some(record) = reader.next_record()? {
        check_width(&record, &mut width, "annotate A")?;
        let (start, end) = virtual_bounds(&record, "annotate A")?;
        let strand = if options.strand == StrandFilter::Any {
            None
        } else {
            Some(record.strand("annotate A")?)
        };
        queries.push(Query {
            bin: source_bin(start, end)?,
            record,
            start,
            end,
            strand,
        });
    }
    let size = queries
        .len()
        .checked_mul(file_count)
        .ok_or_else(|| invalid("annotate result dimensions exceed the addressable domain"))?;
    let mut metrics = vec![Metric::default(); size];
    let mut candidates = Vec::new();
    let mut depths = DepthBuffer::default();
    for (file, input) in annotations.into_iter().enumerate() {
        let label = format!("annotate B{}", file + 1);
        let annotation = AnnotationIndex::load(input, &label, options.strand)?;
        for (id, query) in queries.iter().enumerate() {
            annotation.index.intersection_candidates(
                query.record.chrom(),
                query.record.start(),
                query.record.end(),
                "annotate A",
                &mut candidates,
            )?;
            depths.clear_events();
            let metric = &mut metrics[id * file_count + file];
            for &candidate in &candidates {
                let eligible = match options.strand {
                    StrandFilter::Any => true,
                    StrandFilter::Same => query.strand == Some(annotation.strands[candidate]),
                    StrandFilter::Opposite => query.strand != Some(annotation.strands[candidate]),
                };
                if !eligible {
                    continue;
                }
                metric.count += 1;
                let interval = annotation.index.record(candidate);
                if options.mode != AnnotateMode::Counts && interval.start != interval.end {
                    depths.add(query.start.max(interval.start), query.end.min(interval.end));
                }
            }
            if options.mode != AnnotateMode::Counts {
                depths.reduce(query.start, query.end);
                metric.covered = depths.covered_bases(query.end - query.start);
            }
        }
    }
    let mut order: Vec<usize> = (0..queries.len()).collect();
    order.sort_unstable_by(|&left, &right| {
        queries[left]
            .record
            .chrom()
            .cmp(queries[right].record.chrom())
            .then_with(|| queries[left].bin.cmp(&queries[right].bin))
            .then_with(|| left.cmp(&right))
    });
    let mut output = BufWriter::new(output);
    write_header(&mut output, width.unwrap_or(0), &options)?;
    for id in order {
        let query = &queries[id];
        query.record.write_column(
            &mut output,
            Columns {
                metrics: &metrics[id * file_count..(id + 1) * file_count],
                length: query.end - query.start,
                mode: options.mode,
            },
        )?;
    }
    output.flush().rs_context("flushing annotate output")
}

fn check_width(record: &BedRecord, width: &mut Option<usize>, label: &str) -> Result<()> {
    let count = record.field_count();
    if let Some(expected) = width {
        if *expected != count {
            return Err(invalid(format!(
                "{label} interval {}:{}-{} has {count} fields; expected {expected}",
                record.chrom(),
                record.start(),
                record.end()
            )));
        }
    } else {
        *width = Some(count);
    }
    Ok(())
}

fn source_bin(start: u64, end: u64) -> Result<u64> {
    let mut start = start >> 14;
    let mut end = (end - 1) >> 14;
    for offset in [299503, 37359, 4681, 585, 73, 9, 1, 0] {
        if start == end {
            return Ok(offset + start);
        }
        start >>= 3;
        end >>= 3;
    }
    Err(invalid(
        "annotate interval exceeds the source binning domain",
    ))
}

fn write_header(output: &mut dyn Write, width: usize, options: &AnnotateOptions) -> Result<()> {
    if options.names.is_empty() {
        return Ok(());
    }
    output
        .write_all(b"#")
        .rs_context("writing annotation header")?;
    for _ in 1..width {
        output
            .write_all(b"\t")
            .rs_context("writing annotation header")?;
    }
    for name in &options.names {
        if options.mode == AnnotateMode::Both {
            write!(output, "\t{name}_cnt\t{name}_pct").rs_context("writing annotation header")?;
        } else {
            write!(output, "\t{name}").rs_context("writing annotation header")?;
        }
    }
    output
        .write_all(b"\n")
        .rs_context("writing annotation header")
}
