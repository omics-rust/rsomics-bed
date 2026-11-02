//! Sorted interval mapping with per-column reductions.

mod reduce;
mod sweep;

use std::io::{BufWriter, Read, Write};
use std::num::NonZeroUsize;
use std::str::FromStr;

use crate::bed::{BedReader, invalid};
use crate::{Genome, StrandFilter};
use reduce::{Format, Reducer};
use rsomics_common::{Context, Result, RsomicsError};

/// Reduction applied to selected values from overlapping B records.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Operation {
    /// Sum numeric values.
    Sum,
    /// Arithmetic mean.
    Mean,
    /// Population standard deviation.
    Stdev,
    /// Sample standard deviation; a singleton produces the null token.
    SampleStdev,
    /// Median numeric value.
    Median,
    /// Most frequent byte value, with lexical ties.
    Mode,
    /// Least frequent byte value, with lexical ties.
    Antimode,
    /// Minimum numeric value.
    Min,
    /// Maximum numeric value.
    Max,
    /// Minimum absolute numeric value.
    AbsMin,
    /// Maximum absolute numeric value.
    AbsMax,
    /// Number of overlapping records.
    Count,
    /// Lexically sorted unique byte values.
    Distinct,
    /// Number of unique byte values.
    CountDistinct,
    /// Values occurring once, retaining the source leading-delimiter convention.
    DistinctOnly,
    /// Unique numeric values in ascending order.
    DistinctSortNum,
    /// Unique numeric values in descending order.
    DistinctSortNumDesc,
    /// Delimited values in hit order, retaining duplicates.
    Collapse,
    /// Concatenated values in hit order.
    Concat,
    /// Value frequencies in ascending frequency order.
    FreqAsc,
    /// Value frequencies in descending frequency order.
    FreqDesc,
    /// First value in hit order.
    First,
    /// Last value in hit order.
    Last,
}

impl FromStr for Operation {
    type Err = RsomicsError;

    fn from_str(value: &str) -> Result<Self> {
        Ok(match value {
            "sum" => Self::Sum,
            "mean" => Self::Mean,
            "stdev" => Self::Stdev,
            "sstdev" => Self::SampleStdev,
            "median" => Self::Median,
            "mode" => Self::Mode,
            "antimode" => Self::Antimode,
            "min" => Self::Min,
            "max" => Self::Max,
            "absmin" => Self::AbsMin,
            "absmax" => Self::AbsMax,
            "count" => Self::Count,
            "distinct" => Self::Distinct,
            "count_distinct" => Self::CountDistinct,
            "distinct_only" => Self::DistinctOnly,
            "distinct_sort_num" => Self::DistinctSortNum,
            "distinct_sort_num_desc" => Self::DistinctSortNumDesc,
            "collapse" => Self::Collapse,
            "concat" => Self::Concat,
            "freqasc" => Self::FreqAsc,
            "freqdesc" => Self::FreqDesc,
            "first" => Self::First,
            "last" => Self::Last,
            _ => return Err(invalid(format!("unknown map operation {value:?}"))),
        })
    }
}

/// Column, overlap and reporting options for [`map`].
#[derive(Debug, Clone)]
pub struct MapOptions {
    /// One-based B columns; a single column broadcasts over operations.
    pub columns: Vec<usize>,
    /// Reductions; a single operation broadcasts over columns.
    pub operations: Vec<Operation>,
    /// Significant digits for scalar numeric output; numeric distinct uses six.
    pub precision: NonZeroUsize,
    /// Separator for value-list reductions.
    pub delimiter: Vec<u8>,
    /// No-value token, except counts which produce zero.
    pub null: Vec<u8>,
    /// Required valid BED6 strand relationship.
    pub strand: StrandFilter,
    /// Match BED12 blocks, counting each B record once.
    pub split: bool,
    /// Minimum fraction of A overlapped by B, in (0, 1].
    pub minimum_a: Option<f32>,
    /// Minimum fraction of B overlapped by A, in (0, 1].
    pub minimum_b: Option<f32>,
    /// Copy an explicitly supplied A fraction to B.
    pub reciprocal: bool,
    /// Use OR rather than AND for the two fraction requirements.
    pub either: bool,
    /// Sort hits by effective coordinates, retaining B order for equal ties.
    pub sort_hits: bool,
    /// Preserve the initial A header before data output.
    pub header: bool,
}

impl Default for MapOptions {
    fn default() -> Self {
        Self {
            columns: vec![5],
            operations: vec![Operation::Sum],
            precision: NonZeroUsize::new(10).expect("ten is positive"),
            delimiter: b",".to_vec(),
            null: b".".to_vec(),
            strand: StrandFilter::Any,
            split: false,
            minimum_a: None,
            minimum_b: None,
            reciprocal: false,
            either: false,
            sort_hits: false,
            header: false,
        }
    }
}

impl MapOptions {
    fn validate(&self) -> Result<()> {
        if self.columns.is_empty() || self.operations.is_empty() || self.columns.contains(&0) {
            return Err(invalid(
                "map requires positive columns and at least one operation",
            ));
        }
        if self.columns.len() != 1
            && self.operations.len() != 1
            && self.columns.len() != self.operations.len()
        {
            return Err(invalid(
                "map columns and operations must have equal lengths or one must broadcast",
            ));
        }
        for value in [self.minimum_a, self.minimum_b].into_iter().flatten() {
            if !value.is_finite() || value <= 0.0 || value > 1.0 {
                return Err(invalid(
                    "map overlap fractions must be finite and in (0, 1]",
                ));
            }
        }
        if self.reciprocal && self.minimum_b.is_some() {
            return Err(invalid(
                "map reciprocal selection cannot specify a B fraction",
            ));
        }
        if self.either && self.minimum_b.is_some() && self.minimum_a.is_none() {
            return Err(invalid(
                "map either selection with a B fraction also requires an A fraction",
            ));
        }
        if self.split
            && (self.minimum_a.is_some()
                || self.minimum_b.is_some()
                || self.reciprocal
                || self.either)
        {
            return Err(invalid(
                "map split with explicit fraction selection is not supported",
            ));
        }
        Ok(())
    }
}

/// Append B-column reductions to each sorted A record.
///
/// Both inputs must contain grouped chromosomes and nondecreasing raw starts.
/// A genome file supplies explicit chromosome order, not coordinate clipping.
/// Origin points use effective [-1, 1) bounds without changing output coordinates.
/// Equal-start point candidates are retained regardless of A/B tie ordering.
///
/// # Errors
///
/// Returns errors for incompatible options, empty B, missing columns, inconsistent
/// record widths, sortedness violations, invalid selected values or I/O failure.
pub fn map(
    a: impl Read,
    b: impl Read,
    output: impl Write,
    options: &MapOptions,
    genome: Option<&Genome>,
) -> Result<()> {
    options.validate()?;
    let count = options.columns.len().max(options.operations.len());
    let mut columns = Vec::new();
    let mut reductions = Vec::with_capacity(count);
    for index in 0..count {
        let column = options.columns[if options.columns.len() == 1 { 0 } else { index }] - 1;
        let group = columns
            .iter()
            .position(|&old| old == column)
            .unwrap_or_else(|| {
                columns.push(column);
                columns.len() - 1
            });
        let operation = options.operations[if options.operations.len() == 1 {
            0
        } else {
            index
        }];
        reductions.push((group, operation));
    }
    let mut queries = if options.header {
        BedReader::new(a).preserve_header()
    } else {
        BedReader::new(a)
    };
    let mut current = queries.next_record()?;
    let mut sweep = sweep::Sweep::new(b, &columns, options, genome)?;
    let mut output = BufWriter::new(output);
    if options.header {
        output
            .write_all(queries.header())
            .rs_context("writing map header")?;
    }
    let mut hits = Vec::new();
    while let Some(record) = current {
        sweep.collect_hits(&record, &mut hits)?;
        let values: Vec<Vec<&[u8]>> = columns
            .iter()
            .map(|&column| {
                hits.iter()
                    .map(|&hit| sweep.record(hit).column(column, "map B"))
                    .collect()
            })
            .collect::<Result<_>>()?;
        let mut reducers: Vec<_> = values
            .iter()
            .map(|values| {
                Reducer::new(
                    values,
                    Format {
                        precision: options.precision,
                        delimiter: &options.delimiter,
                        null: &options.null,
                    },
                )
            })
            .collect();
        record.write_prefix(&mut output, record.start(), record.end())?;
        for &(group, operation) in &reductions {
            output
                .write_all(b"\t")
                .rs_context("writing map separator")?;
            reducers[group].write(operation, &mut output)?;
        }
        output.write_all(b"\n").rs_context("writing map record")?;
        current = queries.next_record()?;
    }
    sweep.finish()?;
    output.flush().rs_context("flushing map output")
}
