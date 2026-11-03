//! Streaming chromosome-bounded BED coordinate transforms.

use std::io::{BufWriter, Read, Write};

use rsomics_common::{Context, Result, RsomicsError};

use crate::bed::{BedReader, BedRecord, Genome};

/// Resizing amounts and output policy for slop.
#[derive(Debug, Clone, Copy, Default)]
pub struct SlopOptions {
    /// Amount subtracted from the start; negative values shrink.
    pub left: f32,
    /// Amount added to the end; negative values shrink.
    pub right: f32,
    /// Interpret amounts as fractions of the input feature length.
    pub fraction: bool,
    /// Exchange left/right amounts on the negative strand.
    pub strand_relative: bool,
    /// Preserve initial input headers.
    pub header: bool,
}

/// Strand-selected displacement and output policy for shift.
#[derive(Debug, Clone, Copy, Default)]
pub struct ShiftOptions {
    /// Displacement for all records not marked with a negative strand.
    pub plus: f32,
    /// Displacement for negative-strand records.
    pub minus: f32,
    /// Interpret displacements as fractions of the input feature length.
    pub fraction: bool,
    /// Preserve initial input headers.
    pub header: bool,
}

/// Resize records, clipping to chromosome bounds and preserving trailing fields.
///
/// Negative amounts shrink; crossed endpoints are exchanged and collapsed
/// spans are expanded to one base before restoring point-record coordinates.
///
/// # Errors
/// Returns an error for malformed input, nonfinite amounts, missing or empty
/// chromosomes, unrepresentable coordinates, or input/output failures.
pub fn slop(
    input: impl Read,
    genome: &Genome,
    output: impl Write,
    options: &SlopOptions,
) -> Result<()> {
    finite(options.left, options.right)?;
    stream(input, genome, output, options.header, |record, size| {
        let (start, end, point) = bounds(record);
        let (left, right) = if options.strand_relative && reverse(record) {
            (options.right, options.left)
        } else {
            (options.left, options.right)
        };
        let size = i128::from(size);
        let mut low = (start - amount(left, end - start, options.fraction)?).clamp(0, size);
        let mut high = (end + amount(right, end - start, options.fraction)?).clamp(0, size);
        if low == high {
            if low == size {
                low -= 1;
            } else {
                high += 1;
            }
        } else if low > high {
            std::mem::swap(&mut low, &mut high);
        }
        reported(low, high, point)
    })
}

/// Move records with independent endpoint clipping and unchanged trailing fields.
///
/// Negative strands use `minus`; missing or other strands use `plus`.
///
/// # Errors
/// Returns an error for malformed input, nonfinite displacements, missing or
/// empty chromosomes, unrepresentable coordinates, or input/output failures.
pub fn shift(
    input: impl Read,
    genome: &Genome,
    output: impl Write,
    options: &ShiftOptions,
) -> Result<()> {
    finite(options.plus, options.minus)?;
    stream(input, genome, output, options.header, |record, size| {
        let (start, end, point) = bounds(record);
        let mut offset = f64::from(if reverse(record) {
            options.minus
        } else {
            options.plus
        });
        if options.fraction {
            offset *= (end - start) as f64;
        }
        let low = start as f64 + offset;
        let high = end as f64 + offset;
        let low = if low < 0.0 {
            0
        } else if low > (size - 1) as f64 {
            i128::from(size - 1)
        } else {
            low as i128
        };
        let high = if high <= 0.0 {
            1
        } else if high > size as f64 {
            i128::from(size)
        } else {
            high as i128
        };
        reported(low, high, point)
    })
}

fn stream(
    input: impl Read,
    genome: &Genome,
    output: impl Write,
    header: bool,
    transform: impl Fn(&BedRecord, u64) -> Result<(u64, u64)>,
) -> Result<()> {
    let mut records = if header {
        BedReader::new(input).preserve_header()
    } else {
        BedReader::new(input)
    };
    let mut current = records.next_record()?;
    let mut output = BufWriter::with_capacity(65_536, output);
    if header {
        output
            .write_all(records.header())
            .rs_context("writing BED header")?;
    }
    while let Some(record) = current {
        let size = genome.size(record.chrom()).ok_or_else(|| {
            invalid(format!(
                "chromosome {:?} is absent from the genome file",
                record.chrom()
            ))
        })?;
        if size == 0 {
            return Err(invalid(format!(
                "chromosome {:?} has zero size",
                record.chrom()
            )));
        }
        let (start, end) = transform(&record, size)?;
        if end > size {
            return Err(invalid(format!(
                "coordinate transformation exceeds chromosome {:?} size {size}",
                record.chrom()
            )));
        }
        record.write_with_coords(&mut output, start, end)?;
        current = records.next_record()?;
    }
    output.flush().rs_context("flushing transformed BED")
}

fn bounds(record: &BedRecord) -> (i128, i128, bool) {
    let start = i128::from(record.start());
    let end = i128::from(record.end());
    if start == end {
        (start - 1, end + 1, true)
    } else {
        (start, end, false)
    }
}

fn reverse(record: &BedRecord) -> bool {
    record.optional_column(5) == Some(b"-".as_slice())
}

fn amount(value: f32, length: i128, fraction: bool) -> Result<i128> {
    let value = if fraction {
        value * length as f32
    } else {
        value
    };
    let limit = -(i64::MIN as f32);
    if !value.is_finite() || !(-limit..limit).contains(&value) {
        return Err(invalid(
            "resizing amount exceeds the supported coordinate range",
        ));
    }
    Ok(value as i128)
}

fn reported(mut start: i128, mut end: i128, point: bool) -> Result<(u64, u64)> {
    if point {
        start += 1;
        end -= 1;
    }
    match (u64::try_from(start), u64::try_from(end)) {
        (Ok(start), Ok(end)) if start <= end => Ok((start, end)),
        _ => Err(invalid(
            "coordinate transformation does not yield a valid BED interval",
        )),
    }
}

fn finite(left: f32, right: f32) -> Result<()> {
    if left.is_finite() && right.is_finite() {
        Ok(())
    } else {
        Err(invalid("coordinate amounts must be finite"))
    }
}

fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}
