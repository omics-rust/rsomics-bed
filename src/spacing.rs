//! Previous-record spacing for chromosome-grouped BED streams.

use std::io::{BufWriter, Read, Write};

use rsomics_common::{Context, Result};

use crate::Genome;
use crate::bed::{BedReader, BedRecord};
use crate::order::InputOrder;

/// Header policy for [`spacing`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct SpacingOptions {
    /// Preserve initial input headers, including header-only streams.
    pub header: bool,
}

/// Append distance to the immediately preceding BED record on each chromosome.
///
/// First records receive `.`, qualifying overlaps `-1`, and other records the
/// signed gap. Points use symmetric one-base expansion for distance calculations;
/// overlap requires at least 1e-9 of both expanded spans in binary32 arithmetic.
/// An optional genome supplies chromosome order, not coordinate clipping.
///
/// # Errors
/// Returns an error for malformed or unsorted BED, a chromosome absent from the
/// selected genome, or an input/output failure.
pub fn spacing(
    input: impl Read,
    genome: Option<&Genome>,
    output: impl Write,
    options: SpacingOptions,
) -> Result<()> {
    let mut reader = if options.header {
        BedReader::new(input).preserve_header()
    } else {
        BedReader::new(input)
    };
    let mut current = reader.next_record()?;
    let mut output = BufWriter::with_capacity(65_536, output);
    if options.header {
        output
            .write_all(reader.header())
            .rs_context("writing BED header")?;
    }
    let mut order = InputOrder::new("spacing", genome);
    let mut previous: Option<BedRecord> = None;
    while let Some(record) = current {
        let new_chromosome = order.accept(&record)?;
        let distance = if new_chromosome {
            None
        } else {
            let previous = previous
                .as_ref()
                .expect("same chromosome has a previous record");
            let (start, end) = record.signed_bounds();
            let (previous_start, previous_end) = previous.signed_bounds();
            let overlap = end.min(previous_end) - start.max(previous_start);
            Some(
                if overlap as f32 / (end - start) as f32 >= 1e-9_f32
                    && overlap as f32 / (previous_end - previous_start) as f32 >= 1e-9_f32
                {
                    -1
                } else {
                    start - previous_end
                },
            )
        };
        record.write_prefix(&mut output, record.start(), record.end())?;
        match distance {
            Some(distance) => writeln!(output, "\t{distance}"),
            None => output.write_all(b"\t.\n"),
        }
        .rs_context("writing BED spacing")?;
        previous = Some(record);
        current = reader.next_record()?;
    }
    output.flush().rs_context("flushing BED spacing")
}
