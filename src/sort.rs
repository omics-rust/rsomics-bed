//! Stable BED coordinate sorting.

use std::io::{BufWriter, Read, Write};

use rsomics_common::{Context, Result};

use crate::bed::read_records;

/// Sort BED records by chromosome and start while preserving equal-key order.
///
/// All original BED columns are emitted byte-for-byte apart from normalized
/// line endings. Header, browser, track, comment, and blank lines are omitted.
///
/// # Errors
///
/// Returns an error for malformed BED or an input/output failure.
pub fn sort(input: impl Read, output: impl Write) -> Result<()> {
    let mut records = read_records(input)?;
    records.sort_by(|left, right| {
        left.chrom()
            .cmp(right.chrom())
            .then(left.start().cmp(&right.start()))
    });

    let mut output = BufWriter::new(output);
    for record in records {
        record.write_raw(&mut output)?;
    }
    output.flush().rs_context("flushing sorted BED output")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stable_sort_preserves_columns_and_equal_key_order() {
        let input = b"chr2\t1\t2\tz\nchr1\t5\t9\tfirst\nchr1\t5\t7\tsecond\n";
        let mut output = Vec::new();
        sort(&input[..], &mut output).unwrap();
        assert_eq!(
            output,
            b"chr1\t5\t9\tfirst\nchr1\t5\t7\tsecond\nchr2\t1\t2\tz\n"
        );
    }
}
