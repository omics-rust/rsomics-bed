//! BED subtraction with bedtools-compatible zero-length behavior.

use std::io::{BufWriter, Read, Write};

use rsomics_common::{Context, Result};

use crate::bed::{BedReader, BedRecord};
use crate::overlap_index::IndexedBed;

/// Remove B coverage from A while preserving A's trailing columns.
///
/// # Errors
///
/// Returns an error for malformed BED, I/O failures, coordinates outside the
/// pinned interval-index backend, or zero-length widening outside the `u64`
/// coordinate domain.
pub fn subtract(a_input: impl Read, b_input: impl Read, output: impl Write) -> Result<()> {
    let b = IndexedBed::load_for_subtract(b_input, "B")?;
    let mut a = BedReader::new(a_input);
    let mut output = BufWriter::new(output);
    let mut covers = Vec::new();

    while let Some(record) = a.next_record()? {
        b.ensure_query(&record, "A")?;
        subtract_record(&record, &b, &mut covers, &mut output)?;
    }
    output.flush().rs_context("flushing subtracted BED output")
}

fn subtract_record(
    a: &BedRecord,
    b: &IndexedBed,
    covers: &mut Vec<(u64, u64)>,
    output: &mut dyn Write,
) -> Result<()> {
    b.coverage_overlaps(a, covers);
    if a.start == a.end {
        let target = (a.start - 1, a.end + 1);
        if covers.as_slice() != [target] {
            a.write_raw(output)?;
        }
        return Ok(());
    }

    let mut cursor = a.start;
    for &(low, high) in covers.iter() {
        if low > cursor {
            a.write_with_coords(output, cursor, low)?;
        }
        cursor = cursor.max(high);
        if cursor >= a.end {
            return Ok(());
        }
    }
    if cursor < a.end {
        a.write_with_coords(output, cursor, a.end)?;
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn overlapping_b_intervals_emit_each_gap_once() {
        let a = b"chr1\t0\t100\tA\n";
        let b = b"chr1\t10\t30\nchr1\t20\t40\nchr1\t60\t70\n";
        let mut output = Vec::new();
        subtract(&a[..], &b[..], &mut output).unwrap();
        assert_eq!(
            output,
            b"chr1\t0\t10\tA\nchr1\t40\t60\tA\nchr1\t70\t100\tA\n"
        );
    }

    #[test]
    fn zero_length_b_removes_virtual_footprint() {
        let mut output = Vec::new();
        subtract(
            &b"chr1\t5\t15\tA\n"[..],
            &b"chr1\t10\t10\tB\n"[..],
            &mut output,
        )
        .unwrap();
        assert_eq!(output, b"chr1\t5\t9\tA\nchr1\t11\t15\tA\n");
    }

    #[test]
    fn zero_length_a_is_preserved_only_without_overlap() {
        let mut covered = Vec::new();
        subtract(
            &b"chr1\t10\t10\tA\n"[..],
            &b"chr1\t5\t15\tB\n"[..],
            &mut covered,
        )
        .unwrap();
        assert!(covered.is_empty());

        let mut uncovered = Vec::new();
        subtract(
            &b"chr1\t10\t10\tA\n"[..],
            &b"chr1\t50\t60\tB\n"[..],
            &mut uncovered,
        )
        .unwrap();
        assert_eq!(uncovered, b"chr1\t10\t10\tA\n");
    }

    #[test]
    fn zero_length_a_at_nonzero_b_boundaries_is_preserved() {
        let a = b"chr1\t9\t9\toutside-left\n\
                  chr1\t10\t10\tat-start\n\
                  chr1\t11\t11\tinside-left\n\
                  chr1\t19\t19\tinside-right\n\
                  chr1\t20\t20\tat-end\n\
                  chr1\t21\t21\toutside-right\n";
        let b = b"chr1\t10\t20\tB\n";
        let mut output = Vec::new();
        subtract(&a[..], &b[..], &mut output).unwrap();
        assert_eq!(
            output,
            b"chr1\t9\t9\toutside-left\n\
              chr1\t10\t10\tat-start\n\
              chr1\t20\t20\tat-end\n\
              chr1\t21\t21\toutside-right\n"
        );
    }

    #[test]
    fn zero_length_a_is_removed_only_by_zero_b_at_same_position() {
        let a = b"chr1\t9\t9\tleft\nchr1\t10\t10\tsame\nchr1\t11\t11\tright\n";
        let b = b"chr1\t10\t10\tB\n";
        let mut output = Vec::new();
        subtract(&a[..], &b[..], &mut output).unwrap();
        assert_eq!(output, b"chr1\t9\t9\tleft\nchr1\t11\t11\tright\n");
    }

    #[test]
    fn zero_length_a_between_touching_b_intervals_is_removed() {
        let a = b"chr1\t70\t70\tA\n";
        let b = b"chr1\t48\t70\tleft\nchr1\t70\t98\tright\n";
        let mut output = Vec::new();
        subtract(&a[..], &b[..], &mut output).unwrap();
        assert!(output.is_empty());
    }

    #[test]
    fn zero_length_a_uses_union_of_zero_and_nonzero_b_virtual_coverage() {
        let a = b"chr1\t3\t3\tA\n";
        let b = b"chr1\t2\t2\tleft-zero\nchr1\t3\t4\tright-nonzero\n";
        let mut output = Vec::new();
        subtract(&a[..], &b[..], &mut output).unwrap();
        assert!(output.is_empty());
    }

    #[test]
    fn zero_length_a_uses_union_of_adjacent_zero_b_virtual_coverage() {
        let a = b"chr1\t3\t3\tA\n";
        let b = b"chr1\t2\t2\tleft-zero\nchr1\t4\t4\tright-zero\n";
        let mut output = Vec::new();
        subtract(&a[..], &b[..], &mut output).unwrap();
        assert!(output.is_empty());
    }

    #[test]
    fn maximum_backend_coordinate_is_supported() {
        let limit = i32::MAX as u64;
        let maximum = limit - 1;
        let a = format!("chr1\t{maximum}\t{limit}\tA\n");
        let b = format!("chr1\t{maximum}\t{limit}\tB\n");
        let mut output = Vec::new();
        subtract(a.as_bytes(), b.as_bytes(), &mut output).unwrap();
        assert!(output.is_empty());
    }

    #[test]
    fn first_unrepresentable_coordinate_is_rejected_for_a_and_b() {
        let first_unrepresentable = i32::MAX as u64;
        let invalid = format!(
            "chr1\t{first_unrepresentable}\t{}\n",
            first_unrepresentable + 1
        );
        let b_error = subtract(&b"chr1\t0\t1\n"[..], invalid.as_bytes(), Vec::new()).unwrap_err();
        assert!(b_error.to_string().contains("B interval"), "{b_error}");

        let a_error = subtract(invalid.as_bytes(), &b"chr1\t0\t1\n"[..], Vec::new()).unwrap_err();
        assert!(a_error.to_string().contains("A interval"), "{a_error}");
    }
}
