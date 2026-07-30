//! BED intersection with bedtools-compatible clipping and zero-length behavior.

use std::io::{BufWriter, Read, Write};

use rsomics_common::{Context, Result};

use crate::bed::{BedReader, BedRecord};
use crate::overlap_index::IndexedBed;

/// Emit each clipped A/B intersection in B-file order.
///
/// Trailing columns from A are retained. Header, browser, track, comment, and
/// blank input lines are omitted, matching the default bedtools output policy.
///
/// # Errors
///
/// Returns an error for malformed BED, I/O failures, coordinates outside the
/// pinned interval-index backend, or zero-length widening outside the `u64`
/// coordinate domain.
pub fn intersect(a_input: impl Read, b_input: impl Read, output: impl Write) -> Result<()> {
    let b = IndexedBed::load(b_input, "B")?;
    let mut a = BedReader::new(a_input);
    let mut output = BufWriter::new(output);
    let mut candidate_ids = Vec::new();

    while let Some(record) = a.next_record()? {
        b.ensure_query(&record, "A")?;
        if record.start == record.end {
            emit_zero_length_a(&record, &b, &mut candidate_ids, &mut output)?;
        } else {
            emit_nonzero_a(&record, &b, &mut candidate_ids, &mut output)?;
        }
    }
    output.flush().rs_context("flushing intersect BED output")
}

fn emit_zero_length_a(
    a: &BedRecord,
    b: &IndexedBed,
    candidate_ids: &mut Vec<usize>,
    output: &mut dyn Write,
) -> Result<()> {
    b.intersection_candidates(a, candidate_ids);
    for _ in candidate_ids {
        a.write_raw(output)?;
    }
    Ok(())
}

fn emit_nonzero_a(
    a: &BedRecord,
    b: &IndexedBed,
    candidate_ids: &mut Vec<usize>,
    output: &mut dyn Write,
) -> Result<()> {
    b.intersection_candidates(a, candidate_ids);
    for &id in candidate_ids.iter() {
        let candidate = b.record(id);
        let (low, high) = if candidate.start == candidate.end {
            let (b_low, b_high) = candidate.virtual_bounds();
            (a.start.max(b_low), a.end.min(b_high))
        } else {
            (a.start.max(candidate.start), a.end.min(candidate.end))
        };
        if low < high {
            a.write_with_coords(output, low, high)?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn emits_hits_in_b_file_order() {
        let a = b"chr1\t0\t100\tA\n";
        let b = b"chr1\t50\t60\nchr1\t10\t20\nchr1\t30\t40\n";
        let mut output = Vec::new();
        intersect(&a[..], &b[..], &mut output).unwrap();
        assert_eq!(
            output,
            b"chr1\t50\t60\tA\nchr1\t10\t20\tA\nchr1\t30\t40\tA\n"
        );
    }

    #[test]
    fn duplicate_b_coordinates_preserve_multiplicity_and_file_order() {
        let a = b"chr1\t0\t100\tA\n";
        let b = b"chr1\t10\t20\tfirst\nchr1\t30\t40\tmiddle\nchr1\t10\t20\tlast\n";
        let mut output = Vec::new();
        intersect(&a[..], &b[..], &mut output).unwrap();
        assert_eq!(
            output,
            b"chr1\t10\t20\tA\nchr1\t30\t40\tA\nchr1\t10\t20\tA\n"
        );
    }

    #[test]
    fn zero_length_b_uses_bedtools_virtual_footprint() {
        let mut output = Vec::new();
        intersect(
            &b"chr1\t5\t15\tA\n"[..],
            &b"chr1\t10\t10\tB\n"[..],
            &mut output,
        )
        .unwrap();
        assert_eq!(output, b"chr1\t9\t11\tA\n");
    }

    #[test]
    fn maximum_backend_coordinate_is_supported() {
        let end_exclusive_limit = i32::MAX as u64;
        let maximum = end_exclusive_limit - 1;
        let a = format!("chr1\t{maximum}\t{end_exclusive_limit}\tA\n");
        let b = format!("chr1\t{maximum}\t{end_exclusive_limit}\tB\n");
        let mut output = Vec::new();
        intersect(a.as_bytes(), b.as_bytes(), &mut output).unwrap();
        assert_eq!(
            output,
            format!("chr1\t{maximum}\t{end_exclusive_limit}\tA\n").as_bytes()
        );
    }

    #[test]
    fn first_unrepresentable_backend_coordinate_is_rejected() {
        let first_unrepresentable = i32::MAX as u64;
        let b = format!(
            "chr1\t{first_unrepresentable}\t{}\tB\n",
            first_unrepresentable + 1
        );
        let error = intersect(&b"chr1\t0\t1\tA\n"[..], b.as_bytes(), Vec::new()).unwrap_err();
        assert!(error.to_string().contains("maximum inclusive coordinate"));
    }

    #[test]
    fn first_unrepresentable_a_query_coordinate_is_rejected() {
        let first_unrepresentable = i32::MAX as u64;
        let a = format!(
            "chr1\t{first_unrepresentable}\t{}\tA\n",
            first_unrepresentable + 1
        );
        let error = intersect(a.as_bytes(), &b"chr1\t0\t1\tB\n"[..], Vec::new()).unwrap_err();
        assert!(error.to_string().contains("A interval"), "{error}");
        assert!(error.to_string().contains("maximum inclusive coordinate"));
    }

    #[test]
    fn maximum_safe_zero_length_b_candidate_is_supported() {
        let limit = i32::MAX as u64;
        let maximum = limit - 1;
        let a = format!("chr1\t{}\t{limit}\tA\n", maximum - 1);
        let b = format!("chr1\t{maximum}\t{maximum}\tB\n");
        let mut output = Vec::new();
        intersect(a.as_bytes(), b.as_bytes(), &mut output).unwrap();
        assert_eq!(
            output,
            format!("chr1\t{}\t{limit}\tA\n", maximum - 1).as_bytes()
        );
    }
}
