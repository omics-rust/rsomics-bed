//! Indexed nearest-neighbor queries for BED records.

use std::io::{BufWriter, Read, Write};

use rsomics_common::{Context, Result};

use crate::StrandFilter;
use crate::bed::{BedReader, BedRecord, Strand, virtual_bounds};
use crate::relation_index::RelationBed;

/// Distance representation appended to each closest pair.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum DistanceMode {
    /// Do not append a distance.
    #[default]
    None,
    /// Append the non-negative BEDTools distance.
    Unsigned,
    /// Sign distance by reference-coordinate direction.
    Reference,
    /// Sign distance relative to A's strand.
    A,
    /// Sign distance relative to each B record's strand.
    B,
}

/// Policy for equally close B records.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum TieMode {
    /// Emit every tie.
    #[default]
    All,
    /// Emit the first tie in the selected distance ordering.
    First,
    /// Emit the last tie in the selected distance ordering.
    Last,
}

/// Options for [`closest`].
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub struct ClosestOptions {
    /// Strand relationship required between A and B.
    pub strand: StrandFilter,
    /// Exclude B records with the same BED name as A.
    pub different_name: bool,
    /// Exclude B records that overlap A.
    pub ignore_overlaps: bool,
    /// Distance representation.
    pub distance: DistanceMode,
    /// Equal-distance selection policy.
    pub tie: TieMode,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    Left,
    Overlap,
    Right,
}

#[derive(Debug, Clone, Copy)]
struct Candidate {
    id: usize,
    unsigned: u64,
    signed: i128,
}

struct RequiredBFields {
    strands: Option<Vec<Strand>>,
    names: Option<Vec<Vec<u8>>>,
}

struct Query<'a> {
    fields: &'a RequiredBFields,
    options: ClosestOptions,
    bounds: (u64, u64),
    strand: Option<Strand>,
    name: Option<&'a [u8]>,
}

impl Query<'_> {
    fn eligible(&self, id: usize) -> bool {
        self.fields
            .eligible(id, self.options, self.strand, self.name)
    }

    fn candidate(&self, b: &RelationBed, id: usize) -> Candidate {
        let (unsigned, side) = classify_distance(self.bounds, b.interval(id).virtual_bounds());
        Candidate {
            id,
            unsigned,
            signed: signed_distance(
                unsigned,
                side,
                self.options.distance,
                self.strand,
                self.fields.strand(id),
            ),
        }
    }
}

impl RequiredBFields {
    fn load(b: &RelationBed, options: ClosestOptions) -> Result<Self> {
        let strands = if options.strand != StrandFilter::Any || options.distance == DistanceMode::B
        {
            Some(b.checked_strands("closest B")?)
        } else {
            None
        };
        let names = if options.different_name {
            Some(
                (0..b.len())
                    .map(|id| b.record(id).name("closest B").map(<[u8]>::to_vec))
                    .collect::<Result<Vec<_>>>()?,
            )
        } else {
            None
        };
        Ok(Self { strands, names })
    }

    fn eligible(
        &self,
        id: usize,
        options: ClosestOptions,
        a_strand: Option<Strand>,
        a_name: Option<&[u8]>,
    ) -> bool {
        let strand_ok = match options.strand {
            StrandFilter::Any => true,
            StrandFilter::Same => {
                self.strands.as_ref().expect("B strands were loaded")[id]
                    == a_strand.expect("A strand was loaded")
            }
            StrandFilter::Opposite => {
                self.strands.as_ref().expect("B strands were loaded")[id]
                    != a_strand.expect("A strand was loaded")
            }
        };
        let name_ok = !options.different_name
            || self.names.as_ref().expect("B names were loaded")[id]
                != a_name.expect("A name was loaded");
        strand_ok && name_ok
    }

    fn strand(&self, id: usize) -> Option<Strand> {
        self.strands.as_ref().map(|strands| strands[id])
    }
}

/// Emit the closest eligible B record or records for each A record.
///
/// # Errors
///
/// Returns an error for malformed BED, inconsistent B field widths, missing or
/// invalid fields required by the selected modes, unrepresentable indexed
/// coordinates, or I/O failure.
pub fn closest(
    a_input: impl Read,
    b_input: impl Read,
    output: impl Write,
    options: ClosestOptions,
) -> Result<()> {
    let b = RelationBed::load(b_input, "closest B")?;
    let fields = RequiredBFields::load(&b, options)?;
    let missing = missing_b_fields(b.field_count());
    let mut a = BedReader::new(a_input);
    let mut output = BufWriter::new(output);
    let mut overlap_ids = Vec::new();
    let mut candidates = Vec::new();

    while let Some(record) = a.next_record()? {
        let a_strand = if options.strand != StrandFilter::Any || options.distance == DistanceMode::A
        {
            Some(record.strand("closest A")?)
        } else {
            None
        };
        let a_name = if options.different_name {
            Some(record.name("closest A")?)
        } else {
            None
        };
        let a_bounds = virtual_bounds(&record, "closest A")?;
        let query = Query {
            fields: &fields,
            options,
            bounds: a_bounds,
            strand: a_strand,
            name: a_name,
        };
        b.range_candidates(
            record.chrom(),
            a_bounds.0,
            a_bounds.1,
            "closest A interval index",
            &mut overlap_ids,
        )?;

        candidates.clear();
        if !options.ignore_overlaps {
            candidates.extend(
                overlap_ids
                    .iter()
                    .copied()
                    .filter(|&id| query.eligible(id))
                    .map(|id| Candidate {
                        id,
                        unsigned: 0,
                        signed: 0,
                    }),
            );
        }

        if candidates.is_empty() {
            let mut best = None;
            collect_direction(
                &b,
                b.left_candidates(record.chrom(), a_bounds.0),
                &query,
                &mut best,
                &mut candidates,
            );
            collect_direction(
                &b,
                b.right_candidates(record.chrom(), a_bounds.1),
                &query,
                &mut best,
                &mut candidates,
            );
        }

        order_candidates(&mut candidates, options.distance);
        write_result(&record, &b, &candidates, &missing, &mut output, options)?;
    }
    output.flush().rs_context("flushing closest BED output")
}

fn collect_direction(
    b: &RelationBed,
    ids: impl Iterator<Item = usize>,
    query: &Query<'_>,
    best: &mut Option<u64>,
    candidates: &mut Vec<Candidate>,
) {
    for id in ids {
        let candidate = query.candidate(b, id);
        if best.is_some_and(|distance| candidate.unsigned > distance) {
            break;
        }
        if !query.eligible(id) {
            continue;
        }
        if best.is_none_or(|distance| candidate.unsigned < distance) {
            *best = Some(candidate.unsigned);
            candidates.clear();
        }
        if *best == Some(candidate.unsigned) {
            candidates.push(candidate);
        }
    }
}

fn classify_distance(a: (u64, u64), b: (u64, u64)) -> (u64, Side) {
    if a.0 < b.1 && b.0 < a.1 {
        (0, Side::Overlap)
    } else if b.1 <= a.0 {
        (a.0 - b.1 + 1, Side::Left)
    } else {
        (b.0 - a.1 + 1, Side::Right)
    }
}

fn signed_distance(
    unsigned: u64,
    side: Side,
    mode: DistanceMode,
    a_strand: Option<Strand>,
    b_strand: Option<Strand>,
) -> i128 {
    let unsigned = i128::from(unsigned);
    let reference = match side {
        Side::Left => -unsigned,
        Side::Overlap => 0,
        Side::Right => unsigned,
    };
    match mode {
        DistanceMode::None | DistanceMode::Unsigned | DistanceMode::Reference => reference,
        DistanceMode::A => match a_strand.expect("A strand was loaded") {
            Strand::Forward => reference,
            Strand::Reverse => -reference,
        },
        DistanceMode::B => match b_strand.expect("B strand was loaded") {
            Strand::Forward => -reference,
            Strand::Reverse => reference,
        },
    }
}

fn order_candidates(candidates: &mut [Candidate], mode: DistanceMode) {
    if matches!(
        mode,
        DistanceMode::Reference | DistanceMode::A | DistanceMode::B
    ) {
        candidates.sort_unstable_by_key(|candidate| (candidate.signed, candidate.id));
    } else {
        candidates.sort_unstable_by_key(|candidate| candidate.id);
    }
}

fn write_result(
    a: &BedRecord,
    b: &RelationBed,
    candidates: &[Candidate],
    missing: &str,
    output: &mut dyn Write,
    options: ClosestOptions,
) -> Result<()> {
    if candidates.is_empty() {
        if options.distance == DistanceMode::None {
            return a.write_column(output, missing);
        }
        return a.write_column(output, format_args!("{missing}\t-1"));
    }

    let selected = match options.tie {
        TieMode::All => candidates,
        TieMode::First => &candidates[..1],
        TieMode::Last => &candidates[candidates.len() - 1..],
    };
    for candidate in selected {
        match options.distance {
            DistanceMode::None => {
                a.write_joined_raw(output, b.record(candidate.id).raw())?;
            }
            DistanceMode::Unsigned => {
                a.write_joined_raw_column(
                    output,
                    b.record(candidate.id).raw(),
                    candidate.unsigned,
                )?;
            }
            DistanceMode::Reference | DistanceMode::A | DistanceMode::B => {
                a.write_joined_raw_column(output, b.record(candidate.id).raw(), candidate.signed)?;
            }
        }
    }
    Ok(())
}

fn missing_b_fields(field_count: usize) -> String {
    let mut output = String::from(".");
    for index in 1..field_count {
        output.push('\t');
        if index == 1 || index == 2 || (index == 4 && matches!(field_count, 5 | 6 | 12)) {
            output.push_str("-1");
        } else {
            output.push('.');
        }
    }
    output
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::StrandFilter;

    fn run(a: &[u8], b: &[u8], options: ClosestOptions) -> rsomics_common::Result<Vec<u8>> {
        let mut output = Vec::new();
        closest(a, b, &mut output, options)?;
        Ok(output)
    }

    #[test]
    fn distances_distinguish_overlap_bookends_and_direction() {
        assert_eq!(classify_distance((10, 20), (15, 25)), (0, Side::Overlap));
        assert_eq!(classify_distance((10, 20), (20, 30)), (1, Side::Right));
        assert_eq!(classify_distance((10, 20), (0, 5)), (6, Side::Left));
        assert_eq!(classify_distance((10, 20), (25, 30)), (6, Side::Right));
    }

    #[test]
    fn signed_modes_follow_reference_a_and_b_orientation() {
        assert_eq!(
            signed_distance(6, Side::Left, DistanceMode::Reference, None, None),
            -6
        );
        assert_eq!(
            signed_distance(6, Side::Right, DistanceMode::Reference, None, None),
            6
        );
        assert_eq!(
            signed_distance(6, Side::Left, DistanceMode::A, Some(Strand::Reverse), None),
            6
        );
        assert_eq!(
            signed_distance(6, Side::Left, DistanceMode::B, None, Some(Strand::Forward)),
            6
        );
        assert_eq!(
            signed_distance(6, Side::Left, DistanceMode::B, None, Some(Strand::Reverse)),
            -6
        );
    }

    #[test]
    fn report_modes_preserve_ties_and_signed_order() {
        let a = b"chr1\t10\t20\tA\t0\t+\n";
        let b = b"chr1\t0\t5\tleft-plus\t0\t+\n\
                  chr1\t0\t5\tleft-minus\t0\t-\n\
                  chr1\t25\t30\tright-plus\t0\t+\n";
        let default = run(a, b, ClosestOptions::default()).unwrap();
        assert_eq!(
            default,
            b"chr1\t10\t20\tA\t0\t+\tchr1\t0\t5\tleft-plus\t0\t+\n\
              chr1\t10\t20\tA\t0\t+\tchr1\t0\t5\tleft-minus\t0\t-\n\
              chr1\t10\t20\tA\t0\t+\tchr1\t25\t30\tright-plus\t0\t+\n"
        );

        let signed = run(
            a,
            b,
            ClosestOptions {
                distance: DistanceMode::B,
                ..ClosestOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            signed,
            b"chr1\t10\t20\tA\t0\t+\tchr1\t0\t5\tleft-minus\t0\t-\t-6\n\
              chr1\t10\t20\tA\t0\t+\tchr1\t25\t30\tright-plus\t0\t+\t-6\n\
              chr1\t10\t20\tA\t0\t+\tchr1\t0\t5\tleft-plus\t0\t+\t6\n"
        );

        let first = run(
            a,
            b,
            ClosestOptions {
                distance: DistanceMode::B,
                tie: TieMode::First,
                ..ClosestOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            first,
            b"chr1\t10\t20\tA\t0\t+\tchr1\t0\t5\tleft-minus\t0\t-\t-6\n"
        );

        let last = run(
            a,
            b,
            ClosestOptions {
                distance: DistanceMode::B,
                tie: TieMode::Last,
                ..ClosestOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            last,
            b"chr1\t10\t20\tA\t0\t+\tchr1\t0\t5\tleft-plus\t0\t+\t6\n"
        );
    }

    #[test]
    fn eligibility_precedes_nearest_selection() {
        let a = b"chr1\t10\t20\tsame\t0\t+\n";
        let b = b"chr1\t12\t18\tsame\t0\t+\n\
                  chr1\t25\t30\tother\t0\t-\n";
        let output = run(
            a,
            b,
            ClosestOptions {
                strand: StrandFilter::Opposite,
                different_name: true,
                ..ClosestOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            output,
            b"chr1\t10\t20\tsame\t0\t+\tchr1\t25\t30\tother\t0\t-\n"
        );

        let ignored = run(
            a,
            b,
            ClosestOptions {
                ignore_overlaps: true,
                ..ClosestOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            ignored,
            b"chr1\t10\t20\tsame\t0\t+\tchr1\t25\t30\tother\t0\t-\n"
        );

        let farther_eligible = run(
            a,
            b"chr1\t20\t21\tnear\t0\t+\nchr1\t31\t40\tfar\t0\t-\n",
            ClosestOptions {
                strand: StrandFilter::Opposite,
                ..ClosestOptions::default()
            },
        )
        .unwrap();
        assert_eq!(
            farther_eligible,
            b"chr1\t10\t20\tsame\t0\t+\tchr1\t31\t40\tfar\t0\t-\n"
        );

        let nearer_right = run(
            b"chr1\t10\t20\tA\n",
            b"chr1\t0\t1\tleft\nchr1\t20\t21\tright\n",
            ClosestOptions::default(),
        )
        .unwrap();
        assert_eq!(nearer_right, b"chr1\t10\t20\tA\tchr1\t20\t21\tright\n");
    }

    #[test]
    fn selected_fields_fail_with_side_and_physical_line() {
        let missing_a = run(
            b"chr1\t10\t20\n",
            b"chr1\t0\t5\tB\t0\t+\n",
            ClosestOptions {
                different_name: true,
                ..ClosestOptions::default()
            },
        )
        .unwrap_err();
        assert!(
            missing_a.to_string().contains("closest A BED line 1"),
            "{missing_a}"
        );

        let invalid_b = run(
            b"chr1\t10\t20\tA\t0\t+\n",
            b"chr2\t0\t5\tB\t0\t.\n",
            ClosestOptions {
                distance: DistanceMode::B,
                ..ClosestOptions::default()
            },
        )
        .unwrap_err();
        assert!(
            invalid_b.to_string().contains("closest B BED line 1"),
            "{invalid_b}"
        );
    }

    #[test]
    fn no_hit_uses_b_width_and_optional_distance() {
        let output = run(
            b"chr1\t10\t20\tA\n",
            b"chr2\t0\t5\tB\t7\t+\n",
            ClosestOptions {
                distance: DistanceMode::Unsigned,
                ..ClosestOptions::default()
            },
        )
        .unwrap();
        assert_eq!(output, b"chr1\t10\t20\tA\t.\t-1\t-1\t.\t-1\t.\t-1\n");

        let empty = run(b"chr1\t10\t20\n", b"", ClosestOptions::default()).unwrap();
        assert_eq!(empty, b"chr1\t10\t20\t.\t-1\t-1\n");
    }
}
