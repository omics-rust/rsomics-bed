use std::fmt::Write as _;
use std::io::Cursor;
use std::path::Path;
use std::process::Command;

use rand_chacha::ChaCha8Rng;
use rand_core::{RngCore, SeedableRng};
use rsomics_bed::{StrandFilter, closest, cluster, window};
use tempfile::tempdir;

#[derive(Clone)]
struct Record {
    chrom: String,
    start: u64,
    end: u64,
    name: String,
    strand: char,
}

impl Record {
    fn line(&self, output: &mut String) {
        writeln!(
            output,
            "{}\t{}\t{}\t{}\t0\t{}",
            self.chrom, self.start, self.end, self.name, self.strand
        )
        .unwrap();
    }
}

fn seeded_inputs() -> (Vec<u8>, Vec<u8>) {
    let mut rng = ChaCha8Rng::seed_from_u64(0x5253_4f4d_4943_5342);
    let mut a = Vec::new();
    let mut b = Vec::new();

    for chrom_rank in 1..=5 {
        let chrom = format!("chr{chrom_rank}");
        for index in 0..20 {
            let start = index * 100 + 20 + rng.next_u64() % 10;
            let end = if (index + chrom_rank) % 11 == 0 {
                start
            } else {
                start + 10 + rng.next_u64() % 30
            };
            let strand = if rng.next_u64() & 1 == 0 { '+' } else { '-' };
            a.push(Record {
                chrom: chrom.clone(),
                start,
                end,
                name: format!("n{}", index % 7),
                strand,
            });
        }
    }

    for record in a.iter().filter(|record| record.chrom != "chr5") {
        let index = record.start / 100;
        let opposite = if record.strand == '+' { '-' } else { '+' };
        match index % 4 {
            0 => b.push(Record {
                chrom: record.chrom.clone(),
                start: record.start.saturating_sub(2),
                end: record.end.max(record.start + 1) + 2,
                name: record.name.clone(),
                strand: record.strand,
            }),
            1 => {
                let gap = 2 + rng.next_u64() % 8;
                let left_end = record.start - gap;
                b.push(Record {
                    chrom: record.chrom.clone(),
                    start: left_end - 5,
                    end: left_end,
                    name: format!("left-{index}"),
                    strand: record.strand,
                });
                b.push(Record {
                    chrom: record.chrom.clone(),
                    start: record.end + gap,
                    end: record.end + gap + 5,
                    name: format!("right-{index}"),
                    strand: opposite,
                });
            }
            2 => {
                b.push(Record {
                    chrom: record.chrom.clone(),
                    start: record.start,
                    end: record.start,
                    name: format!("point-{index}"),
                    strand: opposite,
                });
                b.push(Record {
                    chrom: record.chrom.clone(),
                    start: record.end + 3,
                    end: record.end + 3,
                    name: format!("point-right-{index}"),
                    strand: record.strand,
                });
            }
            _ => {
                let nested = Record {
                    chrom: record.chrom.clone(),
                    start: record.start - 3,
                    end: record.end.max(record.start + 1) + 3,
                    name: format!("nested-{index}"),
                    strand: opposite,
                };
                b.push(nested.clone());
                if index % 3 == 0 {
                    b.push(Record {
                        name: format!("duplicate-{index}"),
                        ..nested
                    });
                }
            }
        }
    }

    b.sort_by(|left, right| {
        (&left.chrom, left.start, left.end).cmp(&(&right.chrom, right.start, right.end))
    });
    let mut a_text = String::new();
    let mut b_text = String::new();
    for record in a {
        record.line(&mut a_text);
    }
    for record in b {
        record.line(&mut b_text);
    }
    (a_text.into_bytes(), b_text.into_bytes())
}

fn bedtools_closest(a: &Path, b: &Path, flags: &[&str]) -> Vec<u8> {
    let output = Command::new("bedtools")
        .arg("closest")
        .args(flags)
        .arg("-a")
        .arg(a)
        .arg("-b")
        .arg(b)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

fn bedtools_window(a: &Path, b: &Path, flags: &[&str]) -> Vec<u8> {
    let output = Command::new("bedtools")
        .arg("window")
        .args(flags)
        .arg("-a")
        .arg(a)
        .arg("-b")
        .arg(b)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

fn bedtools_cluster(input: &Path, flags: &[&str]) -> Vec<u8> {
    let output = Command::new("bedtools")
        .arg("cluster")
        .args(flags)
        .arg("-i")
        .arg(input)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}

fn our_closest(a: &[u8], b: &[u8], options: closest::ClosestOptions) -> Vec<u8> {
    let mut output = Vec::new();
    closest::closest(Cursor::new(a), Cursor::new(b), &mut output, options).unwrap();
    output
}

fn our_window(a: &[u8], b: &[u8], options: window::WindowOptions) -> Vec<u8> {
    let mut output = Vec::new();
    window::window(Cursor::new(a), Cursor::new(b), &mut output, options).unwrap();
    output
}

fn our_cluster(input: &[u8], options: cluster::ClusterOptions) -> Vec<u8> {
    let mut output = Vec::new();
    cluster::cluster(Cursor::new(input), &mut output, options).unwrap();
    output
}

#[test]
fn seeded_relation_operations_match_bedtools_231() {
    let version = Command::new("bedtools").arg("--version").output().unwrap();
    assert_eq!(version.stdout, b"bedtools v2.31.1\n");

    let (a, b) = seeded_inputs();
    let directory = tempdir().unwrap();
    let a_path = directory.path().join("a.bed");
    let b_path = directory.path().join("b.bed");
    std::fs::write(&a_path, &a).unwrap();
    std::fs::write(&b_path, &b).unwrap();

    let closest_default = closest::ClosestOptions::default();
    for (label, options, flags) in [
        ("default", closest_default, &[][..]),
        (
            "unsigned",
            closest::ClosestOptions {
                distance: closest::DistanceMode::Unsigned,
                ..closest_default
            },
            &["-d"][..],
        ),
        (
            "reference",
            closest::ClosestOptions {
                distance: closest::DistanceMode::Reference,
                ..closest_default
            },
            &["-D", "ref"][..],
        ),
        (
            "a-oriented",
            closest::ClosestOptions {
                distance: closest::DistanceMode::A,
                ..closest_default
            },
            &["-D", "a"][..],
        ),
        (
            "b-oriented",
            closest::ClosestOptions {
                distance: closest::DistanceMode::B,
                ..closest_default
            },
            &["-D", "b"][..],
        ),
        (
            "same-strand",
            closest::ClosestOptions {
                strand: StrandFilter::Same,
                ..closest_default
            },
            &["-s"][..],
        ),
        (
            "opposite-strand",
            closest::ClosestOptions {
                strand: StrandFilter::Opposite,
                ..closest_default
            },
            &["-S"][..],
        ),
        (
            "different-name",
            closest::ClosestOptions {
                different_name: true,
                ..closest_default
            },
            &["-N"][..],
        ),
        (
            "ignore-overlaps",
            closest::ClosestOptions {
                ignore_overlaps: true,
                ..closest_default
            },
            &["-io"][..],
        ),
        (
            "first-tie",
            closest::ClosestOptions {
                tie: closest::TieMode::First,
                ..closest_default
            },
            &["-t", "first"][..],
        ),
        (
            "last-tie",
            closest::ClosestOptions {
                tie: closest::TieMode::Last,
                ..closest_default
            },
            &["-t", "last"][..],
        ),
    ] {
        assert_eq!(
            our_closest(&a, &b, options),
            bedtools_closest(&a_path, &b_path, flags),
            "closest {label}"
        );
    }

    let window_default = window::WindowOptions {
        left: 12,
        right: 12,
        ..window::WindowOptions::default()
    };
    for (label, options, flags) in [
        ("pairs", window_default, &["-w", "12"][..]),
        (
            "asymmetric",
            window::WindowOptions {
                left: 7,
                right: 3,
                ..window::WindowOptions::default()
            },
            &["-l", "7", "-r", "3"][..],
        ),
        (
            "strand-relative",
            window::WindowOptions {
                left: 7,
                right: 3,
                strand_relative: true,
                strand: StrandFilter::Same,
                report: window::WindowReport::Pairs,
            },
            &["-l", "7", "-r", "3", "-sw", "-sm"][..],
        ),
        (
            "same-strand",
            window::WindowOptions {
                strand: StrandFilter::Same,
                ..window_default
            },
            &["-w", "12", "-sm"][..],
        ),
        (
            "opposite-strand",
            window::WindowOptions {
                strand: StrandFilter::Opposite,
                ..window_default
            },
            &["-w", "12", "-Sm"][..],
        ),
        (
            "any",
            window::WindowOptions {
                report: window::WindowReport::Any,
                ..window_default
            },
            &["-w", "12", "-u"][..],
        ),
        (
            "count",
            window::WindowOptions {
                report: window::WindowReport::Count,
                ..window_default
            },
            &["-w", "12", "-c"][..],
        ),
        (
            "none",
            window::WindowOptions {
                report: window::WindowReport::None,
                ..window_default
            },
            &["-w", "12", "-v"][..],
        ),
    ] {
        assert_eq!(
            our_window(&a, &b, options),
            bedtools_window(&a_path, &b_path, flags),
            "window {label}"
        );
    }

    for (label, options, flags) in [
        ("default", cluster::ClusterOptions::default(), &[][..]),
        (
            "distance",
            cluster::ClusterOptions {
                distance: 7,
                same_strand: false,
            },
            &["-d", "7"][..],
        ),
        (
            "same-strand",
            cluster::ClusterOptions {
                distance: 0,
                same_strand: true,
            },
            &["-s"][..],
        ),
    ] {
        assert_eq!(
            our_cluster(&a, options),
            bedtools_cluster(&a_path, flags),
            "cluster {label}"
        );
    }
}
