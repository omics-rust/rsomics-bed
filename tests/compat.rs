use std::fs::File;
use std::path::{Path, PathBuf};
use std::process::Command;

use rsomics_bed::{cluster, complement, intersect, merge, read_genome, sort, subtract};

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

fn bytes(name: &str) -> Vec<u8> {
    std::fs::read(golden(name)).unwrap()
}

fn sort_output(input: &str) -> Vec<u8> {
    let mut output = Vec::new();
    sort::sort(File::open(golden(input)).unwrap(), &mut output).unwrap();
    output
}

fn merge_output(input: &str) -> Vec<u8> {
    let mut output = Vec::new();
    merge::merge(File::open(golden(input)).unwrap(), &mut output).unwrap();
    output
}

fn intersect_output(a: &str, b: &str) -> Vec<u8> {
    let mut output = Vec::new();
    intersect::intersect(
        File::open(golden(a)).unwrap(),
        File::open(golden(b)).unwrap(),
        &mut output,
    )
    .unwrap();
    output
}

fn subtract_output(a: &str, b: &str) -> Vec<u8> {
    let mut output = Vec::new();
    subtract::subtract(
        File::open(golden(a)).unwrap(),
        File::open(golden(b)).unwrap(),
        &mut output,
    )
    .unwrap();
    output
}

fn complement_output(input: &str) -> Vec<u8> {
    let genome = read_genome(File::open(golden("complement.genome.tsv")).unwrap()).unwrap();
    let mut output = Vec::new();
    complement::complement(File::open(golden(input)).unwrap(), &genome, &mut output).unwrap();
    output
}

fn cluster_output(input: &str, distance: u64, same_strand: bool) -> Vec<u8> {
    let mut output = Vec::new();
    cluster::cluster(
        File::open(golden(input)).unwrap(),
        &mut output,
        cluster::ClusterOptions {
            distance,
            same_strand,
        },
    )
    .unwrap();
    output
}

#[test]
fn committed_golden_outputs_match() {
    assert_eq!(
        cluster_output("cluster.input.bed", 0, false),
        bytes("cluster.default.expected.bed")
    );
    assert_eq!(
        cluster_output("cluster.input.bed", 5, false),
        bytes("cluster.distance.expected.bed")
    );
    assert_eq!(
        cluster_output("cluster.strand.input.bed", 0, true),
        bytes("cluster.strand.expected.bed")
    );
    assert_eq!(sort_output("sort.input.bed"), bytes("sort.expected.bed"));
    assert_eq!(
        sort_output("sort.headers.input.bed"),
        bytes("sort.headers.expected.bed")
    );
    assert_eq!(merge_output("merge.input.bed"), bytes("merge.expected.bed"));
    assert_eq!(
        merge_output("merge.zero.input.bed"),
        bytes("merge.zero.expected.bed")
    );
    assert_eq!(
        intersect_output("a.bed", "b.bed"),
        bytes("intersect.expected.bed")
    );
    assert_eq!(
        intersect_output("intersect.unsorted-a.bed", "intersect.unsorted-b.bed"),
        bytes("intersect.unsorted.expected.bed")
    );
    assert_eq!(
        intersect_output("intersect.zero-a.bed", "intersect.zero-b.bed"),
        bytes("intersect.zero.expected.bed")
    );
    assert_eq!(
        subtract_output("a.bed", "b.bed"),
        bytes("subtract.expected.bed")
    );
    assert_eq!(
        subtract_output(
            "subtract.zero-boundary-a.bed",
            "subtract.zero-boundary-b.bed"
        ),
        bytes("subtract.zero-boundary.expected.bed")
    );
    assert_eq!(
        complement_output("complement.input.bed"),
        bytes("complement.expected.bed")
    );
    assert_eq!(
        complement_output("complement.zero.input.bed"),
        bytes("complement.zero.expected.bed")
    );
}

#[test]
fn inherited_failure_fixtures_fail_loud() {
    let huge = intersect::intersect(
        File::open(golden("intersect.huge-a.bed")).unwrap(),
        File::open(golden("intersect.huge-b.bed")).unwrap(),
        Vec::new(),
    )
    .unwrap_err();
    assert!(huge.to_string().contains("backend maximum"), "{huge}");

    let malformed = intersect::intersect(
        File::open(golden("intersect.malformed-a.bed")).unwrap(),
        File::open(golden("b.bed")).unwrap(),
        Vec::new(),
    )
    .unwrap_err();
    assert!(
        malformed.to_string().contains("greater than end"),
        "{malformed}"
    );

    let intersect_origin = intersect::intersect(
        File::open(golden("indexed.origin-zero-a.bed")).unwrap(),
        File::open(golden("indexed.origin-cover.bed")).unwrap(),
        Vec::new(),
    )
    .unwrap_err();
    assert!(
        intersect_origin
            .to_string()
            .contains("A zero-length interval chr1:0-0"),
        "{intersect_origin}"
    );

    let subtract_origin = subtract::subtract(
        File::open(golden("indexed.origin-zero-a.bed")).unwrap(),
        File::open(golden("indexed.origin-cover.bed")).unwrap(),
        Vec::new(),
    )
    .unwrap_err();
    assert!(
        subtract_origin
            .to_string()
            .contains("A zero-length interval chr1:0-0"),
        "{subtract_origin}"
    );

    let origin_zero = merge::merge(
        File::open(golden("merge.origin-zero.input.bed")).unwrap(),
        Vec::new(),
    )
    .unwrap_err();
    assert!(
        origin_zero.to_string().contains("widens below"),
        "{origin_zero}"
    );

    let genome = read_genome(File::open(golden("complement.genome.tsv")).unwrap()).unwrap();
    for name in [
        "complement.unsorted.bed",
        "complement.absent-chrom.bed",
        "complement.inverted.bed",
    ] {
        let error = complement::complement(File::open(golden(name)).unwrap(), &genome, Vec::new())
            .unwrap_err();
        assert!(!error.to_string().is_empty(), "{name}");
    }
}

#[test]
fn live_bedtools_231_compatibility() {
    let version = Command::new("bedtools")
        .arg("--version")
        .output()
        .expect("bedtools 2.31.1 is required for the compatibility lane");
    assert!(
        version.status.success(),
        "bedtools 2.31.1 is required for the compatibility lane"
    );
    let version = String::from_utf8(version.stdout).unwrap();
    assert!(
        version.contains("v2.31.1"),
        "compatibility oracle changed: {version}"
    );

    let cases = [
        (
            cluster_output("cluster.input.bed", 0, false),
            bedtools(&["cluster", "-i"], &["cluster.input.bed"]),
            "cluster",
        ),
        (
            cluster_output("cluster.input.bed", 5, false),
            bedtools(&["cluster", "-d", "5", "-i"], &["cluster.input.bed"]),
            "cluster distance",
        ),
        (
            cluster_output("cluster.strand.input.bed", 0, true),
            bedtools(&["cluster", "-s", "-i"], &["cluster.strand.input.bed"]),
            "cluster strand",
        ),
        (
            sort_output("sort.input.bed"),
            bedtools(&["sort", "-i"], &["sort.input.bed"]),
            "sort",
        ),
        (
            sort_output("sort.headers.input.bed"),
            bedtools(&["sort", "-i"], &["sort.headers.input.bed"]),
            "sort headers and trailing columns",
        ),
        (
            merge_output("merge.input.bed"),
            bedtools(&["merge", "-i"], &["merge.input.bed"]),
            "merge",
        ),
        (
            merge_output("merge.zero.input.bed"),
            bedtools(&["merge", "-i"], &["merge.zero.input.bed"]),
            "merge zero-length",
        ),
        (
            intersect_output("a.bed", "b.bed"),
            bedtools(&["intersect", "-a"], &["a.bed", "-b", "b.bed"]),
            "intersect",
        ),
        (
            intersect_output("intersect.zero-a.bed", "intersect.zero-b.bed"),
            bedtools(
                &["intersect", "-a"],
                &["intersect.zero-a.bed", "-b", "intersect.zero-b.bed"],
            ),
            "intersect zero-length",
        ),
        (
            subtract_output("a.bed", "b.bed"),
            bedtools(&["subtract", "-a"], &["a.bed", "-b", "b.bed"]),
            "subtract",
        ),
        (
            subtract_output(
                "subtract.zero-boundary-a.bed",
                "subtract.zero-boundary-b.bed",
            ),
            bedtools(
                &["subtract", "-a"],
                &[
                    "subtract.zero-boundary-a.bed",
                    "-b",
                    "subtract.zero-boundary-b.bed",
                ],
            ),
            "subtract zero-length boundaries",
        ),
        (
            complement_output("complement.input.bed"),
            bedtools(
                &["complement", "-i"],
                &["complement.input.bed", "-g", "complement.genome.tsv"],
            ),
            "complement",
        ),
        (
            complement_output("complement.zero.input.bed"),
            bedtools(
                &["complement", "-i"],
                &["complement.zero.input.bed", "-g", "complement.genome.tsv"],
            ),
            "complement zero-length",
        ),
    ];

    for (ours, upstream, operation) in cases {
        assert_eq!(ours, upstream, "{operation} differs from bedtools 2.31.1");
    }

    assert_eq!(
        bedtools(
            &["intersect", "-a"],
            &[
                "indexed.origin-zero-a.bed",
                "-b",
                "indexed.origin-cover.bed"
            ]
        ),
        bytes("indexed.origin-zero-a.bed"),
        "bedtools accepts origin zero-length A for intersect; rsomics-bed intentionally rejects it"
    );
    assert_eq!(
        bedtools(
            &["subtract", "-a"],
            &[
                "indexed.origin-zero-a.bed",
                "-b",
                "indexed.origin-cover.bed"
            ]
        ),
        bytes("indexed.origin-zero-a.bed"),
        "bedtools accepts origin zero-length A for subtract; rsomics-bed intentionally rejects it"
    );
}

fn bedtools(prefix: &[&str], args: &[&str]) -> Vec<u8> {
    let mut command = Command::new("bedtools");
    command.args(prefix);
    for arg in args {
        if arg.starts_with('-') {
            command.arg(arg);
        } else {
            command.arg(golden(arg));
        }
    }
    let output = command.output().unwrap();
    assert!(
        output.status.success(),
        "bedtools failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    output.stdout
}
