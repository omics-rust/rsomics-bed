use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use rsomics_bed::genomecov::{GenomecovMode, GenomecovOptions, genomecov};
use rsomics_bed::read_genome;
use tempfile::tempdir;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/genomecov")
        .join(name)
}

fn run(input: &str, flags: &[&str], upstream: bool) -> std::process::Output {
    let mut command = Command::new(if upstream {
        "bedtools"
    } else {
        env!("CARGO_BIN_EXE_rsomics-bed")
    });
    command.args(["genomecov", "-i"]).arg(fixture(input));
    command.arg("-g").arg(fixture("genome.txt"));
    if upstream {
        let translated: Vec<_> = flags
            .iter()
            .map(|flag| match *flag {
                "--mode" => unreachable!(),
                "--bedgraph" => "-bg",
                "--bedgraph-all" => "-bga",
                "--depth" => "-d",
                "--depth-nonzero" => "-dz",
                "--strand" => "-strand",
                "--five-prime" => "-5",
                "--three-prime" => "-3",
                "--scale" => "-scale",
                "--max" => "-max",
                "--split" => "-split",
                "--track" => "-trackopts",
                value => value,
            })
            .collect();
        command.args(translated);
    } else {
        command.args(flags);
    }
    command.output().unwrap()
}

#[test]
fn genomecov_full_outputs_preserve_points_order_bins_and_blocks() {
    for (input, flags, golden) in [
        ("ordinary.bed", vec![], "histogram.stdout"),
        (
            "ordinary.bed",
            vec!["--bedgraph-all"],
            "bedgraph-all.stdout",
        ),
        (
            "ordinary.bed",
            vec!["--bedgraph", "--max", "1"],
            "max-bedgraph.stdout",
        ),
        (
            "blocks.bed",
            vec!["--bedgraph-all", "--split"],
            "split.stdout",
        ),
    ] {
        let output = run(input, &flags, false);
        assert!(
            output.status.success(),
            "{golden}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(output.stderr.is_empty());
        assert_eq!(
            output.stdout,
            fs::read(fixture(golden)).unwrap(),
            "{golden}"
        );
    }
}

#[test]
fn genomecov_thirteen_workflows_match_live_bedtools_231() {
    assert_eq!(
        Command::new("bedtools")
            .arg("--version")
            .output()
            .unwrap()
            .stdout,
        b"bedtools v2.31.1\n"
    );
    let cases: &[(&str, &[&str])] = &[
        ("ordinary.bed", &[]),
        ("ordinary.bed", &["--bedgraph"]),
        ("ordinary.bed", &["--bedgraph-all"]),
        ("ordinary.bed", &["--depth"]),
        ("ordinary.bed", &["--depth-nonzero"]),
        ("ordinary.bed", &["--bedgraph-all", "--strand", "+"]),
        ("ordinary.bed", &["--bedgraph-all", "--five-prime"]),
        ("ordinary.bed", &["--bedgraph-all", "--three-prime"]),
        ("ordinary.bed", &["--bedgraph", "--scale", "0.1"]),
        ("ordinary.bed", &["--bedgraph", "--track", "name=coverage"]),
        ("ordinary.bed", &["--max", "1"]),
        ("ordinary.bed", &["--bedgraph", "--max", "1"]),
        ("blocks.bed", &["--bedgraph-all", "--split"]),
    ];
    for &(input, flags) in cases {
        let source = run(input, flags, true);
        let native = run(input, flags, false);
        assert!(source.status.success(), "{input} {flags:?}");
        assert!(
            native.status.success(),
            "{}",
            String::from_utf8_lossy(&native.stderr)
        );
        assert!(source.stderr.is_empty() && native.stderr.is_empty());
        assert_eq!(native.stdout, source.stdout, "{input} {flags:?}");
    }
}

#[test]
fn genomecov_filtered_group_order_and_raw_depth_boundaries_match_source() {
    for (input, flags) in [
        (
            "filtered-order.bed",
            vec!["--bedgraph-all", "--strand", "+"],
        ),
        ("tail.bed", vec!["--bedgraph", "--max", "1"]),
        ("touching.bed", vec!["--bedgraph-all"]),
    ] {
        let source = run(input, &flags, true);
        let native = run(input, &flags, false);
        assert!(source.status.success() && native.status.success());
        assert_eq!(native.stdout, source.stdout, "{input}");
    }
}

#[test]
fn genomecov_stdin_json_and_failed_late_group_preserve_named_output() {
    let directory = tempdir().unwrap();
    let output = directory.path().join("track.bg");
    let mut child = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["--json", "genomecov", "-i", "-", "-g"])
        .arg(fixture("genome.txt"))
        .args(["--bedgraph-all", "-o"])
        .arg(&output)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(fixture("ordinary.bed")).unwrap())
        .unwrap();
    let result = child.wait_with_output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let json: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(json["status"], "ok");
    assert_eq!(
        fs::read(&output).unwrap(),
        fs::read(fixture("bedgraph-all.stdout")).unwrap()
    );
    let before = fs::read(&output).unwrap();
    let invalid = directory.path().join("reentry.bed");
    fs::write(&invalid, b"chr1\t1\t3\nchr2\t2\t4\nchr1\t5\t6\n").unwrap();
    let failed = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["genomecov", "-i"])
        .arg(&invalid)
        .arg("-g")
        .arg(fixture("genome.txt"))
        .args(["--bedgraph-all", "-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(!failed.status.success());
    assert_eq!(fs::read(&output).unwrap(), before);
}

#[test]
fn genomecov_split_exact_width_and_point_geometry_match_source() {
    for input in ["point12.bed", "blocks-plus.bed"] {
        let flags = ["--bedgraph-all", "--split"];
        let source = run(input, &flags, true);
        let native = run(input, &flags, false);
        assert!(source.status.success());
        assert!(
            native.status.success(),
            "{input}: {}",
            String::from_utf8_lossy(&native.stderr)
        );
        assert_eq!(native.stdout, source.stdout, "{input}");
    }
}

#[test]
fn genomecov_split_rejects_selected_zero_origin_points() {
    let genome = read_genome(&b"chr1\t20\n"[..]).unwrap();
    assert!(
        genomecov(
            &b"chr1\t0\t0\tpoint\t0\t+\t0\t0\t0\t1\t0\t0\n"[..],
            &genome,
            Vec::new(),
            &GenomecovOptions {
                split: true,
                ..Default::default()
            }
        )
        .is_err()
    );
}

#[test]
fn genomecov_sparse_long_chromosome_keeps_runs_not_per_base_storage() {
    let genome = read_genome(&b"chr1\t1000000000000\n"[..]).unwrap();
    let mut output = Vec::new();
    genomecov(
        &b"chr1\t1\t3\n"[..],
        &genome,
        &mut output,
        &GenomecovOptions {
            mode: GenomecovMode::BedgraphAll,
            ..Default::default()
        },
    )
    .unwrap();
    assert_eq!(
        output,
        b"chr1\t0\t1\t0\nchr1\t1\t3\t1\nchr1\t3\t1000000000000\t0\n"
    );
}

#[test]
fn genomecov_rejects_invalid_coordinates_and_required_strands() {
    let genome = read_genome(&b"chr1\t20\n"[..]).unwrap();
    for input in ["unknown\t1\t2\n", "chr1\t1\t21\n", "chr1\t0\t0\n"] {
        assert!(
            genomecov(
                input.as_bytes(),
                &genome,
                Vec::new(),
                &GenomecovOptions::default()
            )
            .is_err()
        );
    }
    assert!(
        genomecov(
            &b"chr1\t1\t2\n"[..],
            &genome,
            Vec::new(),
            &GenomecovOptions {
                end: Some(rsomics_bed::genomecov::GenomecovEnd::FivePrime),
                ..Default::default()
            }
        )
        .is_err()
    );
}
