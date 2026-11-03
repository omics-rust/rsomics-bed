use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/coordinate")
        .join(name)
}

fn command(operation: &str) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"));
    command
        .arg(operation)
        .arg("--input")
        .arg(fixture("input.bed"))
        .arg("--genome")
        .arg(fixture("genome.tsv"))
        .arg("--header");
    command
}

#[test]
fn coordinate_workflows_match_pinned_upstream_and_preserve_transactional_output() {
    let cases: &[(&str, &[&str], &str)] = &[
        (
            "slop",
            &["--left", "3", "--right", "5", "--strand-relative"],
            "slop-absolute.stdout",
        ),
        (
            "slop",
            &[
                "--left",
                "0.5",
                "--right",
                "1.25",
                "--fraction",
                "--strand-relative",
            ],
            "slop-fraction.stdout",
        ),
        (
            "slop",
            &["--left", "-7", "--right", "-2", "--strand-relative"],
            "slop-shrink.stdout",
        ),
        ("shift", &["--shift", "8"], "shift-absolute.stdout"),
        (
            "shift",
            &["--plus", "0.5", "--minus", "-2.5", "--fraction"],
            "shift-fraction.stdout",
        ),
    ];
    for (operation, flags, expected) in cases {
        let output = command(operation).args(*flags).output().unwrap();
        assert!(
            output.status.success(),
            "{expected}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(output.stderr.is_empty());
        assert_eq!(
            output.stdout,
            fs::read(fixture(expected)).unwrap(),
            "{expected}"
        );
    }

    let directory = tempfile::tempdir().unwrap();
    let expanded = directory.path().join("expanded.bed");
    let output = command("slop")
        .args([
            "--left",
            "3",
            "--right",
            "5",
            "--strand-relative",
            "--json",
            "--output",
        ])
        .arg(&expanded)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(report["status"], "ok");
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["shift", "--shift", "0", "--header", "--input"])
        .arg(&expanded)
        .arg("--genome")
        .arg(fixture("genome.tsv"))
        .output()
        .unwrap();
    assert!(output.status.success());
    assert_eq!(
        output.stdout,
        fs::read(fixture("slop-absolute.stdout")).unwrap()
    );

    let incomplete = directory.path().join("incomplete.genome");
    fs::write(&incomplete, "chr1\t30\n").unwrap();
    let before = fs::read(&expanded).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["slop", "--bases", "1", "--input"])
        .arg(fixture("input.bed"))
        .arg("--genome")
        .arg(&incomplete)
        .arg("--output")
        .arg(&expanded)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert_eq!(fs::read(&expanded).unwrap(), before);

    let large_genome = directory.path().join("large.genome");
    let large_bed = directory.path().join("large.bed");
    fs::write(&large_genome, "chr1\t9007199254740995\n").unwrap();
    fs::write(&large_bed, "chr1\t0\t9007199254740995\n").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["shift", "--shift", "0", "--input"])
        .arg(&large_bed)
        .arg("--genome")
        .arg(&large_genome)
        .arg("--output")
        .arg(&expanded)
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(String::from_utf8_lossy(&output.stderr).contains("exceeds chromosome"));
    assert_eq!(fs::read(&expanded).unwrap(), before);
}
