use std::fs;
use std::io::Write;
use std::process::{Command, Stdio};

use tempfile::tempdir;

const A: &[u8] = b"chr1\t0\t10\tA\t0\t+\nchr2\t0\t3\tnohit\t0\t+\nchr1\t100\t100\tzero\t0\t+\n";
const B: &[u8] = b"chr1\t0\t5\tb1\t2\t+\nchr1\t3\t8\tb2\t3\t-\nchr1\t3\t8\tb2dup\t3\t-\nchr1\t4\t6\tb3\t4\t+\nchr1\t100\t100\tzeroB\t1\t+\n";

fn run(a: &[u8], b: &[u8], flags: &[&str]) -> std::process::Output {
    let directory = tempdir().unwrap();
    let b_path = directory.path().join("b.bed");
    fs::write(&b_path, b).unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["coverage", "-a", "-", "-b"])
        .arg(b_path)
        .args(flags)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(a).unwrap();
    child.wait_with_output().unwrap()
}

#[test]
fn coverage_preserves_multiplicity_breadth_and_a_order() {
    let output = run(A, B, &[]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(output.stdout, b"chr1\t0\t10\tA\t0\t+\t4\t8\t10\t0.8000000\nchr2\t0\t3\tnohit\t0\t+\t0\t0\t3\t0.0000000\nchr1\t100\t100\tzero\t0\t+\t1\t2\t2\t1.0000000\n");
}

#[test]
fn coverage_modes_report_depth_and_histogram_totals() {
    for (flags, expected) in [
        (
            vec!["--counts"],
            "chr1\t0\t10\tA\t0\t+\t4\nchr2\t0\t3\tnohit\t0\t+\t0\nchr1\t100\t100\tzero\t0\t+\t1\n",
        ),
        (
            vec!["--mean"],
            "chr1\t0\t10\tA\t0\t+\t1.7000000\nchr2\t0\t3\tnohit\t0\t+\t0.0000000\nchr1\t100\t100\tzero\t0\t+\t1.0000000\n",
        ),
        (
            vec!["--hist"],
            "chr1\t0\t10\tA\t0\t+\t0\t2\t10\t0.2000000\nchr1\t0\t10\tA\t0\t+\t1\t3\t10\t0.3000000\nchr1\t0\t10\tA\t0\t+\t2\t2\t10\t0.2000000\nchr1\t0\t10\tA\t0\t+\t3\t2\t10\t0.2000000\nchr1\t0\t10\tA\t0\t+\t4\t1\t10\t0.1000000\nchr2\t0\t3\tnohit\t0\t+\t0\t3\t3\t1.0000000\nchr1\t100\t100\tzero\t0\t+\t1\t2\t2\t1.0000000\nall\t0\t5\t15\t0.3333333\nall\t1\t5\t15\t0.3333333\nall\t2\t2\t15\t0.1333333\nall\t3\t2\t15\t0.1333333\nall\t4\t1\t15\t0.0666667\n",
        ),
    ] {
        let output = run(A, B, &flags);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert_eq!(output.stdout, expected.as_bytes());
    }
    let output = run(b"chr1\t0\t10\tA\t0\t+\n", B, &["-d"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let expected = [1, 1, 1, 3, 4, 3, 2, 2, 0, 0]
        .iter()
        .enumerate()
        .map(|(i, depth)| format!("chr1\t0\t10\tA\t0\t+\t{}\t{depth}\n", i + 1))
        .collect::<String>();
    assert_eq!(output.stdout, expected.as_bytes());
}

#[test]
fn coverage_filters_candidates_before_reduction() {
    for (flags, count, covered) in [
        (vec!["-s"], 2, 6),
        (vec!["-S"], 2, 5),
        (vec!["-f", "0.5"], 3, 8),
        (vec!["-f", "0.5", "-r"], 3, 8),
        (vec!["-f", "0.9", "-F", "1", "-e"], 4, 8),
        (vec!["-f", "0.9", "-F", "1"], 0, 0),
    ] {
        let output = run(b"chr1\t0\t10\tA\t0\t+\n", B, &flags);
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let fraction = match covered {
            6 => "0.6000000",
            5 => "0.5000000",
            8 => "0.8000000",
            0 => "0.0000000",
            _ => unreachable!(),
        };
        assert_eq!(
            output.stdout,
            format!("chr1\t0\t10\tA\t0\t+\t{count}\t{covered}\t10\t{fraction}\n").as_bytes()
        );
    }
}

#[test]
fn coverage_split_retains_outer_span_and_counts_block_overlaps() {
    let a = b"chr1\t0\t10\tA\t0\t+\t0\t10\t0\t2\t2,2\t0,8\n";
    let b = b"chr1\t0\t1\tone\t2\t+\nchr1\t8\t10\ttwo\t3\t+\n";
    let output = run(a, b, &["--split"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"chr1\t0\t10\tA\t0\t+\t0\t10\t0\t2\t2,2\t0,8\t2\t3\t10\t0.3000000\n"
    );
    let output = run(a, a, &["--split"]);
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        b"chr1\t0\t10\tA\t0\t+\t0\t10\t0\t2\t2,2\t0,8\t2\t4\t10\t0.4000000\n"
    );
}

#[test]
fn coverage_pooled_b_json_and_transactional_failures() {
    let directory = tempdir().unwrap();
    let a = directory.path().join("a.bed");
    let b = directory.path().join("b.bed");
    let b2 = directory.path().join("b2.bed");
    let output = directory.path().join("out.bed");
    fs::write(&a, A).unwrap();
    fs::write(&b, B).unwrap();
    fs::write(&b2, B).unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["--json", "coverage", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .arg(&b2)
        .arg("-o")
        .arg(&output)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let json: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert!(json.is_object());
    assert_eq!(fs::read(&output).unwrap(), b"chr1\t0\t10\tA\t0\t+\t8\t8\t10\t0.8000000\nchr2\t0\t3\tnohit\t0\t+\t0\t0\t3\t0.0000000\nchr1\t100\t100\tzero\t0\t+\t2\t2\t2\t1.0000000\n");

    fs::write(&output, b"preserved\n").unwrap();
    fs::write(&a, b"chr1\t0\t10\nchr1\tbad\t20\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["coverage", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .arg("-o")
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"preserved\n");

    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["coverage", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .arg(&b2)
        .arg("-o")
        .arg(&b2)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(fs::read(&b2).unwrap(), B);
}

#[test]
fn coverage_rejects_unsupported_and_invalid_selections() {
    for flags in [
        vec!["--split", "-f", "0.5"],
        vec!["-f", "NaN"],
        vec!["-f", "0"],
        vec!["-f", "1.1"],
        vec!["-F", "0.5", "-r"],
        vec!["-F", "0.5", "-e"],
        vec!["--counts", "--mean"],
        vec!["--sorted"],
        vec!["--header"],
    ] {
        let result = run(A, B, &flags);
        assert!(!result.status.success(), "unexpected success for {flags:?}");
    }
    let result = run(
        b"chr1\t0\t10\tA\t0\t+\t0\t10\t0\t2\t2,2\t0\n",
        B,
        &["--split"],
    );
    assert!(!result.status.success());
    let result = run(A, b"chr1\t0\t10\tB\t0\t.\n", &["-s"]);
    assert!(!result.status.success());
}

#[test]
fn coverage_large_sparse_span_does_not_require_per_base_storage() {
    let result = run(b"chr1\t0\t2000000000\n", b"chr1\t10\t20\n", &[]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        result.stdout,
        b"chr1\t0\t2000000000\t1\t10\t2000000000\t0.0000000\n"
    );
}

#[test]
fn coverage_whole_outputs_match_bedtools_231() {
    let version = Command::new("bedtools").arg("--version").output().unwrap();
    assert!(version.status.success());
    assert_eq!(version.stdout, b"bedtools v2.31.1\n");
    let directory = tempdir().unwrap();
    let a = directory.path().join("a.bed");
    let b = directory.path().join("b.bed");
    fs::write(&a, A).unwrap();
    fs::write(&b, B).unwrap();
    for (ours, theirs) in [
        (vec![], vec![]),
        (vec!["--counts"], vec!["-counts"]),
        (vec!["--mean"], vec!["-mean"]),
        (vec!["--hist"], vec!["-hist"]),
        (vec!["-d"], vec!["-d"]),
        (vec!["-s"], vec!["-s"]),
        (vec!["-S"], vec!["-S"]),
        (vec!["-f", "0.5"], vec!["-f", "0.5"]),
        (vec!["-f", "0.5", "-r"], vec!["-f", "0.5", "-r"]),
        (
            vec!["-f", "0.9", "-F", "1", "-e"],
            vec!["-f", "0.9", "-F", "1", "-e"],
        ),
        (vec!["-f", "0.9", "-F", "1"], vec!["-f", "0.9", "-F", "1"]),
    ] {
        let upstream = Command::new("bedtools")
            .args(["coverage", "-a"])
            .arg(&a)
            .arg("-b")
            .arg(&b)
            .args(theirs)
            .output()
            .unwrap();
        assert!(
            upstream.status.success(),
            "{}",
            String::from_utf8_lossy(&upstream.stderr)
        );
        let native = run(A, B, &ours);
        assert!(
            native.status.success(),
            "{}",
            String::from_utf8_lossy(&native.stderr)
        );
        assert_eq!(native.stdout, upstream.stdout, "coverage {ours:?}");
    }
    let split_a = b"chr1\t0\t10\tA\t0\t+\t0\t10\t0\t2\t2,2\t0,8\n";
    let split_b = b"chr1\t0\t1\tone\t2\t+\nchr1\t8\t10\ttwo\t3\t+\n";
    fs::write(&a, split_a).unwrap();
    for b_input in [split_b.as_slice(), split_a.as_slice()] {
        fs::write(&b, b_input).unwrap();
        for (ours, theirs) in [
            (vec!["--split"], vec!["-split"]),
            (vec!["--split", "--hist"], vec!["-split", "-hist"]),
            (vec!["--split", "-d"], vec!["-split", "-d"]),
        ] {
            let upstream = Command::new("bedtools")
                .args(["coverage", "-a"])
                .arg(&a)
                .arg("-b")
                .arg(&b)
                .args(theirs)
                .output()
                .unwrap();
            assert!(
                upstream.status.success(),
                "{}",
                String::from_utf8_lossy(&upstream.stderr)
            );
            let native = run(split_a, b_input, &ours);
            assert!(
                native.status.success(),
                "{}",
                String::from_utf8_lossy(&native.stderr)
            );
            assert_eq!(native.stdout, upstream.stdout, "coverage {ours:?}");
        }
    }
}
