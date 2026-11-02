use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/map")
        .join(name)
}

const ALL: &str = "sum,mean,stdev,sstdev,median,mode,antimode,min,max,absmin,absmax,count,distinct,count_distinct,distinct_only,distinct_sort_num,distinct_sort_num_desc,collapse,concat,freqasc,freqdesc,first,last";

fn run(a: &str, b: &str, options: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["map", "-a"])
        .arg(fixture(a))
        .arg("-b")
        .arg(fixture(b))
        .args(options)
        .output()
        .unwrap()
}

#[test]
fn retained_complete_outputs_cover_selection_order_and_text() {
    let cases: &[(&str, &str, &str, &[&str])] = &[
        (
            "strings-a.bed",
            "strings-b.bed",
            "strings.stdout",
            &[
                "-c",
                "4",
                "--operations",
                "count,distinct,count_distinct,distinct_only,collapse,concat,mode,antimode,freqasc,freqdesc,first,last",
            ],
        ),
        (
            "a.bed",
            "b.bed",
            "encounter-order.stdout",
            &["-c", "4", "--operations", "collapse,first,last"],
        ),
        (
            "a.bed",
            "b.bed",
            "sorted-hits.stdout",
            &[
                "-c",
                "4",
                "--operations",
                "collapse,first,last",
                "--sortout",
            ],
        ),
        (
            "a.bed",
            "b.bed",
            "custom-controls.stdout",
            &[
                "-c",
                "4,5",
                "--operations",
                "distinct,mean",
                "--delimiter",
                "|",
                "--null",
                "NA",
                "--precision",
                "4",
            ],
        ),
        (
            "a.bed",
            "b.bed",
            "broadcast-selected.stdout",
            &[
                "-c",
                "4,5",
                "--operations",
                "count",
                "-s",
                "-f",
                "0.1",
                "-F",
                "0.25",
                "-e",
            ],
        ),
        (
            "split-a.bed",
            "b.bed",
            "ordinary-split.stdout",
            &["-c", "4,5", "--operations", "collapse,count", "--split"],
        ),
        (
            "header-a.bed",
            "header-b.bed",
            "header-positive.stdout",
            &[
                "-c",
                "4,5,5",
                "--operations",
                "collapse,count,sum",
                "--header",
            ],
        ),
        (
            "nonlex-a.bed",
            "nonlex-b.bed",
            "nonlex-aligned.stdout",
            &["-c", "4,5,5", "--operations", "collapse,count,sum"],
        ),
        (
            "lexical-a.bed",
            "lexical-b.bed",
            "lexical-b-exclusive.stdout",
            &["-c", "4,5,5", "--operations", "collapse,count,sum"],
        ),
        (
            "point-zero-a.bed",
            "point-zero-b.bed",
            "point-zero.stdout",
            &["-c", "4,5,5", "--operations", "collapse,count,sum"],
        ),
    ];
    for &(a, b, expected, options) in cases {
        let result = run(a, b, options);
        assert!(
            result.status.success(),
            "{expected}: {}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            result.stdout,
            fs::read(fixture(expected)).unwrap(),
            "{expected}"
        );
    }
    let genome = fixture("genome.txt");
    let result = run(
        "genome-a.bed",
        "genome-b.bed",
        &[
            "-c",
            "4,5",
            "--operations",
            "collapse,sum",
            "-g",
            genome.to_str().unwrap(),
        ],
    );
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        result.stdout,
        fs::read(fixture("genome-nonlex.stdout")).unwrap()
    );
}

#[test]
fn approved_point_loss_repairs_keep_both_equal_start_orders() {
    for b in ["ordinary-point-b.bed", "point-ordinary-b.bed"] {
        let result = run(
            "point-boundary-a.bed",
            b,
            &["-c", "4,5,5", "--operations", "collapse,count,sum"],
        );
        assert!(result.status.success());
        assert_eq!(result.stdout, b"chr1\t0\t10\tA\t0\t+\tpoint\t1\t7\n");
    }
    let result = run(
        "equal-start-a.bed",
        "equal-start-b.bed",
        &["-c", "4,5,5", "--operations", "collapse,count,sum"],
    );
    assert!(result.status.success());
    let source = fs::read_to_string(fixture("equal-start-corrected.stdout")).unwrap();
    let expected = source.replace(
        "point-after\t0\t+\tlong,right,nested\t3\t14",
        "point-after\t0\t+\tlong,left-boundary,right,nested\t4\t17",
    );
    assert_ne!(source, expected);
    assert_eq!(result.stdout, expected.as_bytes());
}

#[test]
fn invalid_selection_and_late_order_errors_preserve_named_output() {
    for options in [
        vec!["-c", "99", "--operations", "count"],
        vec!["--split", "-f", "0.5"],
        vec!["-F", "0.5", "-e"],
        vec!["-f", "0"],
    ] {
        assert!(!run("a.bed", "b.bed", &options).status.success());
    }
    let directory = tempfile::tempdir().unwrap();
    let output = directory.path().join("output.bed");
    for b in ["conflict-b.bed", "nonlex-exclusive-b.bed"] {
        fs::write(&output, b"existing result\n").unwrap();
        let result = run("nonlex-a.bed", b, &["--output", output.to_str().unwrap()]);
        assert!(!result.status.success());
        assert_eq!(fs::read(&output).unwrap(), b"existing result\n");
    }
}

#[test]
fn cli_aliases_json_and_lazy_numeric_errors_use_common_boundaries() {
    let directory = tempfile::tempdir().unwrap();
    let a = directory.path().join("a.bed");
    let b = directory.path().join("b.bed");
    let genome = directory.path().join("genome.txt");
    fs::write(&a, b"chr1\t0\t10\tA\n").unwrap();
    fs::write(&b, b"chr1\t1\t5\t\tinvalid\n").unwrap();
    fs::write(&genome, b"chr1\t100\n").unwrap();
    for input in [&a, &b, &genome] {
        let before = fs::read(input).unwrap();
        let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
            .args(["map", "-a"])
            .arg(&a)
            .arg("-b")
            .arg(&b)
            .arg("--genome")
            .arg(&genome)
            .arg("--output")
            .arg(input)
            .output()
            .unwrap();
        assert!(!result.status.success());
        assert_eq!(fs::read(input).unwrap(), before);
    }
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["--json", "map", "-a"])
        .arg(&a)
        .arg("-b")
        .arg(&b)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert!(result.stdout.is_empty());
    serde_json::from_slice::<serde_json::Value>(&result.stderr).unwrap();
    let output = directory.path().join("mapped.bed");
    let invoke = |options: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
            .args(["--json", "map", "-a"])
            .arg(&a)
            .arg("-b")
            .arg(&b)
            .args(options)
            .arg("--output")
            .arg(&output)
            .output()
            .unwrap()
    };
    let result = invoke(&["-c", "4,5", "--operations", "first,count"]);
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(fs::read(&output).unwrap(), b"chr1\t0\t10\tA\t\t1\n");
    let json: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(json["status"], "ok");
    assert!(result.stderr.is_empty());
    let before = fs::read(&output).unwrap();
    assert!(!invoke(&["-c", "5", "--operations", "sum"]).status.success());
    assert_eq!(fs::read(&output).unwrap(), before);
    fs::write(&b, b"# empty database\n").unwrap();
    assert!(!invoke(&["--operations", "count"]).status.success());
    assert_eq!(fs::read(&output).unwrap(), before);
}

#[test]
fn map_all23_reductions_match_complete_source_output() {
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["map", "-a"])
        .arg(fixture("a.bed"))
        .arg("-b")
        .arg(fixture("b.bed"))
        .args(["-c", "5", "--operations", ALL])
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    assert_eq!(
        result.stdout,
        fs::read(fixture("numeric-all23.stdout")).unwrap()
    );
}
