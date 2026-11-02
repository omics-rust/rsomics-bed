use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use tempfile::tempdir;

struct Case {
    label: &'static str,
    a: &'static str,
    bs: &'static [&'static str],
    flags: &'static [&'static str],
}

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/annotate")
        .join(name)
}

fn cases() -> [Case; 10] {
    [
        Case {
            label: "default",
            a: "a.bed",
            bs: &["b1.bed", "b2.bed"],
            flags: &[],
        },
        Case {
            label: "counts",
            a: "a.bed",
            bs: &["b1.bed", "b2.bed"],
            flags: &["--counts"],
        },
        Case {
            label: "both",
            a: "a.bed",
            bs: &["b1.bed", "b2.bed"],
            flags: &["--both"],
        },
        Case {
            label: "names",
            a: "a.bed",
            bs: &["b1.bed", "b2.bed"],
            flags: &["--names", "first", "second"],
        },
        Case {
            label: "names-both",
            a: "a.bed",
            bs: &["b1.bed", "b2.bed"],
            flags: &["--both", "--names", "first", "second"],
        },
        Case {
            label: "same",
            a: "a.bed",
            bs: &["b1.bed", "b2.bed"],
            flags: &["--both", "-s"],
        },
        Case {
            label: "opposite",
            a: "a.bed",
            bs: &["b1.bed", "b2.bed"],
            flags: &["--both", "-S"],
        },
        Case {
            label: "point-only",
            a: "point-a.bed",
            bs: &["point-b.bed"],
            flags: &["--both"],
        },
        Case {
            label: "empty-names",
            a: "empty.bed",
            bs: &["b1.bed", "b2.bed"],
            flags: &["--both", "--names", "first", "second"],
        },
        Case {
            label: "bed12-names",
            a: "a12.bed",
            bs: &["b1.bed", "b2.bed"],
            flags: &["--both", "--names", "first", "second"],
        },
    ]
}

fn run_case(case: &Case, native: bool) -> std::process::Output {
    let mut command = if native {
        Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
    } else {
        Command::new("bedtools")
    };
    command.args(["annotate", "-i"]).arg(fixture(case.a));
    command.arg(if native { "--files" } else { "-files" });
    command.args(case.bs.iter().map(|name| fixture(name)));
    for flag in case.flags {
        command.arg(if !native && flag.starts_with("--") {
            &flag[1..]
        } else {
            flag
        });
    }
    command.output().unwrap()
}

#[test]
fn annotate_ten_complete_source_witnesses_match() {
    for case in cases() {
        let result = run_case(&case, true);
        assert!(
            result.status.success(),
            "{}: {}",
            case.label,
            String::from_utf8_lossy(&result.stderr)
        );
        assert!(
            result.stderr.is_empty(),
            "{}",
            String::from_utf8_lossy(&result.stderr)
        );
        assert_eq!(
            result.stdout,
            fs::read(fixture(&format!("{}.stdout", case.label))).unwrap(),
            "{}",
            case.label
        );
    }
}

#[test]
fn annotate_ten_whole_outputs_match_live_bedtools_231() {
    let version = Command::new("bedtools").arg("--version").output().unwrap();
    assert!(version.status.success());
    assert_eq!(version.stdout, b"bedtools v2.31.1\n");
    for case in cases() {
        let upstream = run_case(&case, false);
        assert!(
            upstream.status.success(),
            "{}: {}",
            case.label,
            String::from_utf8_lossy(&upstream.stderr)
        );
        assert!(
            upstream.stderr.is_empty(),
            "{}",
            String::from_utf8_lossy(&upstream.stderr)
        );
        let native = run_case(&case, true);
        assert!(
            native.status.success(),
            "{}: {}",
            case.label,
            String::from_utf8_lossy(&native.stderr)
        );
        assert_eq!(native.stdout, upstream.stdout, "{}", case.label);
    }
}

#[test]
fn annotate_stdin_multifile_names_json_and_failure_preservation() {
    let directory = tempdir().unwrap();
    let output = directory.path().join("annotated.bed");
    let mut child = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["--json", "annotate", "-i", "-", "--files"])
        .arg(fixture("b1.bed"))
        .arg(fixture("b2.bed"))
        .args(["--both", "--names", "first", "second", "-o"])
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
        .write_all(&fs::read(fixture("a.bed")).unwrap())
        .unwrap();
    let result = child.wait_with_output().unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let envelope: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(envelope["status"], "ok");
    assert_eq!(
        fs::read(&output).unwrap(),
        fs::read(fixture("names-both.stdout")).unwrap()
    );

    let bad = directory.path().join("bad.bed");
    fs::write(&bad, b"chr1\t0\t10\nchr1\tbad\t20\n").unwrap();
    fs::write(&output, b"previous result\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["annotate", "-i"])
        .arg(fixture("a.bed"))
        .arg("--files")
        .arg(fixture("b1.bed"))
        .arg(&bad)
        .arg("-o")
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"previous result\n");

    let protected_b = directory.path().join("b.bed");
    fs::write(&protected_b, b"chr1\t0\t10\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["annotate", "-i"])
        .arg(fixture("a.bed"))
        .arg("--files")
        .arg(fixture("b1.bed"))
        .arg(&protected_b)
        .arg("-o")
        .arg(&protected_b)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(fs::read(&protected_b).unwrap(), b"chr1\t0\t10\n");

    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["annotate", "-i"])
        .arg(fixture("empty.bed"))
        .arg("--files")
        .arg(&bad)
        .args(["--names", "bad", "-o"])
        .arg(&output)
        .output()
        .unwrap();
    assert!(!result.status.success());
    assert_eq!(fs::read(&output).unwrap(), b"previous result\n");
}

#[test]
fn annotate_rejects_invalid_labels_modes_and_selected_fields() {
    for flags in [
        vec!["--names", "only-one"],
        vec!["--counts", "--both"],
        vec!["--split"],
        vec!["-f", "0.5"],
        vec!["--sorted"],
        vec!["--header"],
        vec!["-s", "-S"],
    ] {
        let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
            .args(["annotate", "-i"])
            .arg(fixture("a.bed"))
            .arg("--files")
            .arg(fixture("b1.bed"))
            .arg(fixture("b2.bed"))
            .args(&flags)
            .output()
            .unwrap();
        assert!(!result.status.success(), "unexpected success for {flags:?}");
    }
    let directory = tempdir().unwrap();
    let input = directory.path().join("input.bed");
    fs::write(&input, b"chr1\t10\t20\tunknown\t0\t.\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["annotate", "-i"])
        .arg(&input)
        .arg("--files")
        .arg(fixture("b1.bed"))
        .arg("-s")
        .output()
        .unwrap();
    assert!(!result.status.success());
    fs::write(&input, b"chr1\t0\t10\nchr1\t20\t30\textra\n").unwrap();
    let result = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["annotate", "-i"])
        .arg(&input)
        .arg("--files")
        .arg(fixture("b1.bed"))
        .output()
        .unwrap();
    assert!(!result.status.success());
}
