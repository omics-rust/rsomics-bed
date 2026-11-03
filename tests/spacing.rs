use std::fs;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

fn fixture(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden/spacing")
        .join(name)
}

fn command(input: &str) -> Command {
    let mut command = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"));
    command.args(["spacing", "--header", "--input"]);
    command.arg(fixture(input));
    command
}

#[test]
fn spacing_whole_stream_matches_upstream_and_preserves_failed_output() {
    for genome in [None, Some("order-good.genome")] {
        let mut command = command("integrated.bed");
        if let Some(genome) = genome {
            command.arg("--genome").arg(fixture(genome));
        }
        let output = command.output().unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(output.stderr.is_empty());
        assert_eq!(
            output.stdout,
            fs::read(fixture("integrated.stdout")).unwrap()
        );
    }
    let output = command("header-only.bed").output().unwrap();
    assert!(output.status.success());
    assert!(output.stderr.is_empty());
    assert_eq!(output.stdout, fs::read(fixture("header-only.bed")).unwrap());

    let directory = tempfile::tempdir().unwrap();
    let destination = directory.path().join("spacing.bed");
    let mut child = Command::new(env!("CARGO_BIN_EXE_rsomics-bed"))
        .args(["--json", "spacing", "--header", "--input", "-", "--genome"])
        .arg(fixture("order-good.genome"))
        .arg("--output")
        .arg(&destination)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(&fs::read(fixture("integrated.bed")).unwrap())
        .unwrap();
    let output = child.wait_with_output().unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert!(output.stderr.is_empty());
    let envelope: serde_json::Value = serde_json::from_slice(&output.stdout).unwrap();
    assert_eq!(envelope["status"], "ok");
    assert_eq!(
        fs::read(&destination).unwrap(),
        fs::read(fixture("integrated.stdout")).unwrap()
    );

    let before = fs::read(&destination).unwrap();
    for (input, genome) in [
        ("descending.bed", None),
        ("reappearance.bed", None),
        ("integrated.bed", Some("order-bad.genome")),
    ] {
        let mut command = command(input);
        command.arg("--json").arg("--output").arg(&destination);
        if let Some(genome) = genome {
            command.arg("--genome").arg(fixture(genome));
        }
        let output = command.output().unwrap();
        assert!(!output.status.success());
        assert!(output.stdout.is_empty());
        let envelope: serde_json::Value = serde_json::from_slice(&output.stderr).unwrap();
        assert_eq!(envelope["error"]["kind"], "InvalidInput");
        assert_eq!(fs::read(&destination).unwrap(), before);
    }
}
