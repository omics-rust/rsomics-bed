use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use tempfile::tempdir;

fn binary() -> PathBuf {
    Path::new(env!("CARGO_BIN_EXE_rsomics-bed")).to_owned()
}

fn golden(name: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("tests/golden")
        .join(name)
}

#[test]
fn top_level_and_nested_help_render() {
    let top = Command::new(binary()).arg("--help").output().unwrap();
    assert!(top.status.success());
    let top = String::from_utf8(top.stdout).unwrap();
    for command in [
        "sort",
        "merge",
        "intersect",
        "subtract",
        "complement",
        "cluster",
        "window",
    ] {
        assert!(top.contains(command), "{top}");
    }
    assert!(top.contains("Global options:"), "{top}");
    assert!(top.contains("--json"), "{top}");
    for absent in ["--threads", "--seed", "--quiet", "--verbose"] {
        assert!(!top.contains(absent), "{top}");
    }

    let nested = Command::new(binary())
        .args(["intersect", "--help"])
        .output()
        .unwrap();
    assert!(nested.status.success());
    let nested = String::from_utf8(nested.stdout).unwrap();
    assert!(nested.contains("--a <BED>"), "{nested}");
    assert!(nested.contains("--b <BED>"), "{nested}");
    assert!(nested.contains("Input/output:"), "{nested}");
    assert!(nested.contains("Global options:"), "{nested}");
    assert!(!nested.contains("rsomics-bed-intersect"), "{nested}");
}

#[test]
fn product_binary_dispatches_window_count() {
    let output = Command::new(binary())
        .args(["window", "-a"])
        .arg(golden("window.a.bed"))
        .arg("-b")
        .arg(golden("window.b.bed"))
        .args(["--window", "5", "--report", "count"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        std::fs::read(golden("window.count.expected.bed")).unwrap()
    );
}

#[test]
fn window_rejects_conflicting_width_and_report_selectors() {
    let widths = Command::new(binary())
        .args(["window", "-a"])
        .arg(golden("window.a.bed"))
        .arg("-b")
        .arg(golden("window.b.bed"))
        .args(["--window", "5", "--left", "2", "--right", "3"])
        .output()
        .unwrap();
    assert_eq!(widths.status.code(), Some(2));

    let reports = Command::new(binary())
        .args(["window", "-a"])
        .arg(golden("window.a.bed"))
        .arg("-b")
        .arg(golden("window.b.bed"))
        .args(["--report", "count", "-u"])
        .output()
        .unwrap();
    assert_eq!(reports.status.code(), Some(2));
}

#[test]
fn product_binary_dispatches_cluster() {
    let output = Command::new(binary())
        .arg("cluster")
        .arg(golden("cluster.input.bed"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        std::fs::read(golden("cluster.default.expected.bed")).unwrap()
    );
}

#[test]
fn cluster_rejects_conflicting_strand_selectors() {
    let output = Command::new(binary())
        .args(["cluster", "--strand", "any", "-s"])
        .arg(golden("cluster.strand.input.bed"))
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(2));
    assert!(output.stdout.is_empty());
}

#[test]
fn product_binary_dispatches_sort() {
    let output = Command::new(binary())
        .arg("sort")
        .arg(golden("sort.input.bed"))
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    assert_eq!(
        output.stdout,
        std::fs::read(golden("sort.expected.bed")).unwrap()
    );
}

#[test]
fn invalid_input_exits_nonzero_with_context() {
    let output = Command::new(binary())
        .args(["intersect", "-a"])
        .arg(golden("intersect.huge-a.bed"))
        .arg("-b")
        .arg(golden("intersect.huge-b.bed"))
        .output()
        .unwrap();
    assert!(!output.status.success());
    assert!(output.stdout.is_empty());
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("backend maximum"), "{stderr}");
}

#[test]
fn named_output_is_transactional_on_failure() {
    let directory = tempdir().unwrap();
    let input = directory.path().join("invalid.bed");
    let output = directory.path().join("output.bed");
    fs::write(&input, b"chr1\t0\t1\nchr1\tbad\t2\n").unwrap();
    fs::write(&output, b"existing output\n").unwrap();

    let result = Command::new(binary())
        .arg("sort")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(1));
    assert_eq!(fs::read(output).unwrap(), b"existing output\n");
}

#[cfg(unix)]
#[test]
fn new_named_output_uses_normal_umask_permissions() {
    use std::os::unix::fs::PermissionsExt;

    let directory = tempdir().unwrap();
    let control = directory.path().join("normal-create-mode");
    let output = directory.path().join("output.bed");
    fs::write(&control, b"control").unwrap();

    let result = Command::new(binary())
        .arg("sort")
        .arg(golden("sort.input.bed"))
        .arg("-o")
        .arg(&output)
        .output()
        .unwrap();
    assert!(result.status.success());
    assert_eq!(
        fs::metadata(output).unwrap().permissions().mode() & 0o777,
        fs::metadata(control).unwrap().permissions().mode() & 0o777
    );
}

#[cfg(unix)]
#[test]
fn replacing_named_output_preserves_existing_permissions() {
    use std::os::unix::fs::PermissionsExt;

    let directory = tempdir().unwrap();
    let output = directory.path().join("output.bed");
    fs::write(&output, b"existing").unwrap();
    fs::set_permissions(&output, fs::Permissions::from_mode(0o666)).unwrap();

    let result = Command::new(binary())
        .arg("sort")
        .arg(golden("sort.input.bed"))
        .arg("-o")
        .arg(&output)
        .output()
        .unwrap();
    assert!(result.status.success());
    assert_eq!(
        fs::metadata(output).unwrap().permissions().mode() & 0o777,
        0o666
    );
}

#[test]
fn all_existing_output_alias_forms_are_rejected_without_data_loss() {
    let directory = tempdir().unwrap();
    let input = directory.path().join("input.bed");
    fs::write(&input, b"chr2\t2\t3\nchr1\t1\t2\n").unwrap();

    for output in [
        input.clone(),
        directory.path().join(".").join("input.bed"),
        {
            let path = directory.path().join("hardlink.bed");
            fs::hard_link(&input, &path).unwrap();
            path
        },
    ] {
        let result = Command::new(binary())
            .arg("sort")
            .arg(&input)
            .arg("-o")
            .arg(&output)
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(2));
        assert_eq!(fs::read(&input).unwrap(), b"chr2\t2\t3\nchr1\t1\t2\n");
    }
}

#[cfg(unix)]
#[test]
fn symlink_output_alias_is_rejected_without_following_and_truncating_target() {
    use std::os::unix::fs::symlink;

    let directory = tempdir().unwrap();
    let input = directory.path().join("input.bed");
    let output = directory.path().join("output.bed");
    fs::write(&input, b"chr2\t2\t3\nchr1\t1\t2\n").unwrap();
    symlink(&input, &output).unwrap();

    let result = Command::new(binary())
        .arg("sort")
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(2));
    assert_eq!(fs::read(input).unwrap(), b"chr2\t2\t3\nchr1\t1\t2\n");
}

#[test]
fn json_is_separate_from_named_bed_output() {
    let directory = tempdir().unwrap();
    let output_path = directory.path().join("sorted.bed");
    let result = Command::new(binary())
        .arg("--json")
        .arg("sort")
        .arg(golden("sort.input.bed"))
        .arg("-o")
        .arg(&output_path)
        .output()
        .unwrap();
    assert!(
        result.status.success(),
        "{}",
        String::from_utf8_lossy(&result.stderr)
    );
    let envelope: serde_json::Value = serde_json::from_slice(&result.stdout).unwrap();
    assert_eq!(envelope["status"], "ok");
    assert_eq!(envelope["tool"], "rsomics-bed");
    assert_eq!(
        fs::read(output_path).unwrap(),
        fs::read(golden("sort.expected.bed")).unwrap()
    );
}

#[test]
fn json_with_bed_stdout_fails_as_configuration_without_mixed_stdout() {
    let result = Command::new(binary())
        .arg("--json")
        .arg("sort")
        .arg(golden("sort.input.bed"))
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(2));
    assert!(result.stdout.is_empty());
    let first_line = result.stderr.split(|&byte| byte == b'\n').next().unwrap();
    let envelope: serde_json::Value = serde_json::from_slice(first_line).unwrap();
    assert_eq!(envelope["status"], "error");
    assert_eq!(envelope["error"]["kind"], "ConfigError");
}

#[test]
fn io_errors_use_the_stable_io_exit_code() {
    let missing = golden("does-not-exist.bed");
    let result = Command::new(binary())
        .arg("sort")
        .arg(missing)
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(4));
}

#[test]
fn complement_rejects_dual_stdin_before_reading() {
    let result = Command::new(binary())
        .args(["complement", "-", "-g", "-"])
        .stdin(Stdio::null())
        .output()
        .unwrap();
    assert_eq!(result.status.code(), Some(2));
    assert!(String::from_utf8_lossy(&result.stderr).contains("both read from standard input"));
}
