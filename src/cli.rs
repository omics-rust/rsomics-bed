//! Command-line parsing and dispatch for the product binary.

use std::path::{Path, PathBuf};
use std::process;

use clap::{Args, Parser, Subcommand};
use rsomics_common::{
    OutputArgs, Result, RsomicsError, ToolMeta, reject_output_alias as reject_path_alias,
    run as run_tool, write_output,
};

use crate::io::open_input;
use crate::{complement, intersect, merge, read_genome, sort, subtract};

const META: ToolMeta = ToolMeta {
    name: "rsomics-bed",
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Debug, Parser)]
#[command(
    name = "rsomics-bed",
    version,
    about = "High-performance BED interval operations",
    arg_required_else_help = true
)]
struct Cli {
    #[command(flatten)]
    output: OutputArgs,

    #[command(subcommand)]
    command: Command,
}

#[derive(Debug, Subcommand)]
enum Command {
    /// Stable coordinate sort preserving all BED columns
    Sort(UnaryArgs),
    /// Merge overlapping or touching intervals from sorted BED input
    Merge(UnaryArgs),
    /// Emit clipped A intervals overlapping B
    Intersect(BinaryArgs),
    /// Remove B coverage from A
    Subtract(BinaryArgs),
    /// Emit regions not covered by sorted BED input
    Complement(ComplementArgs),
}

#[derive(Debug, Args)]
#[command(next_help_heading = "Input/output")]
struct UnaryArgs {
    /// Input BED file; omit or use - for standard input
    #[arg(value_name = "BED")]
    input: Option<PathBuf>,

    /// Output BED file; omit or use - for standard output
    #[arg(short, long, value_name = "BED")]
    output: Option<PathBuf>,
}

#[derive(Debug, Args)]
#[command(next_help_heading = "Input/output")]
struct BinaryArgs {
    /// BED file A; use - for standard input
    #[arg(short = 'a', long, value_name = "BED")]
    a: PathBuf,

    /// BED file B; standard input is not supported
    #[arg(short = 'b', long, value_name = "BED")]
    b: PathBuf,

    /// Output BED file; omit or use - for standard output
    #[arg(short, long, value_name = "BED")]
    output: Option<PathBuf>,
}

#[derive(Debug, Args)]
#[command(next_help_heading = "Input/output")]
struct ComplementArgs {
    /// Input sorted BED file; omit or use - for standard input
    #[arg(value_name = "BED")]
    input: Option<PathBuf>,

    /// Two-column chromosome-size file; use - only with a named BED input
    #[arg(short = 'g', long, value_name = "TSV")]
    genome: PathBuf,

    /// Output BED file; omit or use - for standard output
    #[arg(short, long, value_name = "BED")]
    output: Option<PathBuf>,
}

#[must_use]
pub(crate) fn run() -> process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    run_tool(&output, META, || execute(cli))
}

fn execute(cli: Cli) -> Result<()> {
    let json = cli.output.json;
    match cli.command {
        Command::Sort(args) => {
            require_named_json_output(json, args.output.as_deref())?;
            reject_output_alias(args.output.as_deref(), [args.input.as_deref()])?;
            let input = open_input(args.input.as_deref())?;
            write_output(args.output.as_deref(), |output| sort::sort(input, output))
        }
        Command::Merge(args) => {
            require_named_json_output(json, args.output.as_deref())?;
            reject_output_alias(args.output.as_deref(), [args.input.as_deref()])?;
            let input = open_input(args.input.as_deref())?;
            write_output(args.output.as_deref(), |output| merge::merge(input, output))
        }
        Command::Intersect(args) => {
            require_named_json_output(json, args.output.as_deref())?;
            reject_stdin_b(&args.b)?;
            reject_output_alias(
                args.output.as_deref(),
                [Some(args.a.as_path()), Some(args.b.as_path())],
            )?;
            let a = open_input(Some(&args.a))?;
            let b = open_input(Some(&args.b))?;
            write_output(args.output.as_deref(), |output| {
                intersect::intersect(a, b, output)
            })
        }
        Command::Subtract(args) => {
            require_named_json_output(json, args.output.as_deref())?;
            reject_stdin_b(&args.b)?;
            reject_output_alias(
                args.output.as_deref(),
                [Some(args.a.as_path()), Some(args.b.as_path())],
            )?;
            let a = open_input(Some(&args.a))?;
            let b = open_input(Some(&args.b))?;
            write_output(args.output.as_deref(), |output| {
                subtract::subtract(a, b, output)
            })
        }
        Command::Complement(args) => {
            require_named_json_output(json, args.output.as_deref())?;
            reject_dual_stdin(args.input.as_deref(), &args.genome)?;
            reject_output_alias(
                args.output.as_deref(),
                [args.input.as_deref(), Some(args.genome.as_path())],
            )?;
            let input = open_input(args.input.as_deref())?;
            let genome_input = open_input(Some(&args.genome))?;
            let genome = read_genome(genome_input)?;
            write_output(args.output.as_deref(), |output| {
                complement::complement(input, &genome, output)
            })
        }
    }
}

fn require_named_json_output(json: bool, output: Option<&Path>) -> Result<()> {
    if json && output.is_none_or(|path| path == Path::new("-")) {
        return Err(RsomicsError::ConfigError(
            "--json requires a named --output file so JSON cannot mix with BED stdout".to_owned(),
        ));
    }
    Ok(())
}

fn reject_stdin_b(path: &Path) -> Result<()> {
    if path == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "B must be a file because A is the streaming input".to_owned(),
        ));
    }
    Ok(())
}

fn reject_dual_stdin(input: Option<&Path>, genome: &Path) -> Result<()> {
    let bed_is_stdin = input.is_none_or(|path| path == Path::new("-"));
    if bed_is_stdin && genome == Path::new("-") {
        return Err(RsomicsError::ConfigError(
            "BED input and genome cannot both read from standard input".to_owned(),
        ));
    }
    Ok(())
}

fn reject_output_alias<'a>(
    output: Option<&Path>,
    inputs: impl IntoIterator<Item = Option<&'a Path>>,
) -> Result<()> {
    let Some(output) = output else {
        return Ok(());
    };
    reject_path_alias(output, inputs.into_iter().flatten())
}

#[cfg(test)]
mod tests {
    use std::fs;

    use clap::CommandFactory;
    use tempfile::tempdir;

    use super::*;

    #[test]
    fn command_tree_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn nested_help_is_owned_by_clap() {
        let error = Cli::try_parse_from(["rsomics-bed", "intersect", "--help"]).unwrap_err();
        assert_eq!(error.kind(), clap::error::ErrorKind::DisplayHelp);
        let help = error.to_string();
        assert!(help.contains("--a <BED>"), "{help}");
        assert!(help.contains("--b <BED>"), "{help}");
    }

    #[test]
    fn global_options_are_limited_to_shared_output() {
        let error = Cli::try_parse_from(["rsomics-bed", "--help"]).unwrap_err();
        let help = error.to_string();
        assert!(help.contains("Global options:"), "{help}");
        assert!(help.contains("--json"), "{help}");
        for absent in ["--threads", "--seed", "--quiet", "--verbose"] {
            assert!(!help.contains(absent), "{help}");
        }
    }

    #[test]
    fn output_cannot_alias_input() {
        let path = Path::new("same.bed");
        let error = reject_output_alias(Some(path), [Some(path)]).unwrap_err();
        assert!(error.to_string().contains("also an input"), "{error}");
    }

    #[test]
    fn output_cannot_alias_input_through_path_normalization() {
        let directory = tempdir().unwrap();
        let input = directory.path().join("input.bed");
        fs::write(&input, b"chr1\t0\t1\n").unwrap();
        let output = directory.path().join(".").join("input.bed");
        let error = reject_output_alias(Some(&output), [Some(input.as_path())]).unwrap_err();
        assert!(error.to_string().contains("also an input"), "{error}");
    }

    #[test]
    fn output_cannot_alias_input_through_hard_link() {
        let directory = tempdir().unwrap();
        let input = directory.path().join("input.bed");
        let output = directory.path().join("output.bed");
        fs::write(&input, b"chr1\t0\t1\n").unwrap();
        fs::hard_link(&input, &output).unwrap();
        let error = reject_output_alias(Some(&output), [Some(input.as_path())]).unwrap_err();
        assert!(error.to_string().contains("also an input"), "{error}");
        assert_eq!(fs::read(input).unwrap(), b"chr1\t0\t1\n");
    }

    #[test]
    fn json_requires_named_bed_output() {
        let error = require_named_json_output(true, None).unwrap_err();
        assert!(error.to_string().contains("requires a named"), "{error}");
        let error = require_named_json_output(true, Some(Path::new("-"))).unwrap_err();
        assert!(error.to_string().contains("requires a named"), "{error}");
        require_named_json_output(true, Some(Path::new("result.bed"))).unwrap();
    }

    #[test]
    fn complement_rejects_two_stdin_consumers() {
        let error = reject_dual_stdin(None, Path::new("-")).unwrap_err();
        assert!(error.to_string().contains("both read"), "{error}");
        reject_dual_stdin(Some(Path::new("input.bed")), Path::new("-")).unwrap();
    }
}
