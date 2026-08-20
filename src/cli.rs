//! Command-line parsing and dispatch for the product binary.

use std::path::{Path, PathBuf};
use std::process;

use clap::{Args, Parser, Subcommand, ValueEnum};
use rsomics_common::{
    OutputArgs, Result, RsomicsError, ToolMeta, reject_output_alias as reject_path_alias,
    run as run_tool, write_output,
};

use crate::io::open_input;
use crate::{closest, cluster, complement, intersect, merge, read_genome, sort, subtract, window};

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
    /// Assign cluster IDs to overlapping or nearby sorted intervals
    Cluster(ClusterArgs),
    /// Report B intervals in a configurable neighborhood of A
    Window(WindowArgs),
    /// Report the closest eligible B interval or intervals for each A
    Closest(ClosestArgs),
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

#[derive(Debug, Args)]
#[command(next_help_heading = "Input/output")]
struct ClusterArgs {
    /// Input sorted BED file; omit or use - for standard input
    #[arg(value_name = "BED")]
    input: Option<PathBuf>,

    /// Output BED file; omit or use - for standard output
    #[arg(short, long, value_name = "BED")]
    output: Option<PathBuf>,

    /// Maximum gap between records in one cluster
    #[arg(
        short = 'd',
        long,
        default_value_t = 0,
        value_name = "BP",
        help_heading = "Clustering"
    )]
    distance: u64,

    /// Strand policy
    #[arg(long, value_enum, value_name = "MODE", help_heading = "Clustering")]
    strand: Option<ClusterStrand>,

    /// Cluster forward and reverse strands independently
    #[arg(short = 's', conflicts_with = "strand", help_heading = "Clustering")]
    same_strand: bool,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum ClusterStrand {
    Any,
    Same,
}

#[derive(Debug, Args)]
#[command(next_help_heading = "Input/output")]
struct WindowArgs {
    /// BED file A; use - for standard input
    #[arg(short = 'a', long, value_name = "BED")]
    a: PathBuf,

    /// BED file B; standard input is not supported
    #[arg(short = 'b', long, value_name = "BED")]
    b: PathBuf,

    /// Output BED file; omit or use - for standard output
    #[arg(short, long, value_name = "BED")]
    output: Option<PathBuf>,

    /// Symmetric window size
    #[arg(
        short = 'w',
        long,
        value_name = "BP",
        conflicts_with_all = ["left", "right"],
        help_heading = "Window"
    )]
    window: Option<u64>,

    /// Bases added to the reference-left side of A
    #[arg(
        short = 'l',
        long,
        value_name = "BP",
        requires = "right",
        help_heading = "Window"
    )]
    left: Option<u64>,

    /// Bases added to the reference-right side of A
    #[arg(
        short = 'r',
        long,
        value_name = "BP",
        requires = "left",
        help_heading = "Window"
    )]
    right: Option<u64>,

    /// Interpret left and right relative to A's strand
    #[arg(long, help_heading = "Window")]
    strand_relative: bool,

    /// Strand relationship between A and B
    #[arg(
        long,
        value_enum,
        default_value = "any",
        value_name = "MODE",
        help_heading = "Window"
    )]
    strand: WindowStrand,

    /// Output reduction
    #[arg(long, value_enum, value_name = "MODE", help_heading = "Reporting")]
    report: Option<WindowReport>,

    /// Emit A once when any B record matches
    #[arg(
        short = 'u',
        conflicts_with_all = ["report", "report_count", "report_none"],
        help_heading = "Reporting"
    )]
    report_any: bool,

    /// Append the number of matching B records
    #[arg(
        short = 'c',
        conflicts_with_all = ["report", "report_any", "report_none"],
        help_heading = "Reporting"
    )]
    report_count: bool,

    /// Emit A only when no B record matches
    #[arg(
        short = 'v',
        conflicts_with_all = ["report", "report_any", "report_count"],
        help_heading = "Reporting"
    )]
    report_none: bool,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum WindowStrand {
    Any,
    Same,
    Opposite,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum WindowReport {
    Pairs,
    Any,
    Count,
    None,
}

#[derive(Debug, Args)]
#[command(next_help_heading = "Input/output")]
struct ClosestArgs {
    /// BED file A; use - for standard input
    #[arg(short = 'a', long, value_name = "BED")]
    a: PathBuf,

    /// BED file B; standard input is not supported
    #[arg(short = 'b', long, value_name = "BED")]
    b: PathBuf,

    /// Output BED file; omit or use - for standard output
    #[arg(short, long, value_name = "BED")]
    output: Option<PathBuf>,

    /// Strand relationship between A and B
    #[arg(long, value_enum, value_name = "MODE", help_heading = "Selection")]
    strand: Option<ClosestStrand>,

    /// Require matching strands
    #[arg(
        short = 's',
        conflicts_with_all = ["strand", "opposite_strand"],
        help_heading = "Selection"
    )]
    same_strand: bool,

    /// Require opposing strands
    #[arg(
        short = 'S',
        conflicts_with_all = ["strand", "same_strand"],
        help_heading = "Selection"
    )]
    opposite_strand: bool,

    /// Exclude B records with the same BED name as A
    #[arg(short = 'N', long, help_heading = "Selection")]
    different_name: bool,

    /// Exclude B records that overlap A
    #[arg(long, help_heading = "Selection")]
    ignore_overlaps: bool,

    /// Distance representation
    #[arg(long, value_enum, value_name = "MODE", help_heading = "Reporting")]
    distance: Option<ClosestDistance>,

    /// Append the non-negative distance
    #[arg(
        short = 'd',
        conflicts_with_all = ["distance", "signed_distance"],
        help_heading = "Reporting"
    )]
    unsigned_distance: bool,

    /// Append a signed distance
    #[arg(
        short = 'D',
        value_enum,
        value_name = "ORIENTATION",
        conflicts_with_all = ["distance", "unsigned_distance"],
        help_heading = "Reporting"
    )]
    signed_distance: Option<SignedDistance>,

    /// Equal-distance ties
    #[arg(
        short = 't',
        long,
        value_enum,
        default_value = "all",
        value_name = "MODE",
        help_heading = "Reporting"
    )]
    tie: ClosestTie,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum ClosestStrand {
    Any,
    Same,
    Opposite,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum ClosestDistance {
    None,
    Unsigned,
    Reference,
    A,
    B,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum SignedDistance {
    #[value(name = "ref")]
    Reference,
    A,
    B,
}

#[derive(Debug, Clone, Copy, ValueEnum)]
enum ClosestTie {
    All,
    First,
    Last,
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
        Command::Cluster(args) => {
            require_named_json_output(json, args.output.as_deref())?;
            reject_output_alias(args.output.as_deref(), [args.input.as_deref()])?;
            let input = open_input(args.input.as_deref())?;
            let same_strand = args.same_strand || matches!(args.strand, Some(ClusterStrand::Same));
            write_output(args.output.as_deref(), |output| {
                cluster::cluster(
                    input,
                    output,
                    cluster::ClusterOptions {
                        distance: args.distance,
                        same_strand,
                    },
                )
            })
        }
        Command::Window(args) => {
            require_named_json_output(json, args.output.as_deref())?;
            reject_stdin_b(&args.b)?;
            reject_output_alias(
                args.output.as_deref(),
                [Some(args.a.as_path()), Some(args.b.as_path())],
            )?;
            let (left, right) = match (args.window, args.left, args.right) {
                (Some(distance), None, None) => (distance, distance),
                (None, Some(left), Some(right)) => (left, right),
                (None, None, None) => (1000, 1000),
                _ => unreachable!("Clap validates window width combinations"),
            };
            let strand = match args.strand {
                WindowStrand::Any => crate::StrandFilter::Any,
                WindowStrand::Same => crate::StrandFilter::Same,
                WindowStrand::Opposite => crate::StrandFilter::Opposite,
            };
            let report = if args.report_any {
                window::WindowReport::Any
            } else if args.report_count {
                window::WindowReport::Count
            } else if args.report_none {
                window::WindowReport::None
            } else {
                match args.report.unwrap_or(WindowReport::Pairs) {
                    WindowReport::Pairs => window::WindowReport::Pairs,
                    WindowReport::Any => window::WindowReport::Any,
                    WindowReport::Count => window::WindowReport::Count,
                    WindowReport::None => window::WindowReport::None,
                }
            };
            let a = open_input(Some(&args.a))?;
            let b = open_input(Some(&args.b))?;
            write_output(args.output.as_deref(), |output| {
                window::window(
                    a,
                    b,
                    output,
                    window::WindowOptions {
                        left,
                        right,
                        strand_relative: args.strand_relative,
                        strand,
                        report,
                    },
                )
            })
        }
        Command::Closest(args) => {
            require_named_json_output(json, args.output.as_deref())?;
            reject_stdin_b(&args.b)?;
            reject_output_alias(
                args.output.as_deref(),
                [Some(args.a.as_path()), Some(args.b.as_path())],
            )?;
            let strand = if args.same_strand {
                crate::StrandFilter::Same
            } else if args.opposite_strand {
                crate::StrandFilter::Opposite
            } else {
                match args.strand.unwrap_or(ClosestStrand::Any) {
                    ClosestStrand::Any => crate::StrandFilter::Any,
                    ClosestStrand::Same => crate::StrandFilter::Same,
                    ClosestStrand::Opposite => crate::StrandFilter::Opposite,
                }
            };
            let distance = if args.unsigned_distance {
                closest::DistanceMode::Unsigned
            } else if let Some(mode) = args.signed_distance {
                match mode {
                    SignedDistance::Reference => closest::DistanceMode::Reference,
                    SignedDistance::A => closest::DistanceMode::A,
                    SignedDistance::B => closest::DistanceMode::B,
                }
            } else {
                match args.distance.unwrap_or(ClosestDistance::None) {
                    ClosestDistance::None => closest::DistanceMode::None,
                    ClosestDistance::Unsigned => closest::DistanceMode::Unsigned,
                    ClosestDistance::Reference => closest::DistanceMode::Reference,
                    ClosestDistance::A => closest::DistanceMode::A,
                    ClosestDistance::B => closest::DistanceMode::B,
                }
            };
            let tie = match args.tie {
                ClosestTie::All => closest::TieMode::All,
                ClosestTie::First => closest::TieMode::First,
                ClosestTie::Last => closest::TieMode::Last,
            };
            let a = open_input(Some(&args.a))?;
            let b = open_input(Some(&args.b))?;
            write_output(args.output.as_deref(), |output| {
                closest::closest(
                    a,
                    b,
                    output,
                    closest::ClosestOptions {
                        strand,
                        different_name: args.different_name,
                        ignore_overlaps: args.ignore_overlaps,
                        distance,
                        tie,
                    },
                )
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

        let error = Cli::try_parse_from(["rsomics-bed", "cluster", "--help"]).unwrap_err();
        let help = error.to_string();
        assert!(help.contains("Input/output:"), "{help}");
        assert!(help.contains("Clustering:"), "{help}");
        assert!(help.contains("--strand <MODE>"), "{help}");

        let error = Cli::try_parse_from(["rsomics-bed", "closest", "--help"]).unwrap_err();
        let help = error.to_string();
        assert!(help.contains("Input/output:"), "{help}");
        assert!(help.contains("Selection:"), "{help}");
        assert!(help.contains("Reporting:"), "{help}");
        assert!(help.contains("--distance <MODE>"), "{help}");
        assert!(help.contains("--ignore-overlaps"), "{help}");
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
