use std::path::PathBuf;

use clap::Args;
use rsomics_common::{Result, write_output};

use crate::coordinate::{self, ShiftOptions, SlopOptions};
use crate::{io::open_input, read_genome};

#[derive(Debug, Args)]
#[command(next_help_heading = "Input/output")]
struct InputArgs {
    /// Input BED file; use - for standard input
    #[arg(short, long, value_name = "BED")]
    input: PathBuf,
    /// Tab-delimited chromosome sizes; use - only with a named BED input
    #[arg(short, long, value_name = "GENOME")]
    genome: PathBuf,
    /// Output BED file; omit or use - for standard output
    #[arg(short, long, value_name = "BED")]
    output: Option<PathBuf>,
    /// Preserve initial BED header lines
    #[arg(long, help_heading = "Reporting")]
    header: bool,
}

#[derive(Debug, Args)]
pub(super) struct SlopArgs {
    #[command(flatten)]
    io: InputArgs,
    /// Expansion on both sides; negative values shrink intervals
    #[arg(short = 'b', long, value_name = "WIDTH", allow_hyphen_values = true, required_unless_present = "left", conflicts_with_all = ["left", "right"], help_heading = "Coordinates")]
    bases: Option<f32>,
    /// Left expansion; requires --right
    #[arg(
        short = 'l',
        long,
        value_name = "WIDTH",
        allow_hyphen_values = true,
        requires = "right",
        conflicts_with = "bases",
        help_heading = "Coordinates"
    )]
    left: Option<f32>,
    /// Right expansion; requires --left
    #[arg(
        short = 'r',
        long,
        value_name = "WIDTH",
        allow_hyphen_values = true,
        requires = "left",
        conflicts_with = "bases",
        help_heading = "Coordinates"
    )]
    right: Option<f32>,
    /// Interpret widths as fractions of each interval's length
    #[arg(long, help_heading = "Coordinates")]
    fraction: bool,
    /// Swap left and right widths for reverse-strand records
    #[arg(long, help_heading = "Coordinates")]
    strand_relative: bool,
}

#[derive(Debug, Args)]
pub(super) struct ShiftArgs {
    #[command(flatten)]
    io: InputArgs,
    /// Signed shift applied to every interval
    #[arg(long, value_name = "DISTANCE", allow_hyphen_values = true, required_unless_present = "plus", conflicts_with_all = ["plus", "minus"], help_heading = "Coordinates")]
    shift: Option<f32>,
    /// Shift for forward or unknown strands; requires --minus
    #[arg(
        long,
        value_name = "DISTANCE",
        allow_hyphen_values = true,
        requires = "minus",
        conflicts_with = "shift",
        help_heading = "Coordinates"
    )]
    plus: Option<f32>,
    /// Shift for reverse-strand records; requires --plus
    #[arg(
        long,
        value_name = "DISTANCE",
        allow_hyphen_values = true,
        requires = "plus",
        conflicts_with = "shift",
        help_heading = "Coordinates"
    )]
    minus: Option<f32>,
    /// Interpret shifts as fractions of each interval's length
    #[arg(long, help_heading = "Coordinates")]
    fraction: bool,
}

impl InputArgs {
    fn validate(&self, json: bool) -> Result<()> {
        super::require_named_json_output(json, self.output.as_deref())?;
        super::reject_dual_stdin(Some(&self.input), &self.genome)?;
        super::reject_output_alias(
            self.output.as_deref(),
            [Some(self.input.as_path()), Some(self.genome.as_path())],
        )
    }
}

pub(super) fn execute_slop(args: SlopArgs, json: bool) -> Result<()> {
    args.io.validate(json)?;
    let (left, right) = match (args.bases, args.left, args.right) {
        (Some(width), None, None) => (width, width),
        (None, Some(left), Some(right)) => (left, right),
        _ => unreachable!("Clap validates slop width combinations"),
    };
    let options = SlopOptions {
        left,
        right,
        fraction: args.fraction,
        strand_relative: args.strand_relative,
        header: args.io.header,
    };
    let genome = read_genome(open_input(Some(&args.io.genome))?)?;
    let input = open_input(Some(&args.io.input))?;
    write_output(args.io.output.as_deref(), |output| {
        coordinate::slop(input, &genome, output, &options)
    })
}

pub(super) fn execute_shift(args: ShiftArgs, json: bool) -> Result<()> {
    args.io.validate(json)?;
    let (plus, minus) = match (args.shift, args.plus, args.minus) {
        (Some(distance), None, None) => (distance, distance),
        (None, Some(plus), Some(minus)) => (plus, minus),
        _ => unreachable!("Clap validates shift distance combinations"),
    };
    let options = ShiftOptions {
        plus,
        minus,
        fraction: args.fraction,
        header: args.io.header,
    };
    let genome = read_genome(open_input(Some(&args.io.genome))?)?;
    let input = open_input(Some(&args.io.input))?;
    write_output(args.io.output.as_deref(), |output| {
        coordinate::shift(input, &genome, output, &options)
    })
}
