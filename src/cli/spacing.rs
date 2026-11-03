use std::path::PathBuf;

use clap::Args as ClapArgs;
use rsomics_common::{Result, write_output};

use crate::{
    io::open_input,
    read_genome,
    spacing::{self, SpacingOptions},
};

#[derive(Debug, ClapArgs)]
#[command(next_help_heading = "Input/output")]
pub(super) struct Args {
    /// Sorted BED file; omit or use - for standard input
    #[arg(short, long, value_name = "BED")]
    input: Option<PathBuf>,
    /// Two-column chromosome sizes defining input order
    #[arg(short, long, value_name = "GENOME")]
    genome: Option<PathBuf>,
    /// Output BED file; omit or use - for standard output
    #[arg(short, long, value_name = "BED")]
    output: Option<PathBuf>,
    /// Preserve initial BED header lines
    #[arg(long, help_heading = "Reporting")]
    header: bool,
}

pub(super) fn execute(args: Args, json: bool) -> Result<()> {
    super::require_named_json_output(json, args.output.as_deref())?;
    if let Some(genome) = &args.genome {
        super::reject_dual_stdin(args.input.as_deref(), genome)?;
    }
    super::reject_output_alias(
        args.output.as_deref(),
        [args.input.as_deref(), args.genome.as_deref()],
    )?;
    let genome = args
        .genome
        .as_deref()
        .map(|path| read_genome(open_input(Some(path))?))
        .transpose()?;
    let input = open_input(args.input.as_deref())?;
    write_output(args.output.as_deref(), |output| {
        spacing::spacing(
            input,
            genome.as_ref(),
            output,
            SpacingOptions {
                header: args.header,
            },
        )
    })
}
