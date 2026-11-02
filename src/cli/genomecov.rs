use std::path::PathBuf;

use clap::Args as ClapArgs;
use rsomics_common::{Result, write_output};

use crate::genomecov::{self, GenomecovEnd, GenomecovMode, GenomecovOptions, GenomecovStrand};
use crate::{io::open_input, read_genome};

#[derive(Debug, ClapArgs)]
#[command(next_help_heading = "Input/output")]
pub(super) struct Args {
    /// Chromosome-grouped BED file; use - for standard input
    #[arg(short, long, value_name = "BED")]
    input: PathBuf,
    /// Tab-delimited chromosome sizes
    #[arg(short, long, value_name = "GENOME")]
    genome: PathBuf,
    /// Output file; omit or use - for standard output
    #[arg(short, long, value_name = "OUTPUT")]
    output: Option<PathBuf>,
    /// Emit nonzero bedGraph depth runs
    #[arg(long, group = "mode", help_heading = "Reporting")]
    bedgraph: bool,
    /// Emit bedGraph runs including uncovered regions
    #[arg(long, group = "mode", help_heading = "Reporting")]
    bedgraph_all: bool,
    /// Emit depth at every base with one-based positions
    #[arg(long, group = "mode", help_heading = "Reporting")]
    depth: bool,
    /// Emit only covered bases with zero-based positions
    #[arg(long, group = "mode", help_heading = "Reporting")]
    depth_nonzero: bool,
    /// Maximum histogram bin; also caps nonfinal bedGraph runs
    #[arg(long, help_heading = "Reporting")]
    max: Option<u32>,
    /// Finite depth multiplier; not valid with histogram output
    #[arg(long, allow_hyphen_values = true, help_heading = "Reporting")]
    scale: Option<f32>,
    /// Browser track options; requires bedGraph output
    #[arg(long, value_name = "OPTIONS", help_heading = "Reporting")]
    track: Option<String>,
    /// Include only valid BED6 + or - strands
    #[arg(long, value_parser = ["+", "-"], help_heading = "Selection")]
    strand: Option<String>,
    /// Count strand-relative five-prime bases; requires BED6
    #[arg(long, conflicts_with_all = ["three_prime", "split"], help_heading = "Selection")]
    five_prime: bool,
    /// Count strand-relative three-prime bases; requires BED6
    #[arg(long, conflicts_with_all = ["five_prime", "split"], help_heading = "Selection")]
    three_prime: bool,
    /// Count BED12 blocks rather than the outer span
    #[arg(long, help_heading = "Selection")]
    split: bool,
}

pub(super) fn execute(args: Args, json: bool) -> Result<()> {
    super::require_named_json_output(json, args.output.as_deref())?;
    super::reject_dual_stdin(Some(&args.input), &args.genome)?;
    super::reject_output_alias(
        args.output.as_deref(),
        [Some(args.input.as_path()), Some(args.genome.as_path())],
    )?;
    let genome = read_genome(open_input(Some(&args.genome))?)?;
    let input = open_input(Some(&args.input))?;
    let options = GenomecovOptions {
        mode: if args.bedgraph {
            GenomecovMode::Bedgraph
        } else if args.bedgraph_all {
            GenomecovMode::BedgraphAll
        } else if args.depth {
            GenomecovMode::Depth
        } else if args.depth_nonzero {
            GenomecovMode::DepthNonzero
        } else {
            GenomecovMode::Histogram
        },
        strand: args.strand.map(|strand| {
            if strand == "+" {
                GenomecovStrand::Forward
            } else {
                GenomecovStrand::Reverse
            }
        }),
        end: if args.five_prime {
            Some(GenomecovEnd::FivePrime)
        } else if args.three_prime {
            Some(GenomecovEnd::ThreePrime)
        } else {
            None
        },
        split: args.split,
        max: args.max,
        scale: args.scale,
        track: args.track,
    };
    write_output(args.output.as_deref(), |output| {
        genomecov::genomecov(input, &genome, output, &options)
    })
}
