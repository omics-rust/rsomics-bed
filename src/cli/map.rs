use std::num::NonZeroUsize;
use std::path::{Path, PathBuf};

use clap::Args as ClapArgs;
use rsomics_common::{Result, RsomicsError, write_output};

use crate::io::open_input;
use crate::map::{MapOptions, Operation};
use crate::{StrandFilter, read_genome};

#[derive(Debug, ClapArgs)]
#[command(next_help_heading = "Input/output")]
pub(super) struct Args {
    /// Sorted BED queries; use - for standard input
    #[arg(short, long, value_name = "BED")]
    a: PathBuf,
    /// Sorted BED annotations; standard input is not supported
    #[arg(short, long, value_name = "BED")]
    b: PathBuf,
    /// Output BED file; omit or use - for standard output
    #[arg(short, long, value_name = "BED")]
    output: Option<PathBuf>,
    /// Chromosome order from a genome file
    #[arg(short, long, value_name = "GENOME")]
    genome: Option<PathBuf>,
    /// One-based B columns, comma separated
    #[arg(
        short = 'c',
        long,
        value_delimiter = ',',
        default_value = "5",
        help_heading = "Mapping"
    )]
    columns: Vec<usize>,
    /// Reductions, comma separated; one column or operation broadcasts
    #[arg(
        long,
        value_delimiter = ',',
        default_value = "sum",
        help_heading = "Mapping"
    )]
    operations: Vec<Operation>,
    /// Significant digits for scalar results; numeric distinct uses six
    #[arg(long, default_value = "10", help_heading = "Mapping")]
    precision: NonZeroUsize,
    /// Separator between reduced list values
    #[arg(long, default_value = ",", help_heading = "Mapping")]
    delimiter: String,
    /// Token for absent results; counts still produce zero
    #[arg(long, default_value = ".", help_heading = "Mapping")]
    null: String,
    /// Require matching valid BED6 strands
    #[arg(
        short = 's',
        long,
        conflicts_with = "opposite_strand",
        help_heading = "Selection"
    )]
    same_strand: bool,
    /// Require opposing valid BED6 strands
    #[arg(
        short = 'S',
        long,
        conflicts_with = "same_strand",
        help_heading = "Selection"
    )]
    opposite_strand: bool,
    /// Minimum fraction of A overlapped by a B candidate
    #[arg(short = 'f', long, value_name = "FRACTION", help_heading = "Selection")]
    minimum_a: Option<f32>,
    /// Minimum fraction of B overlapped by A
    #[arg(short = 'F', long, value_name = "FRACTION", help_heading = "Selection")]
    minimum_b: Option<f32>,
    /// Copy an explicit A fraction to B
    #[arg(
        short = 'r',
        long,
        conflicts_with = "minimum_b",
        help_heading = "Selection"
    )]
    reciprocal: bool,
    /// Require either fraction rather than both
    #[arg(short = 'e', long, help_heading = "Selection")]
    either: bool,
    /// Match BED12 blocks; explicit fraction selection is not supported
    #[arg(long, help_heading = "Selection")]
    split: bool,
    /// Preserve initial A header lines
    #[arg(long, help_heading = "Reporting")]
    header: bool,
    /// Sort hits by effective coordinates with stable B-order ties
    #[arg(long, help_heading = "Reporting")]
    sortout: bool,
}

pub(super) fn execute(args: Args, json: bool) -> Result<()> {
    super::require_named_json_output(json, args.output.as_deref())?;
    super::reject_stdin_b(&args.b)?;
    if args.genome.as_deref() == Some(Path::new("-")) {
        return Err(RsomicsError::ConfigError(
            "map genome must be a named file".into(),
        ));
    }
    super::reject_output_alias(
        args.output.as_deref(),
        [
            Some(args.a.as_path()),
            Some(args.b.as_path()),
            args.genome.as_deref(),
        ],
    )?;
    let genome = args
        .genome
        .as_deref()
        .map(|path| read_genome(open_input(Some(path))?))
        .transpose()?;
    let options = MapOptions {
        columns: args.columns,
        operations: args.operations,
        precision: args.precision,
        delimiter: args.delimiter.into_bytes(),
        null: args.null.into_bytes(),
        strand: if args.same_strand {
            StrandFilter::Same
        } else if args.opposite_strand {
            StrandFilter::Opposite
        } else {
            StrandFilter::Any
        },
        split: args.split,
        minimum_a: args.minimum_a,
        minimum_b: args.minimum_b,
        reciprocal: args.reciprocal,
        either: args.either,
        sort_hits: args.sortout,
        header: args.header,
    };
    let a = open_input(Some(&args.a))?;
    let b = open_input(Some(&args.b))?;
    write_output(args.output.as_deref(), |output| {
        crate::map::map(a, b, output, &options, genome.as_ref())
    })
}
