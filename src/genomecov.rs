//! Genome-wide coverage tracks and depth histograms from grouped BED input.

use std::collections::BTreeMap;
use std::io::{BufWriter, Read, Write};
use std::num::NonZeroUsize;

use rsomics_common::{Context, Result};

use crate::Genome;
use crate::bed::{BedReader, BedRecord, Strand, invalid, virtual_bounds};
use crate::depth::{DepthBuffer, DepthRun};
use crate::numeric::significant;

/// Genome coverage output representation.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum GenomecovMode {
    /// Per-chromosome and genome-wide depth histograms.
    #[default]
    Histogram,
    /// Nonzero depth runs in zero-based half-open coordinates.
    Bedgraph,
    /// All depth runs, including uncovered regions.
    BedgraphAll,
    /// Depth at every base, with one-based positions.
    Depth,
    /// Depth at covered bases only, with zero-based positions.
    DepthNonzero,
}

/// BED6 strand to include.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GenomecovStrand {
    /// Positive strand.
    Forward,
    /// Negative strand.
    Reverse,
}

/// Strand-relative endpoint to count instead of the entire interval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GenomecovEnd {
    /// Five-prime base.
    FivePrime,
    /// Three-prime base.
    ThreePrime,
}

/// Selection and reporting contract for [`genomecov`].
#[derive(Clone, Debug, Default)]
pub struct GenomecovOptions {
    /// Output representation.
    pub mode: GenomecovMode,
    /// Require a valid BED6 strand and include only the selected strand.
    pub strand: Option<GenomecovStrand>,
    /// Count a strand-relative endpoint; requires valid BED6 strands.
    pub end: Option<GenomecovEnd>,
    /// Count BED12 blocks rather than the outer span.
    pub split: bool,
    /// Maximum histogram bin; retains BEDTools' nonfinal bedGraph-run capping.
    pub max: Option<u32>,
    /// Finite binary32 multiplier for depth or bedGraph output.
    pub scale: Option<f32>,
    /// Browser track options, or an empty string for the default track line.
    pub track: Option<String>,
}

impl GenomecovOptions {
    fn validate(&self) -> Result<()> {
        if self.end.is_some() && self.split {
            return Err(invalid(
                "genomecov split cannot be combined with endpoint selection",
            ));
        }
        if self.scale.is_some_and(|scale| !scale.is_finite()) {
            return Err(invalid("genomecov scale must be finite"));
        }
        if self.scale.is_some() && self.mode == GenomecovMode::Histogram {
            return Err(invalid("genomecov scale requires depth or bedGraph output"));
        }
        if self.track.is_some() && !self.bedgraph() {
            return Err(invalid("genomecov track requires bedGraph output"));
        }
        if self
            .track
            .as_ref()
            .is_some_and(|text| text.contains(['\r', '\n']))
        {
            return Err(invalid("genomecov track options must fit on one line"));
        }
        if self.max.is_some_and(|max| max > i32::MAX as u32) {
            return Err(invalid(
                "genomecov max exceeds the supported signed 32-bit depth",
            ));
        }
        Ok(())
    }

    fn bedgraph(&self) -> bool {
        matches!(
            self.mode,
            GenomecovMode::Bedgraph | GenomecovMode::BedgraphAll
        )
    }

    fn max_depth(&self) -> usize {
        self.max.unwrap_or(i32::MAX as u32) as usize
    }

    fn accepts(&self, record: &BedRecord) -> Result<bool> {
        let Some(wanted) = self.strand else {
            return Ok(true);
        };
        let wanted = match wanted {
            GenomecovStrand::Forward => Strand::Forward,
            GenomecovStrand::Reverse => Strand::Reverse,
        };
        Ok(record.strand("genomecov")? == wanted)
    }

    fn add(&self, record: &BedRecord, size: u64, depth: &mut DepthBuffer) -> Result<()> {
        if record.end() > size {
            return Err(invalid(format!(
                "genomecov interval on {} exceeds chromosome size {size}",
                record.chrom()
            )));
        }
        let bounds = virtual_bounds(record, "genomecov")?;
        if self.split && record.field_count() == 12 {
            for (start, end) in record.blocks_in("genomecov", bounds)? {
                depth.add(start, end.min(size));
            }
        } else {
            let (mut start, mut end) = bounds;
            if let Some(endpoint) = self.end {
                let strand = record.strand("genomecov")?;
                let left = matches!(
                    (endpoint, strand),
                    (GenomecovEnd::FivePrime, Strand::Forward)
                        | (GenomecovEnd::ThreePrime, Strand::Reverse)
                );
                if left {
                    end = start + 1;
                } else {
                    start = end - 1;
                }
            }
            depth.add(start, end.min(size));
        }
        Ok(())
    }
}

/// Write coverage for chromosome-grouped BED input and a declared genome.
///
/// Chromosome groups retain encounter order after strand selection;
/// unobserved chromosomes follow in genome-file order. Genome lookup and
/// geometry checks apply to selected records. Starts within a group need not
/// be sorted. Memory scales with one chromosome's interval events, not length.
///
/// # Errors
///
/// Returns parse/I/O errors, invalid option combinations, chromosome reentry,
/// unknown chromosomes, out-of-bounds intervals, invalid required strands or
/// points whose virtual bounds underflow or overflow.
pub fn genomecov(
    input: impl Read,
    genome: &Genome,
    output: impl Write,
    options: &GenomecovOptions,
) -> Result<()> {
    options.validate()?;
    let mut output = BufWriter::new(output);
    if let Some(track) = &options.track {
        write!(output, "track type=bedGraph")?;
        if !track.is_empty() {
            write!(output, " {track}")?;
        }
        writeln!(output)?;
    }
    let chromosomes: Vec<_> = genome.chromosomes().collect();
    let mut visited = vec![false; chromosomes.len()];
    let total = chromosomes.iter().try_fold(0u64, |sum, (_, size)| {
        sum.checked_add(*size)
            .ok_or_else(|| invalid("genomecov genome size overflows u64"))
    })?;
    let mut histogram = BTreeMap::new();
    let mut depth = DepthBuffer::default();
    let mut current = None;
    let mut reader = BedReader::new(input);
    while let Some(record) = reader.next_record()? {
        if !options.accepts(&record)? {
            continue;
        }
        let rank = genome.rank(record.chrom()).ok_or_else(|| {
            invalid(format!(
                "genomecov chromosome {:?} is absent from the genome",
                record.chrom()
            ))
        })?;
        if current != Some(rank) {
            if visited[rank] {
                return Err(invalid(format!(
                    "genomecov chromosome {:?} reappears after another group",
                    record.chrom()
                )));
            }
            if let Some(previous) = current {
                report(
                    chromosomes[previous],
                    &mut depth,
                    options,
                    &mut histogram,
                    &mut output,
                )?;
            }
            depth.clear_events();
            visited[rank] = true;
            current = Some(rank);
        }
        options.add(&record, chromosomes[rank].1, &mut depth)?;
    }
    if let Some(rank) = current {
        report(
            chromosomes[rank],
            &mut depth,
            options,
            &mut histogram,
            &mut output,
        )?;
    }
    depth.clear_events();
    for (rank, chromosome) in chromosomes.iter().enumerate() {
        if !visited[rank] {
            report(
                *chromosome,
                &mut depth,
                options,
                &mut histogram,
                &mut output,
            )?;
        }
    }
    if options.mode == GenomecovMode::Histogram {
        write_histogram("genome", total, &histogram, &mut output)?;
    }
    output.flush().rs_context("flushing genome coverage")
}

fn report(
    (chromosome, size): (&str, u64),
    depth: &mut DepthBuffer,
    options: &GenomecovOptions,
    genome_histogram: &mut BTreeMap<usize, u64>,
    output: &mut impl Write,
) -> Result<()> {
    depth.reduce(0, size);
    if options.mode == GenomecovMode::Histogram {
        let mut histogram = BTreeMap::new();
        for (&depth, &bases) in depth.histogram() {
            *histogram.entry(depth.min(options.max_depth())).or_default() += bases;
            *genome_histogram
                .entry(depth.min(options.max_depth()))
                .or_default() += bases;
        }
        return write_histogram(chromosome, size, &histogram, output);
    }
    let mut pending: Option<DepthRun> = None;
    for &run in depth.runs() {
        match &mut pending {
            Some(previous) if previous.depth == run.depth => previous.end = run.end,
            _ => {
                if let Some(previous) = pending {
                    write_run(chromosome, size, previous, options, output)?;
                }
                pending = Some(run);
            }
        }
    }
    if let Some(run) = pending {
        write_run(chromosome, size, run, options, output)?;
    }
    Ok(())
}

#[allow(clippy::cast_precision_loss)]
fn write_histogram(
    chromosome: &str,
    size: u64,
    histogram: &BTreeMap<usize, u64>,
    output: &mut impl Write,
) -> Result<()> {
    for (depth, bases) in histogram {
        let fraction = f64::from(*bases as f32 / size as f32);
        writeln!(
            output,
            "{chromosome}\t{depth}\t{bases}\t{size}\t{}",
            significant(fraction, NonZeroUsize::new(6).unwrap())
        )?;
    }
    Ok(())
}

#[allow(clippy::cast_precision_loss)]
fn write_run(
    chromosome: &str,
    size: u64,
    run: DepthRun,
    options: &GenomecovOptions,
    output: &mut impl Write,
) -> Result<()> {
    if run.depth == 0
        && matches!(
            options.mode,
            GenomecovMode::Bedgraph | GenomecovMode::DepthNonzero
        )
    {
        return Ok(());
    }
    let depth = if options.bedgraph() && run.end != size {
        run.depth.min(options.max_depth())
    } else {
        run.depth
    };
    let scale = options.scale.unwrap_or(1.0);
    let value = depth as f32 * scale;
    let value = if scale == 1.0 {
        format!("{value:.0}")
    } else {
        significant(f64::from(value), NonZeroUsize::new(6).unwrap())
    };
    if options.bedgraph() {
        writeln!(output, "{chromosome}\t{}\t{}\t{value}", run.start, run.end)?;
    } else {
        let offset = u64::from(options.mode == GenomecovMode::Depth);
        for position in run.start..run.end {
            writeln!(output, "{chromosome}\t{}\t{value}", position + offset)?;
        }
    }
    Ok(())
}
