//! Complement of sorted BED coverage against declared chromosome sizes.

use std::collections::HashSet;
use std::io::{BufWriter, Read, Write};

use rsomics_common::{Context, Result, RsomicsError};

use crate::bed::{BedReader, Genome, virtual_bounds};

/// Emit regions not covered by sorted BED input.
///
/// Chromosomes absent from the BED input are emitted in full and output follows
/// genome-file order.
///
/// # Errors
///
/// Returns an error for malformed BED, unknown chromosomes, coordinates beyond
/// declared chromosome sizes, input outside genome order, input/output
/// failures, or zero-length widening outside chromosome bounds.
pub fn complement(input: impl Read, genome: &Genome, output: impl Write) -> Result<()> {
    let mut reader = BedReader::new(input);
    let mut by_chrom = vec![Vec::new(); genome.chromosomes().count()];
    let mut closed = HashSet::new();
    let mut current_chrom: Option<String> = None;
    let mut current_rank: Option<usize> = None;
    let mut last_start = 0_u64;

    while let Some(record) = reader.next_record()? {
        let rank = genome.rank(&record.chrom).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "complement input chromosome {:?} is absent from the genome file",
                record.chrom
            ))
        })?;
        let size = genome.size(&record.chrom).ok_or_else(|| {
            RsomicsError::InvalidInput(format!(
                "complement input chromosome {:?} is absent from the genome file",
                record.chrom
            ))
        })?;
        if record.end > size {
            return Err(RsomicsError::InvalidInput(format!(
                "complement interval {}:{}-{} exceeds chromosome size {size}",
                record.chrom, record.start, record.end
            )));
        }

        match current_chrom.as_deref() {
            Some(chrom) if chrom == record.chrom => {
                if record.start < last_start {
                    return Err(RsomicsError::InvalidInput(format!(
                        "complement input is not sorted: {}:{} follows start {}",
                        record.chrom, record.start, last_start
                    )));
                }
            }
            Some(chrom) => {
                closed.insert(chrom.to_owned());
                if closed.contains(&record.chrom) {
                    return Err(RsomicsError::InvalidInput(format!(
                        "complement input chromosome {:?} reappears",
                        record.chrom
                    )));
                }
                let previous_rank = current_rank.ok_or_else(|| {
                    RsomicsError::InvalidInput(
                        "internal complement chromosome-order state is incomplete".to_owned(),
                    )
                })?;
                if rank < previous_rank {
                    return Err(RsomicsError::InvalidInput(format!(
                        "complement input chromosome {:?} is out of genome-file order",
                        record.chrom
                    )));
                }
                current_chrom = Some(record.chrom.clone());
                current_rank = Some(rank);
            }
            None => {
                current_chrom = Some(record.chrom.clone());
                current_rank = Some(rank);
            }
        }
        last_start = record.start;

        let (low, high) = virtual_bounds(&record, "complement")?;
        if high > size {
            return Err(RsomicsError::InvalidInput(format!(
                "widened zero-length interval {}:{}-{} exceeds chromosome size {size}",
                record.chrom, record.start, record.end
            )));
        }
        match by_chrom[rank].last_mut() {
            Some((_, previous_end)) if low <= *previous_end => {
                *previous_end = (*previous_end).max(high);
            }
            _ => by_chrom[rank].push((low, high)),
        }
    }

    let mut output = BufWriter::new(output);
    for (rank, (chrom, size)) in genome.chromosomes().enumerate() {
        let mut cursor = 0_u64;
        for &(start, end) in &by_chrom[rank] {
            if start > cursor {
                writeln!(output, "{chrom}\t{cursor}\t{start}")
                    .rs_context("writing complement BED record")?;
            }
            cursor = cursor.max(end);
        }
        if cursor < size {
            writeln!(output, "{chrom}\t{cursor}\t{size}")
                .rs_context("writing complement BED record")?;
        }
    }
    output.flush().rs_context("flushing complement BED output")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::read_genome;

    #[test]
    fn emits_missing_chromosomes() {
        let genome = read_genome(&b"chr1\t100\nchr2\t50\n"[..]).unwrap();
        let mut output = Vec::new();
        complement(&b"chr1\t10\t20\n"[..], &genome, &mut output).unwrap();
        assert_eq!(output, b"chr1\t0\t10\nchr1\t20\t100\nchr2\t0\t50\n");
    }

    #[test]
    fn interval_beyond_chromosome_fails() {
        let genome = read_genome(&b"chr1\t100\n"[..]).unwrap();
        let error = complement(&b"chr1\t90\t101\n"[..], &genome, Vec::new()).unwrap_err();
        assert!(
            error.to_string().contains("exceeds chromosome size"),
            "{error}"
        );
    }
}
