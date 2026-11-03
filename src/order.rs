use std::collections::HashSet;

use rsomics_common::Result;

use crate::Genome;
use crate::bed::{BedRecord, invalid};

pub(crate) struct InputOrder<'a> {
    operation: &'static str,
    genome: Option<&'a Genome>,
    chromosome: Option<String>,
    last_start: u64,
    last_rank: Option<usize>,
    closed: HashSet<String>,
}

impl<'a> InputOrder<'a> {
    pub(crate) fn new(operation: &'static str, genome: Option<&'a Genome>) -> Self {
        Self {
            operation,
            genome,
            chromosome: None,
            last_start: 0,
            last_rank: None,
            closed: HashSet::new(),
        }
    }

    pub(crate) fn accept(&mut self, record: &BedRecord) -> Result<bool> {
        if self.chromosome.as_deref() == Some(record.chrom()) {
            if record.start() < self.last_start {
                return Err(invalid(format!(
                    "{} input is not sorted: {}:{} follows start {}",
                    self.operation,
                    record.chrom(),
                    record.start(),
                    self.last_start
                )));
            }
            self.last_start = record.start();
            return Ok(false);
        }
        if let Some(chromosome) = &self.chromosome {
            self.closed.insert(chromosome.clone());
        }
        if self.closed.contains(record.chrom()) {
            return Err(invalid(format!(
                "{} input chromosome {:?} reappears",
                self.operation,
                record.chrom()
            )));
        }
        if let Some(genome) = self.genome {
            let rank = genome.rank(record.chrom()).ok_or_else(|| {
                invalid(format!(
                    "{} chromosome {:?} is absent from the genome file",
                    self.operation,
                    record.chrom()
                ))
            })?;
            if self.last_rank.is_some_and(|previous| rank < previous) {
                return Err(invalid(format!(
                    "{} input chromosome {:?} is out of genome order",
                    self.operation,
                    record.chrom()
                )));
            }
            self.last_rank = Some(rank);
        }
        self.chromosome = Some(record.chrom().to_owned());
        self.last_start = record.start();
        Ok(true)
    }
}
