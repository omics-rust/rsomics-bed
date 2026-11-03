use std::collections::HashMap;
use std::io::{BufReader, Read};

use rsomics_common::Result;

use super::MapOptions;
use crate::bed::{BedReader, BedRecord, Strand, invalid};
use crate::{Genome, StrandFilter};

type Bounds = (i128, i128);

fn blocks(record: &BedRecord, label: &str) -> Result<Vec<Bounds>> {
    if record.field_count() < 12 {
        return Ok(vec![record.signed_bounds()]);
    }
    Ok(record
        .blocks(label)?
        .into_iter()
        .map(|(start, end)| (i128::from(start), i128::from(end)))
        .collect())
}

fn overlap(a: Bounds, b: Bounds) -> i128 {
    a.1.min(b.1) - a.0.max(b.0)
}

#[derive(Default)]
struct FileOrder {
    chromosomes: HashMap<String, usize>,
    previous: Option<String>,
    start: u64,
    width: Option<usize>,
}

struct Order<'a> {
    files: [FileOrder; 2],
    genome: Option<&'a Genome>,
    lexico_disproven: bool,
    lexico_assumed: bool,
}

impl<'a> Order<'a> {
    fn new(genome: Option<&'a Genome>) -> Self {
        Self {
            files: Default::default(),
            genome,
            lexico_disproven: false,
            lexico_assumed: false,
        }
    }

    fn observe(&mut self, file: usize, record: &BedRecord) -> Result<()> {
        let label = if file == 0 { "map A" } else { "map B" };
        let width = record.field_count();
        if self.files[file].width.is_some_and(|old| old != width) {
            return Err(invalid(format!("{label}: inconsistent BED field count")));
        }
        self.files[file].width = Some(width);
        let chrom = record.chrom();
        let previous = self.files[file].previous.as_deref();
        if previous == Some(chrom) {
            if record.start() < self.files[file].start {
                return Err(invalid(format!(
                    "{label}: starts are not sorted on {chrom}"
                )));
            }
            self.files[file].start = record.start();
            return Ok(());
        }
        if self.files[file].chromosomes.contains_key(chrom) {
            return Err(invalid(format!(
                "{label}: chromosome {chrom} reappears after another chromosome"
            )));
        }
        if let Some(genome) = self.genome {
            let rank = genome.rank(chrom).ok_or_else(|| {
                invalid(format!("{label}: chromosome {chrom} is absent from genome"))
            })?;
            if previous.is_some_and(|previous| {
                rank < genome
                    .rank(previous)
                    .expect("previous chromosome was validated")
            }) {
                return Err(invalid(format!(
                    "{label}: chromosome {chrom} is out of genome order"
                )));
            }
        } else if let Some(previous) = previous {
            let other = &self.files[1 - file].chromosomes;
            if let (Some(current), Some(old)) = (other.get(chrom), other.get(previous))
                && current < old
            {
                return Err(invalid("map input chromosome orders are inconsistent"));
            }
            if chrom < previous {
                if self.lexico_assumed {
                    return Err(invalid(
                        "map unspecified chromosome order is not lexicographic; supply --genome",
                    ));
                }
                self.lexico_disproven = true;
            }
        }
        let state = &mut self.files[file];
        state
            .chromosomes
            .insert(chrom.to_owned(), state.chromosomes.len());
        state.previous = Some(chrom.to_owned());
        state.start = record.start();
        Ok(())
    }

    fn query_after(&mut self, query: &BedRecord, database: &BedRecord) -> bool {
        if let Some(genome) = self.genome {
            return genome
                .rank(query.chrom())
                .expect("query chromosome was validated")
                > genome
                    .rank(database.chrom())
                    .expect("database chromosome was validated");
        }
        let track = &self.files[0].chromosomes;
        if let Some(&database_order) = track.get(database.chrom()) {
            return track[query.chrom()] > database_order;
        }
        if self.lexico_disproven {
            return false;
        }
        let after = query.chrom() > database.chrom();
        self.lexico_assumed |= after;
        after
    }

    fn finish(&self) -> Result<()> {
        if self.genome.is_none()
            && self.lexico_disproven
            && let Some(chrom) = self.files[1]
                .chromosomes
                .keys()
                .find(|chrom| !self.files[0].chromosomes.contains_key(*chrom))
        {
            return Err(invalid(format!(
                "map B chromosome {chrom} is absent from nonlexical A order; supply --genome"
            )));
        }
        Ok(())
    }
}

struct Entry {
    record: BedRecord,
    bounds: Bounds,
    strand: Option<Strand>,
    blocks: Vec<Bounds>,
    ordinal: usize,
}

pub(super) struct Sweep<'a, R: Read> {
    input: BedReader<BufReader<R>>,
    columns: &'a [usize],
    options: &'a MapOptions,
    order: Order<'a>,
    next: Option<Entry>,
    cache: Vec<Entry>,
    chromosome: Option<String>,
    ordinal: usize,
}

impl<'a, R: Read> Sweep<'a, R> {
    pub(super) fn new(
        input: R,
        columns: &'a [usize],
        options: &'a MapOptions,
        genome: Option<&'a Genome>,
    ) -> Result<Self> {
        let mut sweep = Self {
            input: BedReader::new(input),
            columns,
            options,
            order: Order::new(genome),
            next: None,
            cache: Vec::new(),
            chromosome: None,
            ordinal: 0,
        };
        sweep.next = sweep.read_next()?;
        if sweep.next.is_none() {
            return Err(invalid("map B contains no BED records"));
        }
        Ok(sweep)
    }

    fn read_next(&mut self) -> Result<Option<Entry>> {
        let Some(record) = self.input.next_record()? else {
            return Ok(None);
        };
        self.order.observe(1, &record)?;
        for &column in self.columns {
            record.column(column, "map B")?;
        }
        let strand = if self.options.strand == StrandFilter::Any {
            None
        } else {
            Some(record.strand("map B")?)
        };
        let blocks = if self.options.split {
            blocks(&record, "map B")?
        } else {
            Vec::new()
        };
        let entry = Entry {
            bounds: record.signed_bounds(),
            record,
            strand,
            blocks,
            ordinal: self.ordinal,
        };
        self.ordinal += 1;
        Ok(Some(entry))
    }

    pub(super) fn collect_hits(&mut self, query: &BedRecord, hits: &mut Vec<usize>) -> Result<()> {
        self.order.observe(0, query)?;
        if self.chromosome.as_deref() != Some(query.chrom()) {
            self.cache.clear();
            self.chromosome = Some(query.chrom().to_owned());
        }
        // A later equal-start point can still overlap B ending at this raw start.
        self.cache
            .retain(|entry| entry.bounds.1 >= i128::from(query.start()));
        let query_bounds = query.signed_bounds();
        let strand = if self.options.strand == StrandFilter::Any {
            None
        } else {
            Some(query.strand("map A")?)
        };
        let query_blocks = if self.options.split {
            blocks(query, "map A")?
        } else {
            Vec::new()
        };
        while let Some(entry) = &self.next {
            if entry.record.chrom() == query.chrom()
                || !self.order.query_after(query, &entry.record)
            {
                break;
            }
            self.next = self.read_next()?;
        }
        while let Some(entry) = &self.next {
            if entry.record.chrom() != query.chrom()
                || i128::from(entry.record.start()) > query_bounds.1
            {
                break;
            }
            let entry = self.next.take().expect("lookahead is present");
            if entry.bounds.1 >= i128::from(query.start()) {
                self.cache.push(entry);
            }
            self.next = self.read_next()?;
        }
        hits.clear();
        for (index, entry) in self.cache.iter().enumerate() {
            let matches_strand = match self.options.strand {
                StrandFilter::Any => true,
                StrandFilter::Same => strand == entry.strand,
                StrandFilter::Opposite => strand != entry.strand,
            };
            if !matches_strand {
                continue;
            }
            let amount = overlap(query_bounds, entry.bounds);
            if amount <= 0 {
                continue;
            }
            if self.options.split {
                if !query_blocks
                    .iter()
                    .any(|&a| entry.blocks.iter().any(|&b| overlap(a, b) > 0))
                {
                    continue;
                }
            } else if !self.accepts_fractions(amount, query_bounds, entry.bounds) {
                continue;
            }
            hits.push(index);
        }
        if self.options.sort_hits {
            hits.sort_by_key(|&index| (self.cache[index].bounds, self.cache[index].ordinal));
        }
        Ok(())
    }

    #[allow(clippy::cast_precision_loss)]
    fn accepts_fractions(&self, amount: i128, a: Bounds, b: Bounds) -> bool {
        let minimum_a = self.options.minimum_a.unwrap_or(0.0);
        let minimum_b = if self.options.reciprocal {
            minimum_a
        } else {
            self.options.minimum_b.unwrap_or(0.0)
        };
        let a = amount as f32 / (a.1 - a.0) as f32 >= minimum_a;
        let b = amount as f32 / (b.1 - b.0) as f32 >= minimum_b;
        if self.options.either { a || b } else { a && b }
    }

    pub(super) fn record(&self, hit: usize) -> &BedRecord {
        &self.cache[hit].record
    }

    pub(super) fn finish(&mut self) -> Result<()> {
        while self.next.is_some() {
            self.next = self.read_next()?;
        }
        self.order.finish()
    }
}
