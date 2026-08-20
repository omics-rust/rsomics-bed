use std::collections::HashMap;
use std::fmt::Display;
use std::io::{BufRead, BufReader, Read, Write};

use rsomics_common::{Context, Result, RsomicsError};
use rsomics_intervals::Interval;

#[derive(Debug, Clone, PartialEq, Eq)]
pub(crate) struct BedRecord {
    interval: Interval,
    raw: Vec<u8>,
    suffix_start: Option<usize>,
    line_number: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum Strand {
    Forward,
    Reverse,
}

impl BedRecord {
    pub(crate) fn chrom(&self) -> &str {
        self.interval.chrom()
    }

    pub(crate) fn start(&self) -> u64 {
        self.interval.start()
    }

    pub(crate) fn end(&self) -> u64 {
        self.interval.end()
    }

    pub(crate) fn line_number(&self) -> usize {
        self.line_number
    }

    pub(crate) fn field_count(&self) -> usize {
        self.raw.split(|&byte| byte == b'\t').count()
    }

    pub(crate) fn name(&self, label: &str) -> Result<&[u8]> {
        self.field(3, label, "name")
    }

    pub(crate) fn strand(&self, label: &str) -> Result<Strand> {
        match self.field(5, label, "strand")? {
            b"+" => Ok(Strand::Forward),
            b"-" => Ok(Strand::Reverse),
            value => Err(invalid(format!(
                "{label} BED line {}: invalid strand {:?}",
                self.line_number,
                String::from_utf8_lossy(value)
            ))),
        }
    }

    fn field(&self, index: usize, label: &str, field: &str) -> Result<&[u8]> {
        let value = self
            .raw
            .split(|&byte| byte == b'\t')
            .nth(index)
            .ok_or_else(|| {
                invalid(format!(
                    "{label} BED line {}: missing {field}",
                    self.line_number
                ))
            })?;
        if value.is_empty() {
            return Err(invalid(format!(
                "{label} BED line {}: empty {field}",
                self.line_number
            )));
        }
        Ok(value)
    }

    pub(crate) fn write_raw(&self, output: &mut dyn Write) -> Result<()> {
        output
            .write_all(&self.raw)
            .rs_context("writing BED record")?;
        output.write_all(b"\n").rs_context("writing BED record")
    }

    pub(crate) fn write_with_coords(
        &self,
        output: &mut dyn Write,
        start: u64,
        end: u64,
    ) -> Result<()> {
        write!(output, "{}\t{start}\t{end}", self.chrom()).rs_context("writing BED record")?;
        if let Some(index) = self.suffix_start {
            output
                .write_all(&self.raw[index..])
                .rs_context("writing BED record")?;
        }
        output.write_all(b"\n").rs_context("writing BED record")
    }

    pub(crate) fn write_joined(&self, output: &mut dyn Write, other: &Self) -> Result<()> {
        output
            .write_all(&self.raw)
            .rs_context("writing BED record")?;
        output.write_all(b"\t").rs_context("writing BED record")?;
        output
            .write_all(&other.raw)
            .rs_context("writing BED record")?;
        output.write_all(b"\n").rs_context("writing BED record")
    }

    pub(crate) fn write_column(&self, output: &mut dyn Write, value: impl Display) -> Result<()> {
        output
            .write_all(&self.raw)
            .rs_context("writing BED record")?;
        writeln!(output, "\t{value}").rs_context("writing BED record")
    }
}

pub(crate) struct BedReader<R: BufRead> {
    input: R,
    line: Vec<u8>,
    line_number: usize,
}

pub(crate) struct BedCoordinates<'a> {
    interval: Interval<&'a str>,
}

impl BedCoordinates<'_> {
    pub(crate) fn chrom(&self) -> &str {
        self.interval.chrom()
    }

    pub(crate) fn start(&self) -> u64 {
        self.interval.start()
    }

    pub(crate) fn end(&self) -> u64 {
        self.interval.end()
    }
}

struct ParsedRecord<'a> {
    interval: Interval<&'a str>,
    raw: &'a [u8],
    suffix_start: Option<usize>,
    line_number: usize,
}

impl<R: Read> BedReader<BufReader<R>> {
    pub(crate) fn new(input: R) -> Self {
        Self {
            input: BufReader::new(input),
            line: Vec::with_capacity(256),
            line_number: 0,
        }
    }
}

impl<R: BufRead> BedReader<R> {
    pub(crate) fn next_record(&mut self) -> Result<Option<BedRecord>> {
        Ok(self.next_parsed()?.map(|record| {
            let interval = Interval::new(
                record.interval.chrom().to_string(),
                record.interval.start(),
                record.interval.end(),
            )
            .expect("parsed BED interval remains valid when its chromosome is owned");
            BedRecord {
                interval,
                raw: record.raw.to_vec(),
                suffix_start: record.suffix_start,
                line_number: record.line_number,
            }
        }))
    }

    pub(crate) fn next_coordinates(&mut self) -> Result<Option<BedCoordinates<'_>>> {
        Ok(self.next_parsed()?.map(|record| BedCoordinates {
            interval: record.interval,
        }))
    }

    fn next_parsed(&mut self) -> Result<Option<ParsedRecord<'_>>> {
        loop {
            self.line.clear();
            if self
                .input
                .read_until(b'\n', &mut self.line)
                .rs_context("reading BED input")?
                == 0
            {
                return Ok(None);
            }
            self.line_number += 1;
            while matches!(self.line.last(), Some(b'\n' | b'\r')) {
                self.line.pop();
            }
            if is_skippable(&self.line) {
                continue;
            }
            return parse_record(&self.line, self.line_number).map(Some);
        }
    }
}

pub(crate) fn read_records(input: impl Read) -> Result<Vec<BedRecord>> {
    let mut reader = BedReader::new(input);
    let mut records = Vec::new();
    while let Some(record) = reader.next_record()? {
        records.push(record);
    }
    Ok(records)
}

fn is_skippable(line: &[u8]) -> bool {
    line.is_empty() || line[0] == b'#' || line.starts_with(b"track") || line.starts_with(b"browser")
}

fn parse_record(line: &[u8], line_number: usize) -> Result<ParsedRecord<'_>> {
    let first_tab = line.iter().position(|&byte| byte == b'\t').ok_or_else(|| {
        invalid(format!(
            "BED line {line_number}: expected at least three tab-separated columns"
        ))
    })?;
    let second_rel = line[first_tab + 1..]
        .iter()
        .position(|&byte| byte == b'\t')
        .ok_or_else(|| {
            invalid(format!(
                "BED line {line_number}: expected at least three tab-separated columns"
            ))
        })?;
    let second_tab = first_tab + 1 + second_rel;
    let third_tab = line[second_tab + 1..]
        .iter()
        .position(|&byte| byte == b'\t')
        .map(|offset| second_tab + 1 + offset);

    let chrom_bytes = &line[..first_tab];
    if chrom_bytes.is_empty() {
        return Err(invalid(format!("BED line {line_number}: empty chromosome")));
    }
    let chrom = std::str::from_utf8(chrom_bytes).map_err(|error| {
        invalid(format!(
            "BED line {line_number}: chromosome is not UTF-8: {error}"
        ))
    })?;
    let context = format!("BED line {line_number}");
    let start = parse_u64(&line[first_tab + 1..second_tab], &context, "start")?;
    let end_bytes = third_tab.map_or(&line[second_tab + 1..], |index| {
        &line[second_tab + 1..index]
    });
    let end = parse_u64(end_bytes, &context, "end")?;
    let interval = Interval::new(chrom, start, end).map_err(|_| {
        invalid(format!(
            "BED line {line_number}: start {start} is greater than end {end}"
        ))
    })?;

    Ok(ParsedRecord {
        interval,
        raw: line,
        suffix_start: third_tab,
        line_number,
    })
}

fn parse_u64(bytes: &[u8], context: &str, field: &str) -> Result<u64> {
    if bytes.is_empty() {
        return Err(invalid(format!("{context}: empty {field}")));
    }
    let mut value = 0_u64;
    for &byte in bytes {
        let digit = byte.wrapping_sub(b'0');
        if digit > 9 {
            return Err(invalid(format!(
                "{context}: invalid {field} {:?}",
                String::from_utf8_lossy(bytes)
            )));
        }
        value = value
            .checked_mul(10)
            .and_then(|current| current.checked_add(u64::from(digit)))
            .ok_or_else(|| invalid(format!("{context}: {field} overflows u64")))?;
    }
    Ok(value)
}

/// Chromosome names, sizes, and declared ordering from a genome file.
#[derive(Debug, Clone)]
pub struct Genome {
    entries: Vec<(String, u64)>,
    ranks: HashMap<String, usize>,
}

impl Genome {
    /// Iterate over chromosomes in genome-file order.
    pub fn chromosomes(&self) -> impl Iterator<Item = (&str, u64)> {
        self.entries
            .iter()
            .map(|(chrom, size)| (chrom.as_str(), *size))
    }

    /// Return the declared size of `chrom`.
    #[must_use]
    pub fn size(&self, chrom: &str) -> Option<u64> {
        self.rank(chrom).map(|rank| self.entries[rank].1)
    }

    pub(crate) fn rank(&self, chrom: &str) -> Option<usize> {
        self.ranks.get(chrom).copied()
    }
}

/// Read a two-column, tab-delimited chromosome-size file.
///
/// Blank lines and lines beginning with `#` are ignored. Chromosome order is
/// retained for complement output.
///
/// # Errors
///
/// Returns an error for malformed rows, invalid UTF-8 chromosome names,
/// duplicate chromosomes, overflowing sizes, or an empty genome file.
pub fn read_genome(input: impl Read) -> Result<Genome> {
    let mut reader = BufReader::new(input);
    let mut line = Vec::with_capacity(128);
    let mut line_number = 0_usize;
    let mut entries = Vec::new();
    let mut ranks = HashMap::new();

    loop {
        line.clear();
        if reader
            .read_until(b'\n', &mut line)
            .rs_context("reading chromosome sizes")?
            == 0
        {
            break;
        }
        line_number += 1;
        while matches!(line.last(), Some(b'\n' | b'\r')) {
            line.pop();
        }
        if line.is_empty() || line[0] == b'#' {
            continue;
        }
        let mut fields = line.split(|&byte| byte == b'\t');
        let chrom_bytes = fields.next().expect("split yields one field");
        let size_bytes = fields.next().ok_or_else(|| {
            invalid(format!(
                "genome line {line_number}: expected chromosome and size"
            ))
        })?;
        if fields.next().is_some() {
            return Err(invalid(format!(
                "genome line {line_number}: expected exactly two columns"
            )));
        }
        let chrom = std::str::from_utf8(chrom_bytes).map_err(|error| {
            invalid(format!(
                "genome line {line_number}: chromosome is not UTF-8: {error}"
            ))
        })?;
        if chrom.is_empty() {
            return Err(invalid(format!(
                "genome line {line_number}: empty chromosome"
            )));
        }
        let size = parse_u64(
            size_bytes,
            &format!("genome line {line_number}"),
            "chromosome size",
        )?;
        if ranks.contains_key(chrom) {
            return Err(invalid(format!(
                "genome line {line_number}: duplicate chromosome {chrom:?}"
            )));
        }
        let rank = entries.len();
        entries.push((chrom.to_owned(), size));
        ranks.insert(chrom.to_owned(), rank);
    }

    if entries.is_empty() {
        return Err(invalid("genome file contains no chromosome sizes"));
    }
    Ok(Genome { entries, ranks })
}

pub(crate) fn invalid(message: impl Into<String>) -> RsomicsError {
    RsomicsError::InvalidInput(message.into())
}

pub(crate) fn virtual_bounds(record: &BedRecord, operation: &str) -> Result<(u64, u64)> {
    if record.start() != record.end() {
        return Ok((record.start(), record.end()));
    }
    let low = record.start().checked_sub(1).ok_or_else(|| {
        invalid(format!(
            "{operation} zero-length interval {}:{}-{} widens below coordinate zero",
            record.chrom(),
            record.start(),
            record.end()
        ))
    })?;
    let high = record.end().checked_add(1).ok_or_else(|| {
        invalid(format!(
            "{operation} zero-length interval {}:{}-{} widens beyond u64",
            record.chrom(),
            record.start(),
            record.end()
        ))
    })?;
    Ok((low, high))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parser_preserves_trailing_columns() {
        let records = read_records(&b"chr1\t10\t20\tname\t0\t+\r\n"[..]).unwrap();
        assert_eq!(records.len(), 1);
        let mut output = Vec::new();
        records[0].write_with_coords(&mut output, 12, 18).unwrap();
        assert_eq!(output, b"chr1\t12\t18\tname\t0\t+\n");
    }

    #[test]
    fn parser_reports_physical_line_number() {
        let error = read_records(&b"# header\nchr1\t10\t20\nchr1\tbad\t30\n"[..]).unwrap_err();
        assert!(error.to_string().contains("line 3"), "{error}");
    }

    #[test]
    fn optional_fields_are_checked_only_when_requested() {
        let mut records =
            read_records(&b"# header\nchr1\t10\t20\nchr1\t20\t30\tname\t0\t-\n"[..]).unwrap();
        let named = records.pop().unwrap();
        let bed3 = records.pop().unwrap();

        let error = bed3.strand("cluster").unwrap_err();
        assert!(
            error
                .to_string()
                .contains("cluster BED line 2: missing strand"),
            "{error}"
        );
        assert_eq!(named.strand("closest A").unwrap(), Strand::Reverse);
    }

    #[test]
    fn names_are_checked_only_when_requested() {
        let records = read_records(&b"chr1\t10\t20\nchr1\t30\t40\t\n"[..]).unwrap();
        let missing = records[0].name("closest A").unwrap_err();
        assert!(
            missing.to_string().contains("closest A BED line 1"),
            "{missing}"
        );
        let empty = records[1].name("closest B").unwrap_err();
        assert!(
            empty.to_string().contains("closest B BED line 2"),
            "{empty}"
        );

        let named = read_records(&b"chr1\t50\t60\twanted\n"[..]).unwrap();
        assert_eq!(named[0].name("closest B").unwrap(), b"wanted");
    }

    #[test]
    fn invalid_strands_fail_with_physical_lines() {
        let records = read_records(&b"chr1\t20\t30\tname\t0\t.\n"[..]).unwrap();
        let strand_error = records[0].strand("window B").unwrap_err();
        assert!(
            strand_error
                .to_string()
                .contains("window B BED line 1: invalid strand"),
            "{strand_error}"
        );
    }

    #[test]
    fn appended_columns_preserve_original_fields() {
        let records = read_records(&b"chr1\t10\t20\ta\t0\t+\r\n"[..]).unwrap();
        let mut appended = Vec::new();
        records[0].write_column(&mut appended, 42).unwrap();
        assert_eq!(appended, b"chr1\t10\t20\ta\t0\t+\t42\n");
    }

    #[test]
    fn browser_and_track_prefixes_are_skipped_like_bedtools() {
        let records = read_records(
            &b"track1\t1\t2\n\
              track\t3\t4\n\
              browser_alt\t5\t6\n\
              browser\t7\t8\n\
              Track\t9\t10\n"[..],
        )
        .unwrap();
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].chrom(), "Track");
    }

    #[test]
    fn genome_rejects_duplicates() {
        let error = read_genome(&b"chr1\t10\nchr1\t20\n"[..]).unwrap_err();
        assert!(
            error.to_string().contains("duplicate chromosome"),
            "{error}"
        );
    }
}
