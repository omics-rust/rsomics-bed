use std::collections::BTreeMap;
use std::io::Write;
use std::num::NonZeroUsize;

use rsomics_common::{Context, Result, RsomicsError};

use super::Operation;

pub(super) struct Format<'a> {
    pub precision: NonZeroUsize,
    pub delimiter: &'a [u8],
    pub null: &'a [u8],
}

pub(super) struct Reducer<'a> {
    values: &'a [&'a [u8]],
    format: Format<'a>,
    numbers: Option<Vec<f64>>,
    sorted: Option<Vec<f64>>,
    frequencies: Option<BTreeMap<&'a [u8], usize>>,
}

impl<'a> Reducer<'a> {
    pub fn new(values: &'a [&'a [u8]], format: Format<'a>) -> Self {
        Self {
            values,
            format,
            numbers: None,
            sorted: None,
            frequencies: None,
        }
    }

    fn numbers(&mut self) -> Result<&[f64]> {
        if self.numbers.is_none() {
            let numbers = self
                .values
                .iter()
                .map(|value| {
                    let number = std::str::from_utf8(value)?.parse::<f64>()?;
                    if !number.is_finite() {
                        return Err(RsomicsError::InvalidInput(
                            "map numeric value must be finite".into(),
                        ));
                    }
                    Ok(number)
                })
                .collect::<Result<Vec<_>>>()?;
            self.numbers = Some(numbers);
        }
        Ok(self.numbers.as_deref().unwrap())
    }

    fn sorted(&mut self) -> Result<&[f64]> {
        if self.sorted.is_none() {
            let mut numbers = self.numbers()?.to_vec();
            numbers.sort_by(|a, b| a.partial_cmp(b).unwrap());
            self.sorted = Some(numbers);
        }
        Ok(self.sorted.as_deref().unwrap())
    }

    fn frequencies(&mut self) -> &BTreeMap<&'a [u8], usize> {
        self.frequencies.get_or_insert_with(|| {
            let mut counts = BTreeMap::new();
            for &value in self.values {
                *counts.entry(value).or_insert(0) += 1;
            }
            counts
        })
    }

    pub fn write(&mut self, op: Operation, out: &mut impl Write) -> Result<()> {
        use Operation::*;
        let delimiter = self.format.delimiter;
        if self.values.is_empty() {
            return out
                .write_all(if matches!(op, Count | CountDistinct) {
                    b"0"
                } else {
                    self.format.null
                })
                .rs_context("writing map reduction");
        }
        let scalar = match op {
            Count => Some(self.values.len().to_string().into_bytes()),
            CountDistinct => Some(self.frequencies().len().to_string().into_bytes()),
            Sum | Mean | Stdev | SampleStdev => {
                let values = self.numbers()?;
                let sum: f64 = values.iter().sum();
                #[allow(clippy::cast_precision_loss)]
                let divisor = f64::from(values.len() as f32);
                let mean = sum / divisor;
                let value = match op {
                    Sum => sum,
                    Mean => mean,
                    _ => {
                        let squared: f64 = values.iter().map(|n| (n - mean) * (n - mean)).sum();
                        (squared
                            / if op == SampleStdev {
                                divisor - 1.0
                            } else {
                                divisor
                            })
                        .sqrt()
                    }
                };
                Some(self.scalar(value))
            }
            Median => {
                let values = self.sorted()?;
                let value = if values.len() % 2 == 0 {
                    (values[values.len() / 2 - 1] + values[values.len() / 2]) / 2.0
                } else {
                    values[values.len() / 2]
                };
                Some(self.scalar(value))
            }
            Min | Max => {
                let values = self.numbers()?;
                let value = values
                    .iter()
                    .copied()
                    .reduce(|a, b| {
                        if (op == Min && b < a) || (op == Max && b > a) {
                            b
                        } else {
                            a
                        }
                    })
                    .unwrap();
                Some(self.scalar(value))
            }
            AbsMin | AbsMax => {
                let values = self.numbers()?;
                let value = values
                    .iter()
                    .map(|n| n.abs())
                    .reduce(if op == AbsMin { f64::min } else { f64::max })
                    .unwrap();
                Some(self.scalar(value))
            }
            _ => None,
        };
        if let Some(value) = scalar {
            return out.write_all(&value).rs_context("writing map reduction");
        }
        match op {
            First | Last => out.write_all(if op == First {
                self.values[0]
            } else {
                self.values[self.values.len() - 1]
            })?,
            Collapse | Concat => {
                for (index, value) in self.values.iter().enumerate() {
                    if index != 0 && op == Collapse {
                        out.write_all(delimiter)?;
                    }
                    out.write_all(value)?;
                }
            }
            Mode | Antimode => {
                let mut selected = None;
                for (&value, &count) in self.frequencies() {
                    if selected
                        .is_none_or(|(_, old)| if op == Mode { count > old } else { count < old })
                    {
                        selected = Some((value, count));
                    }
                }
                out.write_all(selected.unwrap().0)?;
            }
            Distinct | DistinctOnly => {
                for (index, (&value, &count)) in self.frequencies().iter().enumerate() {
                    if op == DistinctOnly && count != 1 {
                        continue;
                    }
                    if index != 0 {
                        out.write_all(delimiter)?;
                    }
                    out.write_all(value)?;
                }
            }
            FreqAsc | FreqDesc => {
                let mut entries: Vec<_> = self.frequencies().iter().collect();
                entries.sort_by(|(a, ac), (b, bc)| {
                    (if op == FreqAsc {
                        ac.cmp(bc)
                    } else {
                        bc.cmp(ac)
                    })
                    .then_with(|| a.cmp(b))
                });
                for (index, (value, count)) in entries.iter().enumerate() {
                    if index != 0 {
                        out.write_all(delimiter)?;
                    }
                    out.write_all(value)?;
                    write!(out, ":{count}")?;
                }
            }
            DistinctSortNum | DistinctSortNumDesc => {
                let mut values = self.sorted()?.to_vec();
                values.dedup();
                if op == DistinctSortNumDesc {
                    values.reverse();
                }
                for (index, value) in values.iter().enumerate() {
                    if index != 0 {
                        out.write_all(delimiter)?;
                    }
                    out.write_all(significant(*value, NonZeroUsize::new(6).unwrap()).as_bytes())?;
                }
            }
            _ => unreachable!("scalar operation already emitted"),
        }
        Ok(())
    }

    fn scalar(&self, value: f64) -> Vec<u8> {
        if value.is_nan() {
            self.format.null.to_vec()
        } else {
            significant(value, self.format.precision).into_bytes()
        }
    }
}

fn significant(value: f64, precision: NonZeroUsize) -> String {
    if !value.is_finite() || value == 0.0 {
        return value.to_string();
    }
    let scientific = format!("{:.*e}", precision.get() - 1, value);
    let (mantissa, exponent) = scientific.split_once('e').unwrap();
    let exponent: i32 = exponent.parse().unwrap();
    if exponent < -4 || usize::try_from(exponent).is_ok_and(|e| e >= precision.get()) {
        let mantissa = if mantissa.contains('.') {
            mantissa.trim_end_matches('0').trim_end_matches('.')
        } else {
            mantissa
        };
        return format!("{mantissa}e{exponent:+03}");
    }
    let negative = mantissa.starts_with('-');
    let digits = mantissa.trim_start_matches('-').replace('.', "");
    let position = exponent + 1;
    let mut output = if position <= 0 {
        format!(
            "0.{}{}",
            "0".repeat(usize::try_from(-position).unwrap()),
            digits
        )
    } else {
        let position = usize::try_from(position).unwrap();
        if position >= digits.len() {
            format!("{digits}{}", "0".repeat(position - digits.len()))
        } else {
            format!("{}.{}", &digits[..position], &digits[position..])
        }
    };
    if output.contains('.') {
        output.truncate(output.trim_end_matches('0').trim_end_matches('.').len());
    }
    if negative {
        output.insert(0, '-');
    }
    output
}

#[cfg(test)]
mod tests {
    use std::num::NonZeroUsize;

    use super::*;
    use crate::map::Operation;

    const ALL: [Operation; 23] = [
        Operation::Sum,
        Operation::Mean,
        Operation::Stdev,
        Operation::SampleStdev,
        Operation::Median,
        Operation::Mode,
        Operation::Antimode,
        Operation::Min,
        Operation::Max,
        Operation::AbsMin,
        Operation::AbsMax,
        Operation::Count,
        Operation::Distinct,
        Operation::CountDistinct,
        Operation::DistinctOnly,
        Operation::DistinctSortNum,
        Operation::DistinctSortNumDesc,
        Operation::Collapse,
        Operation::Concat,
        Operation::FreqAsc,
        Operation::FreqDesc,
        Operation::First,
        Operation::Last,
    ];

    fn render(
        values: &[&[u8]],
        operations: &[Operation],
        precision: usize,
    ) -> rsomics_common::Result<Vec<u8>> {
        let format = Format {
            precision: NonZeroUsize::new(precision).unwrap(),
            delimiter: b",",
            null: b".",
        };
        let mut reducer = Reducer::new(values, format);
        let mut output = Vec::new();
        for (index, &operation) in operations.iter().enumerate() {
            if index != 0 {
                output.push(b'\t');
            }
            reducer.write(operation, &mut output)?;
        }
        Ok(output)
    }

    #[test]
    fn all23_finite_reducers_match_complete_source_row() {
        let values: &[&[u8]] = &[
            b"2",
            b"-3",
            b"2",
            b"1.234567890123",
            b"1.234567890123",
            b"4",
            b"4",
            b"0",
        ];
        assert_eq!(render(values, &ALL, 10).unwrap(), b"11.46913578\t1.433641973\t2.109670629\t2.255332769\t1.617283945\t1.234567890123\t-3\t-3\t4\t0\t4\t8\t-3,0,1.234567890123,2,4\t5\t-3,0\t-3,0,1.23457,2,4\t4,2,1.23457,0,-3\t2,-3,2,1.234567890123,1.234567890123,4,4,0\t2-321.2345678901231.234567890123440\t-3:1,0:1,1.234567890123:2,2:2,4:2\t1.234567890123:2,2:2,4:2,-3:1,0:1\t2\t0");
    }

    #[test]
    fn lexical_reducers_preserve_source_ties_and_distinct_only_separator() {
        use Operation::*;
        let operations = [
            Count,
            Distinct,
            CountDistinct,
            DistinctOnly,
            Collapse,
            Concat,
            Mode,
            Antimode,
            FreqAsc,
            FreqDesc,
            First,
            Last,
        ];
        assert_eq!(
            render(&[b"a", b"a", b"c", b"b", b"c"], &operations, 10).unwrap(),
            b"5\ta,b,c\t3\t,b\ta,a,c,b,c\taacbc\ta\tb\tb:1,a:2,c:2\ta:2,c:2,b:1\ta\tc"
        );
        assert_eq!(
            render(&[b"\xff", b""], &[Collapse, CountDistinct], 10).unwrap(),
            b"\xff,\t2"
        );
    }

    #[test]
    fn empty_and_singleton_columns_use_operation_specific_results() {
        let empty = ALL
            .iter()
            .map(|op| {
                if matches!(op, Operation::Count | Operation::CountDistinct) {
                    "0"
                } else {
                    "."
                }
            })
            .collect::<Vec<_>>()
            .join("\t");
        assert_eq!(render(&[], &ALL, 10).unwrap(), empty.as_bytes());
        assert_eq!(
            render(&[b"7"], &ALL, 10).unwrap(),
            b"7\t7\t0\t.\t7\t7\t7\t7\t7\t7\t7\t1\t7\t1\t7\t7\t7\t7\t7\t7:1\t7:1\t7\t7"
        );
    }

    #[test]
    fn numeric_selection_rejects_invalid_and_nonfinite_values() {
        for token in [b"invalid".as_slice(), b"NaN", b"inf", b"1e999", b" 1"] {
            assert!(matches!(
                render(&[b"1", token], &[Operation::Sum], 10),
                Err(RsomicsError::InvalidInput(_))
            ));
            assert!(render(&[token], &[Operation::First], 10).is_ok());
        }
    }

    #[test]
    fn numeric_formats_keep_scalar_and_sorted_precision_distinct() {
        assert_eq!(
            render(
                &[b"1.234567890123"],
                &[Operation::Sum, Operation::DistinctSortNum],
                10
            )
            .unwrap(),
            b"1.23456789\t1.23457"
        );
        assert_eq!(
            render(&[b"0.0000123456789"], &[Operation::Sum], 6).unwrap(),
            b"1.23457e-05"
        );
        assert_eq!(
            render(&[b"999999.9"], &[Operation::Sum], 6).unwrap(),
            b"1e+06"
        );
        assert_eq!(
            render(&[b"1e308", b"1e308"], &[Operation::Sum], 10).unwrap(),
            b"inf"
        );
    }

    #[test]
    fn configured_format_preserves_raw_delimiter_and_null_bytes() {
        let format = Format {
            precision: NonZeroUsize::new(3).unwrap(),
            delimiter: b"|",
            null: b"\xff",
        };
        let mut reducer = Reducer::new(&[b"1", b"2"], format);
        let mut output = Vec::new();
        reducer.write(Operation::Distinct, &mut output).unwrap();
        assert_eq!(output, b"1|2");
        let format = Format {
            precision: NonZeroUsize::new(3).unwrap(),
            delimiter: b"|",
            null: b"\xff",
        };
        let mut singleton = Reducer::new(&[b"7"], format);
        output.clear();
        singleton
            .write(Operation::SampleStdev, &mut output)
            .unwrap();
        assert_eq!(output, b"\xff");
    }
}
