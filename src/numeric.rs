use std::num::NonZeroUsize;

pub(crate) fn significant(value: f64, precision: NonZeroUsize) -> String {
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
