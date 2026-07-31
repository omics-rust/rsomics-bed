use std::fs::File;
use std::io::{self, Read, Write};
use std::path::Path;

use rsomics_common::{Context, Result, write_atomic};

pub(crate) fn open_input(path: Option<&Path>) -> Result<Box<dyn Read>> {
    match path {
        None => Ok(Box::new(io::stdin())),
        Some(path) if path == Path::new("-") => Ok(Box::new(io::stdin())),
        Some(path) => File::open(path)
            .rs_with_context(|| format!("opening input {}", path.display()))
            .map(|file| Box::new(file) as Box<dyn Read>),
    }
}

pub(crate) fn with_output(
    path: Option<&Path>,
    operation: impl FnOnce(&mut dyn Write) -> Result<()>,
) -> Result<()> {
    match path {
        None => operation(&mut io::stdout().lock()),
        Some(path) if path == Path::new("-") => operation(&mut io::stdout().lock()),
        Some(path) => write_atomic(path, |output| operation(output)),
    }
}
