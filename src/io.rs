use std::fs::File;
use std::io::{self, Read};
use std::path::Path;

use rsomics_common::{Context, Result};

pub(crate) fn open_input(path: Option<&Path>) -> Result<Box<dyn Read>> {
    match path {
        None => Ok(Box::new(io::stdin())),
        Some(path) if path == Path::new("-") => Ok(Box::new(io::stdin())),
        Some(path) => File::open(path)
            .rs_with_context(|| format!("opening input {}", path.display()))
            .map(|file| Box::new(file) as Box<dyn Read>),
    }
}
