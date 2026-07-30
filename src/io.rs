use std::fs::{self, File};
use std::io::{self, Read, Write};
use std::path::Path;

use rsomics_common::{Context, Result, RsomicsError};
#[cfg(unix)]
use std::os::unix::fs::PermissionsExt;
use tempfile::Builder;

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
        Some(path) => {
            let existing_permissions = match fs::metadata(path) {
                Ok(metadata) => Some(metadata.permissions()),
                Err(error) if error.kind() == io::ErrorKind::NotFound => None,
                Err(error) => {
                    return Err(RsomicsError::Io(io::Error::new(
                        error.kind(),
                        format!(
                            "reading existing output metadata {}: {error}",
                            path.display()
                        ),
                    )));
                }
            };
            let parent = path
                .parent()
                .filter(|parent| !parent.as_os_str().is_empty())
                .unwrap_or_else(|| Path::new("."));
            let mut builder = Builder::new();
            builder.prefix(".rsomics-bed-");
            #[cfg(unix)]
            if existing_permissions.is_none() {
                builder.permissions(fs::Permissions::from_mode(0o666));
            }
            if let Some(permissions) = existing_permissions.as_ref() {
                builder.permissions(permissions.clone());
            }
            let mut temporary = builder.tempfile_in(parent).rs_with_context(|| {
                format!(
                    "creating temporary output beside destination {}",
                    path.display()
                )
            })?;
            if let Some(permissions) = existing_permissions {
                temporary
                    .as_file()
                    .set_permissions(permissions)
                    .rs_with_context(|| {
                        format!("preserving permissions for output {}", path.display())
                    })?;
            }
            operation(temporary.as_file_mut())?;
            temporary
                .as_file_mut()
                .flush()
                .rs_context("flushing temporary BED output")?;
            temporary
                .as_file_mut()
                .sync_all()
                .rs_context("syncing temporary BED output")?;
            temporary.persist(path).map_err(|error| {
                let kind = error.error.kind();
                RsomicsError::Io(io::Error::new(
                    kind,
                    format!(
                        "atomically persisting output {}: {}",
                        path.display(),
                        error.error
                    ),
                ))
            })?;
            Ok(())
        }
    }
}
