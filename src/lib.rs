//! BED parsing and interval-algebra operations for the `rsomics-bed` product.

#![forbid(unsafe_code)]
#![warn(missing_docs)]

mod bed;
mod cli;
pub mod closest;
pub mod cluster;
pub mod complement;
pub mod intersect;
mod interval_index;
mod io;
pub mod merge;
mod overlap_index;
mod relation_index;
pub mod sort;
pub mod subtract;
pub mod window;

pub use bed::{Genome, read_genome};

/// Strand relationship required between two BED records.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum StrandFilter {
    /// Ignore strand fields.
    #[default]
    Any,
    /// Require matching strands.
    Same,
    /// Require opposing strands.
    Opposite,
}

/// Execute the product command-line entry point.
///
/// This exists only so the package binary can remain a separate Cargo target;
/// CLI parser types and dispatch details are not public API.
#[doc(hidden)]
#[must_use]
pub fn run_binary() -> std::process::ExitCode {
    cli::run()
}
