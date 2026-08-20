# rsomics-bed

`rsomics-bed` consolidates BED interval algebra into one product binary:

```text
rsomics-bed sort
rsomics-bed merge
rsomics-bed intersect
rsomics-bed subtract
rsomics-bed complement
rsomics-bed cluster
rsomics-bed window
rsomics-bed closest
```

Each operation accepts BED input from a file or standard input and writes BED
to standard output unless `--output` is supplied. Run `rsomics-bed --help` or
`rsomics-bed <subcommand> --help` for the complete command tree.

Named outputs are written to a same-directory temporary file and atomically
persisted only after the operation succeeds. Existing outputs therefore remain
unchanged on parse, compatibility, or I/O failure. Paths that identify an input
through spelling, normalization, symlink, or hard link are rejected. A new
output uses normal `0666 & !umask` creation permissions; replacing an existing
output preserves its permission bits.

`--json` uses the `rsomics-common` result and error envelope. Because BED is
itself a stdout format, JSON mode requires a named `--output`: BED goes to that
file and the success envelope goes to stdout. Invalid input, configuration, and
I/O failures use the stable common exit codes 1, 2, and 4 respectively.

## Current slice

- `sort` performs a stable lexicographic chromosome/start sort while preserving
  all original BED columns. bedtools does not preserve input order for every
  exact coordinate tie; stable tie ordering is an intentional product
  guarantee, so compatibility checks compare coordinate order and record
  multiplicity for that case rather than claiming byte-identical tie order.
- `merge` requires grouped chromosomes and nondecreasing starts and emits BED3
  merged spans. A zero-length feature at coordinate zero is rejected when it
  would define a new cluster's negative widened start; it remains valid when a
  preceding interval at the same start already fixes the cluster start at zero.
- `intersect` emits clipped A intervals with A's trailing columns.
- `subtract` emits the uncovered fragments of A with A's trailing columns.
- `complement` requires a chromosome-size file and sorted input in that file's
  chromosome order.
- `cluster` appends one-based cluster IDs to grouped, start-sorted BED input.
  It supports a checked maximum gap and either unstranded or same-strand
  clustering; same-strand mode requires valid BED6 records.
- `window` joins A with B records in symmetric or asymmetric neighborhoods.
  Windows may be strand-relative, B eligibility may require the same or
  opposite strand, and pair, any, count, and no-hit reports are available.
- `closest` emits the nearest eligible B record or records for each A. The
  single-B contract includes strand and different-name filters, overlap
  exclusion, unsigned or signed distances, and all/first/last tie selection.

All operations skip blank lines, `#` comments, and every line whose bytes begin
with lowercase `track` or `browser`, matching bedtools 2.31.1. Consequently,
chromosome names with either lowercase prefix are interpreted as header lines.

`rsomics-intervals 0.3` supplies the validated half-open coordinate model used
by the BED parser. `intersect`, `window`, and `closest` own a private checked
COITrees core because its coordinate limit and BED zero-length policy are not
generic interval contracts. Coordinates above the backend's last safe
inclusive coordinate, `i32::MAX - 1`, are reported instead of truncated or
panicked. Relation operations retain B bytes in one contiguous buffer and add
only the file and coordinate orderings their output contracts consume.
`subtract` builds a separate merged coverage map and supports the full nonzero
`u64` coordinate domain without constructing an unused overlap tree. All
operations reject zero-length widening outside the supported `u64` domain.

`intersect` and `subtract` also intentionally reject a zero-length A or B
record at `0-0`. bedtools 2.31.1 accepts such records, but its virtual interval
would begin at `-1`, which cannot be represented safely by this product's
nonnegative `u64` model. This is a fail-loud compatibility divergence, not an
undocumented omission.

`intersect` indexes each distinct non-empty B coordinate pair once and retains
the original B record IDs separately, so it emits hits in B-file order and
preserves duplicates without a chromosome-wide scan. `subtract` additionally
pre-merges B's virtual coverage per chromosome and queries those disjoint
spans, preventing dense overlapping B records from expanding into a
per-A-record candidate list.

`window` deliberately emits multiple hits in stable B-file order. The default
BEDTools executable instead exposes its internal UCSC-bin traversal order at
some bin boundaries. Hit membership and multiplicity remain compatible, but
rsomics does not reproduce that incidental order instability. `closest` also
accepts arbitrarily ordered B records and applies its documented tie ordering
after eligibility and distance selection.

## Compatibility and performance

The test suite requires the real `bedtools v2.31.1` binary; absence or a
different version is a test failure. CI builds the pinned upstream release
archive after verifying its SHA-256 digest on native Linux and macOS runners for
both `x86_64` and `aarch64`.

`cargo bench --bench operations` compares all eight operations with bedtools on
deterministic fixtures containing real overlap, output, duplicates, merge and
cluster groups, nearest-neighbor ties, and dense relation queries. The
representative release gate uses one million primary relation records plus a
separate dense-count workload and verifies complete output hashes before
timing. Measurements, inputs, flags, and machine provenance are recorded in
[PERFORMANCE.md](PERFORMANCE.md). A new release must retain a strict throughput
or resource-use advantage on its declared hot paths.

## Origin

This product consolidates team-owned historical rsomics implementations listed
in [MIGRATION.md](MIGRATION.md). Compatibility behavior is checked against
bedtools 2.31.1 and committed golden output.

bedtools and COITrees are distributed under the MIT license. This independent
Rust product is licensed under MIT OR Apache-2.0.
