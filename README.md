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
rsomics-bed coverage
rsomics-bed annotate
rsomics-bed map
rsomics-bed genomecov
```

Query input accepts a file or standard input; secondary annotation inputs
require files. Each operation writes BED
to standard output unless `--output` is supplied. Run `rsomics-bed --help` or
`rsomics-bed <subcommand> --help` for the complete command tree.

## Installation and delivery status

The twelve-command product is available as [GitHub v0.4.0](https://github.com/omics-rust/rsomics-bed/releases/tag/v0.4.0),
adding `genomecov`. Release 432026368 was published at 2026-11-02T22:18:19Z from
`5e5b709d00c07977ce81deb0f7f0f70fece91c4a`.
The registry sparse index still lists only the earlier five-command 0.1.0;
registry authorization remains held and has not been retried.
The older [v0.2.0](https://github.com/omics-rust/rsomics-bed/releases/tag/v0.2.0)
and [v0.3.0](https://github.com/omics-rust/rsomics-bed/releases/tag/v0.3.0)
tags and assets remain unchanged. Published delivery guarantees refer to the
release tag rather than the current checkout.

The four native Linux/macOS x86_64/aarch64 archives contain one `rsomics-bed`
executable, documentation, original license/notice bodies, and `BUILD.json`,
with a `SHA256SUMS` companion. Exact-head CI 41532198240 passed 126 tests per
debug/release/packaged profile and the 15-call, twelve-operation installed
lifecycle on all four native targets. A fresh anonymous public Mac ARM download
passed the same lifecycle; its source crate matched all 266 Git bodies and modes
without another source compilation.
BUILD records actual native host/toolchain, source head, lock, binary and
package identities; it is not a signing, minimum-OS, or cross-build guarantee.

To build the reviewed source checkout on this project's local Mac, use
`scripts/guarded-cargo.sh install --path . --locked --root` with an explicitly
external installation directory. The guard verifies mounted external paths
and rejects Cargo storage overrides; its trusted `--exec` mode is not a sandbox
for arbitrary tools or their output arguments. Other users may build the
complete Cargo source package with `cargo install --path . --locked` in their
own environment. [MIGRATION.md](MIGRATION.md) records the delivery boundary.

Named outputs are written to a same-directory temporary file and atomically
replaced only after the operation succeeds. Parse, compatibility, and I/O
failures before that replacement leave an existing output unchanged. A later
directory-sync error can return a nonzero exit after the new file has already
replaced it; failure does not imply rollback after commit. Paths that identify
an input through spelling, normalization, symlink, or hard link are rejected. A new
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

- `coverage` reports overlap counts, covered bases/fractions or depth
  histograms, with strand, block and fraction selection.
- `annotate` appends independent counts, fractions or both for multiple B
  files; files are not pooled into one annotation.
- `map` reduces selected columns of sorted B records overlapping each sorted
  A record. Its23 operations include numeric summaries, counts and text lists.
  Columns and operations broadcast when either list has one item. Empty B is
  rejected; selected invalid/non-finite numeric values are errors.
- `genomecov` reports genome-wide depth histograms, nonzero or zero-inclusive
  bedGraph runs, or per-base depth against a required genome file. It supports
  BED6 strand/endpoint selection, exact-width BED12 blocks, depth scaling and
  histogram depth bins, using sparse events rather than chromosome-sized arrays.

`genomecov` requires chromosome-grouped BED, but starts within a group need not
be sorted. Selected chromosome groups retain encounter order; untouched
chromosomes follow genome-file order. Under `--split`, BED12+ retains its outer
span rather than using BED12 blocks. Selected origin points and spans outside
the declared genome fail explicitly. `--max` retains BEDTools' uncapped final
bedGraph-run behavior. BAM input and paired-read fragment modes are excluded.

`coverage` and `map` currently reject split-block selection combined with an
explicit fraction, reciprocal or either-fraction flag. Ordinary split-block
selection is supported; the combined aggregate-fraction policy is not claimed.
`map -o` names the output file, unlike BEDTools where it selects reductions;
use `--operations` for reductions.

Sorted `map` preserves boundary point hits independent of equal-start A/B
record order. BEDTools2.31.1 can lose these hits; the two deliberate repairs
are `BEDMAP-2311-POINT-CACHE-LOSS` and
`BEDMAP-2311-POINT-LOOKAHEAD-LOSS`. Map uses private widened signed virtual
bounds, so a point at0-0 is supported without emitting negative coordinates.

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
indexed operations retain their documented nonnegative virtual-bound limits;
the map-private signed virtual model does not change those contracts.

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

Those historical measurements belong to source
`d7b1507b178053a087862255d84a244e4921f192` and executable SHA-256
`9a55972f96e11bf087515b61799575505291d12d67902b4c965b8b9e32b51ff8`,
not to a new delivery binary. Closest retains its RSS disadvantage, merge's
1.11x margin is narrow, and the rejected closest candidate remains rejected.
Delivery verification is not another benchmark or a universal memory claim.

New full RefSeq-derived workflows on Apple M2 retain complete output
comparisons. Coverage summary/count/histogram improves wall and process RSS;
dual-file annotate both-mode has median native/source wall0.392/37.403s and
RSS56.9/175.0MB. Sorted map has median wall0.577/0.470s and RSS4.31/4.98MB:
it is about23percent slower but uses13.5percent less peak process memory.
Map qualifies a resource-use advantage, not a throughput advantage. See
[PERFORMANCE.md](PERFORMANCE.md) for source and executable bindings and limits.

The full-genome `genomecov` campaign on shared Linux Xeon Gold 6238R uses all
584,846 exon records and the complete 93-contig hg19 genome. Histogram median
native/source wall is 0.21/113.77s; zero-inclusive bedGraph is 0.36/20.80s.
Measured maximum-process RSS is 5.51/1,998.32MB in both modes. These are scoped
full-workflow results bound to the c6e35 implementation and measured Linux
executables, not measurements of the 0.4.0 release binaries or other modes.

## Origin

This product consolidates team-owned historical rsomics implementations listed
in [MIGRATION.md](MIGRATION.md). Compatibility behavior is checked against
bedtools 2.31.1 and committed golden output.

The original upstream and dependency terms, identities and scope are in
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md),
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md), and [licenses](licenses/).
No BEDTools oracle executable is redistributed. This independent Rust product
is licensed under MIT OR Apache-2.0.
