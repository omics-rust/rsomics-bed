# Historical implementation migration

## 0.6.0 GitHub delivery: previous-record spacing

The release adds `spacing` as the fifteenth command; published v0.5.0 and its
fourteen-command archives remain immutable. Historical
`rsomics-bed-spacing@fb4ede5c9a18bfdffcaafd1e45e369d7b6402e95` is retained as
fixture/reference only (classification 3), not merged production code.
BEDTools 2.31.1 supplies the observed previous-record, signed widened-point and
binary32 `1e-9` overlap semantics. Raw start ties are allowed; chromosome blocks
cannot revisit. Optional genome order checks membership/order, not clipping.
Initial headers, including header-only input, and opaque suffixes are retained.

`spacing::{spacing, SpacingOptions}` provides the streaming library interface;
the unified CLI retains stdin, JSON/named output, alias rejection, transactional
output and typed errors. Identical signed point bounds are shared privately by
map, coordinate and spacing; cluster and spacing share the raw-order validator,
with optional genome ranking for spacing only. No Layer A API or dependency is
added, and other operations retain their order policies.

Local formatting, strict Clippy and one clean whole-feature review pass.
[GitHub v0.6.0](https://github.com/omics-rust/rsomics-bed/releases/tag/v0.6.0),
release432289022, was published at2026-11-03T04:53:07Z from
`35d656f00f0db19877b16aec23453db14520b8aa`. Exact-head
[CI41577594195](https://github.com/omics-rust/rsomics-bed/actions/runs/41577594195)
passes all four native targets:128 tests per debug/release/packaged profile,
18 installed calls/fifteen operations per host. Fresh anonymous public Mac ARM
installation passes17 successes and one expected error preserving prior output;
the public source crate matches287 Git bodies/modes and289 members without
compilation. The
[public receipt](https://github.com/omics-rust/rsomics-world/blob/main/.autopilot/oracles/bed-public-delivery-2026-11-03-v06.json)
binds actual delivery. PERFORMANCE.md records final-CI spacing throughput,
higher RSS and slower native warmup separately.
Older measurements retain their original source/binary bindings and losses.
Flank, makewindows and reldist remain unimplemented with separate policy holds;
registry authorization remains held without retry.

## 0.5.0 GitHub delivery: bounded coordinate transforms

Fourteen commands add `slop` and `shift` to the immutable twelve-command
v0.4.0 release. `flank` remains unimplemented pending its point-output policy;
this is not completion of the three-operation coordinate wave.

Historical `rsomics-bed-slop@5fc4c997` and `rsomics-bed-shift@06ae5853`
are fixture/reference assets only (classification 3), not production merges.
The implementation uses the existing BED/genome readers and record output,
a private streaming coordinate kernel, and the unified product CLI. No
micro-crate, dependency or public foundation is added. BEDTools 2.31.1 revision
`705ccfdf2c9a77d71560c8adcece0663c2f5e18e` is the behavioral reference.

The library exposes `coordinate::{slop, shift, SlopOptions, ShiftOptions}`.
Slop accepts finite left/right binary32 amounts, fractional and strand-relative
selection; shift accepts finite plus/minus binary32 displacements and fractional
selection. Both preserve initial headers when requested. Clipping, truncation,
crossed/collapsed slop spans and expanded-point reporting are operation-specific.
Negative strand selects reversed slop widths or minus shift; missing/unknown
strand uses forward widths or plus shift. Invalid resulting BED coordinates,
missing/empty chromosomes and unsupported amounts return typed `InvalidInput`
errors; parse and I/O failures propagate. Named output retains the existing
transactional replacement contract. No general point-boundary equivalence is
inferred from the five captured slop/shift source goldens.

Feature `446b144ed5e0b6c35be1057a8ee297022c5ff09d` and installed-lifecycle
head `161961c688a5c2ecedfe28e3dc152abd8a56a011` pass native-four CI
41569005243 with 127 tests per debug/release/packaged profile and 17 installed
calls across fourteen operations per host. Root adoption checks 748 packet
pins, all 276 Git bodies/modes and 278 source-crate members. This package still
reports 0.4.0 and is not a formal release. The measured Mac ARM CI artifact's
slop/shift throughput advantage and higher RSS are bound in PERFORMANCE.md.
The 0.5.0 metadata change leaves these hot paths unchanged.

[GitHub v0.5.0](https://github.com/omics-rust/rsomics-bed/releases/tag/v0.5.0),
release 432257319, was published at 2026-11-03T03:50:43Z from
`91f2bf9498341b1d1fd4cf27fd244c2820f73ad4`. Exact-head
[CI 41573159869](https://github.com/omics-rust/rsomics-bed/actions/runs/41573159869)
passes all four native targets with 127 tests per debug/release/packaged profile
and 17 installed calls across fourteen operations per host. Fresh anonymous
public Mac ARM installation passes the same lifecycle (16 successes and one
expected error preserving prior output); the public source crate matches all
276 Git bodies/modes, without compilation. The
[public-delivery receipt](https://github.com/omics-rust/rsomics-world/blob/main/.autopilot/oracles/bed-public-delivery-2026-11-03-v05.json)
records the independently adopted public originals. Registry authorization
remains held without retry. Earlier tags, assets and performance records remain
unchanged; later documentation does not replace packaged documentation or assets.

## 0.4.0 GitHub delivery: genome coverage tracks

The twelve-operation release adds `genomecov` to the eleven operations in
immutable GitHub v0.3.0 at bc40084ccc2d3837384da1377d4d57f0600fcb75.
The older v0.2.0 and v0.3.0 tags/assets are not replaced. `unionbedg`, BAM
coverage and paired-read fragment modes are not part of this release.

`rsomics-bed-genomecov@0087cf10fc15` is retained as fixture and algorithmic
reference only (classification3). The product operation uses the existing
strict BED/genome reader, product-private sparse depth events and atomic output,
not the old all-chromosome loader or permissive strand parser. No micro-crate,
dependency or foundation is added. The private significant-digit formatter is
shared with map without changing its implementation.

BEDTools2.31.1 owns chromosome encounter order after strand selection, binary32
fractions/scaling, point expansion and exact-width BED12 splitting. BED12+ is
treated as an outer span here; other product operations keep their existing
block contracts. The final bedGraph run retains the upstream uncapped-max
quirk. Selected zero-origin points, unknown contigs, spans beyond declared
bounds and missing endpoint strands fail explicitly rather than inheriting
unsafe indexing or implicit strand assumptions. Implementation c6e35 passes
four-native correctness:126 tests per debug/release/packaged profile and a
15-call installed lifecycle across twelve operations. Its accepted full-genome
histogram/zero-inclusive-bedGraph workflow retains exact output and measured
Linux executable identities in PERFORMANCE.md, not measurements of the 0.4.0
release binaries.

GitHub release 432026368 was published at 2026-11-02T22:18:19Z as
[v0.4.0](https://github.com/omics-rust/rsomics-bed/releases/tag/v0.4.0), at
`5e5b709d00c07977ce81deb0f7f0f70fece91c4a`. Exact-head CI 41532198240 passed 126
tests per debug/release/packaged profile and 15 installed calls across twelve
operations on all four native targets. Fresh anonymous public Mac ARM
installation passed the same 15-call lifecycle; the downloaded source crate
matched 266 Git bodies and modes without compilation. The registry sparse index
still lists only 0.1.0; authorization remains held without a retry.

## 0.2 GitHub delivery

The delivery-only change starts from clean owner
`02b85a1a348c271e485cca629dc4e71fa075388a`, tree
`49a916836bc119696a6845d97c03c992b1d44dec` (91 files). Every original source,
test, golden, benchmark, Cargo manifest/lock and license body remains unchanged.
Only product-local delivery/guard/notices and installation/CI documentation
are added or updated; no API, operation or dependency is introduced.

The new native gate exercises the complete pinned BEDTools 2.31.1 oracle on
all four native target classes, verifies a clean exact-head Cargo source
package, installs its binary, and checks one eleven-call installed lifecycle:
version, help, the eight declared operations, and a malformed-input failure.
Stdin sort output is reloaded into a JSON/named-output merge; golden bytes
retain trailing fields, duplicate hits, signed closest/no-hit behavior and
known ordering differences. The failure must preserve a pre-existing output.

Source packages retain the complete Git bodies/modes, original manifest,
normalized Cargo metadata, public library/docs and goldens. Native bundles
contain one executable and the original enabled normal/build dependency terms;
[licenses/provenance.json](licenses/provenance.json) distinguishes shipping,
source/dev-only and upstream-attribution scopes. Clean exact-head native CI
41348074250, complete original adoption, GitHub v0.2.0 and a fresh public
Mac ARM eleven-call lifecycle pass for the immutable delivery commit
`bcbee5f87aeda2e7664c279b32c5a3098635ca9b`. Registry publication remains
separately held and unattempted; later documentation does not replace the
release tag or assets.

The retained performance decision stays bound to `d7b1507b178053a087862255d84a244e4921f192`
and measured binary `9a55972f96e11bf087515b61799575505291d12d67902b4c965b8b9e32b51ff8`.
Closest's higher RSS, stable sort/window ordering differences, merge's narrow
win and the rejected candidate remain unchanged. No new delivery binary is
claimed to have been timed. See [PERFORMANCE.md](PERFORMANCE.md).

## Historical source disposition

The initial product slice was reconstructed from these clean local revisions on
2026-07-30:

| Operation | Historical repository | Revision | Disposition |
|---|---|---|---|
| sort | `rsomics-bed-sort` | `5dce75d5a7e28667134585aae129cb8735af4c53` | parser/output policy refactored into shared internal modules |
| merge | `rsomics-bed-merge` | `c9e39b3e760ea55d08f7c56694ea8cc7b039b18c` | sweep algorithm retained; sortedness and checked arithmetic strengthened |
| intersect | `rsomics-bed-intersect` | `e7dc1e7e462e5bc1c40d58f1515c5e1d84a246f4` | index strategy retained; shared parser and explicit backend guard used |
| subtract | `rsomics-bed-subtract` | `927c257de23bef10ab71e6821c04318e6de72a5c` | index strategy retained; shared parser and explicit backend guard used |
| complement | `rsomics-bed-complement` | `94ea9f5526208f841932b02b21337af10f1c78ac` | genome-order sweep retained; genome and bound validation strengthened |

All five historical worktrees were clean during inspection. Their repositories
remain unchanged.

The 0.2 relation slice used these additional clean historical revisions:

| Operation | Historical repository | Revision | Disposition |
|---|---|---|---|
| cluster | `rsomics-bed-cluster` | `b63b75567ba729c016a4baabbbc3bb28bad0718e` | retained the sorted sweep and oracle cases; replaced its parser, CLI, chromosome allocation, and comments |
| window | `rsomics-bed-window` | `875459ee2f793505d8256d958bb634e36a4ab19a` | retained fixtures and benchmark shape only; rebuilt parsing, indexing, output modes, and CLI |
| closest | `rsomics-bed-closest` | `e85ed1339165d2552f86223190975175cbe4318a` | retained corrected distance, zero-length, tie, and no-hit cases; replaced the permissive parser and full B scan |

## Evidence retained

- sort: unsorted input and committed bedtools-sorted output;
- merge: basic and in-domain zero-length bedtools output, plus a fail-loud
  origin-widening fixture that prevents negative BED output;
- intersect: ordinary, unsorted-B, zero-length-A, malformed, and
  out-of-backend-range fixtures, plus an explicit origin-zero divergence test;
- subtract: ordinary A/B fixtures, exact zero-length boundary behavior, and
  mixed zero/nonzero coverage-union behavior, bedtools comparison, and an
  explicit origin-zero divergence test;
- complement: ordinary, zero-length, unknown-chromosome, inverted, and
  unsorted fixtures.
- cluster: distance, book-end, same-strand ordering, malformed strand, and
  sortedness fixtures;
- window: symmetric, asymmetric, strand-relative, strand-filtered, duplicate,
  zero-length, and all four report-mode fixtures;
- closest: overlap, left/right tie, signed orientation, filtering, placeholder,
  zero-length, and all/first/last fixtures;
- all three relation operations: deterministic generated differentials over
  every declared option family and one-million-record release workloads;
- all five operations: seeded random differentials with duplicates,
  zero-length features, unordered inputs, and multiple chromosomes, plus a
  one-million-record end-to-end output and performance comparison with the
  pinned bedtools 2.31.1 oracle; see `PERFORMANCE.md`.

The old single-command `rsomics-help` specifications were not copied. Clap owns
the product/subcommand command tree and renders nested help directly.

## Dependency boundary

`rsomics-intervals 0.3` supplies only the validated, generic interval value
used by the BED parser. BED parsing, headers, column preservation, indexing,
zero-length policy, sortedness, chromosome-size validation, and output
formatting remain private to this product. `intersect`, `window`, and `closest`
share one checked product-private COITrees core; the relation wrapper adds
contiguous raw B storage and directional record order without exposing BED
policy through a foundation. `subtract` uses a separate merged `u64` coverage
map because it does not need an overlap tree.

Distinct non-empty B coordinate pairs are indexed once and expanded back to
their original record IDs for intersect queries. This retains duplicate
multiplicity and B-file output order without the historical per-A
chromosome-wide scan. Subtract builds a separate per-chromosome union of B's
virtual coverage, so dense distinct overlaps are merged once rather than
expanded and re-sorted for every A record.

Relation B input is retained as one byte buffer rather than one allocated
record object per row. `window` adds file-order range queries; `closest` adds
start/end directional orders and stops each direction after the best eligible
distance is exceeded. This replaces the historical eager strings and per-A
full scans without creating another public crate.

## Intentional strictness

bedtools can emit a negative start when a zero-length `0-0` feature begins a
merge cluster or is complemented. `rsomics-bed` rejects that widening because
its public BED model uses nonnegative `u64` coordinates and every successful
output must be readable by the same parser. A later `0-0` record absorbed into
a cluster whose start is already zero remains valid. This scoped divergence is
covered by fail-loud fixtures; all representable zero-length cases remain
checked against bedtools 2.31.1.

The indexed `intersect` and `subtract` operations also reject any zero-length A
or B record at `0-0`. bedtools 2.31.1 accepts origin-zero A records, but their
virtual interval begins at `-1`; neither the product's nonnegative `u64` model
nor the current published index backend can represent that query safely. Live
oracle tests record bedtools' behavior while failure tests lock in the explicit
rsomics-bed divergence.
