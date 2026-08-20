# Historical implementation migration

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
