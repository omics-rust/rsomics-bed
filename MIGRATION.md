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
- all five operations: seeded random differentials with duplicates,
  zero-length features, unordered inputs, and multiple chromosomes, plus a
  one-million-record end-to-end output and performance comparison with the
  pinned bedtools 2.31.1 oracle; see `PERFORMANCE.md`.

The old single-command `rsomics-help` specifications were not copied. Clap owns
the product/subcommand command tree and renders nested help directly.

## Dependency boundary

`rsomics-intervals 0.2.0` is consumed only for its generic interval model and
overlap index. BED parsing, headers, column preservation, zero-length policy,
sortedness, chromosome-size validation, and output formatting remain private to
this product. Because the published backend currently exposes infallible calls,
the product temporarily rejects indexed intervals beyond its last safe
inclusive coordinate (`i32::MAX - 1`) before build or query. A future published
release with recoverable checked build/query APIs should replace this local
guard.

Distinct non-empty B coordinate pairs are indexed once and expanded back to
their original record IDs for intersect queries. This retains duplicate
multiplicity and B-file output order without the historical per-A
chromosome-wide scan. Subtract builds a separate per-chromosome union of B's
virtual coverage, so dense distinct overlaps are merged once rather than
expanded and re-sorted for every A record.

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
