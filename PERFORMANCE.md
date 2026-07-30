# Performance evidence

The release benchmark is `cargo bench --bench operations`. It compares every
current `rsomics-bed` operation with bedtools 2.31.1 on deterministic
multi-chromosome fixtures while discarding timed output.

The fixture is designed to exercise work rather than process startup alone:

- sort receives reverse-coordinate BED4 records on ten chromosomes;
- merge receives groups of five overlapping intervals;
- every A record in intersect has a B hit, and every fiftieth B interval is
  duplicated;
- subtract emits two fragments per A record;
- complement emits the real gaps between merged groups.

Results are recorded only after the complete output of every operation matches
the pinned oracle and the full test suite passes. Each release record includes
the exact source and binary identities, fixture hashes, command flags, timing
distribution, peak-memory method, and a per-operation decision.

## 2026-07-30 representative Linux gate

Source and tools:

- `rsomics-bed` Git revision
  `ed415eeebd9d6a3bcb34cc9cf15bcfc5f7c587cd`;
- release binary SHA-256
  `8f4808840fa5a8cab88079aa617f60cef2c1196fb5dbf0edd2dd8a3f9a29ae17`;
- `rustc 1.91.0 (f8297e351 2025-10-28)`;
- bedtools `v2.31.1`, built from the release archive whose SHA-256 is
  `fc7e660c2279b1e008b80aca0165a4a157daf4994d08a533ee925d73ce732b97`;
- bedtools binary SHA-256
  `40c465f999a58dc8ff42959ff8635ede4f54dce51a4309e0fa43afa7c9a38e56`.

The host was `dell-Precision-7920-Tower`, Ubuntu 22.04, Linux
6.8.0-90-generic, `x86_64`, with two Intel Xeon Gold 6238R CPUs. Commands were
bound to physical cores 48-51 on one socket. Source, Cargo home, target output,
temporary files, fixtures, and results remained under
`/data1/liangjy/rsomics-linux-x86_64-20260730/bed-gate`.

The formal fixture has 1,000,000 A or primary records across ten chromosomes.
B has 1,020,000 records because of the deliberate duplicates. Fixture
SHA-256 values are:

| Input | Bytes | SHA-256 |
|---|---:|---|
| reverse sort | 30,766,670 | `60a5feaa6296923bbca5045df9b8b5b93b60f4ecb44e79c62b5b42a3eb6e3e7f` |
| merge/complement | 29,877,790 | `e62523d7f64e7a8d40e81da71cd45cd7842322c5d80577b3acffcf8e1da2632b` |
| A | 30,766,670 | `bbdf49df43933b13171b201c0a00c6b71b0f301ac104f189949f8fe4665aa3d2` |
| B | 31,381,970 | `b90a8954bfe3e358f679882dbd9ee4d766b436bcdced98e13a534df45143ee46` |
| genome | 150 | `36fa2190ac75b78680667f540c8db147444fbb4426c9555cfe951584c5027095` |

Before timing, each complete output matched bedtools byte for byte:

| Operation | Output lines | Output bytes | SHA-256 |
|---|---:|---:|---|
| sort | 1,000,000 | 30,766,670 | `ff81f1cb965c57d8dfe282a76382e54579b5fdf7f24929837a71220ac1a8c209` |
| merge | 200,000 | 4,177,770 | `c26e988fae5a0ba8e6835117609ee65d2983d5b32ee24f8b9396e65d20d1b038` |
| intersect | 1,020,000 | 31,381,970 | `40a55d578b0b58e58c3c124e39850d31ac9a551c5a4f77a491891707e37c94b6` |
| subtract | 2,000,000 | 61,533,350 | `6bbd97179ccfbd6a6d545b77f3c56cdea1a7661ad60b4dbcbef01a233c361772` |
| complement | 200,010 | 4,177,940 | `f66397a8bdd2109271c667715f0e7a18a82f7cfe99618475c883216b936964fc` |

GNU `time` ran one warmup followed by ten measurements for each command.
Elapsed values are mean and sample standard deviation; peak RSS is the median
of the ten process high-water marks.

| Operation | rsomics-bed | bedtools 2.31.1 | Speedup | RSS KiB, ours / bedtools |
|---|---:|---:|---:|---:|
| sort | 0.409 ± 0.006 s | 1.071 ± 0.014 s | 2.62x | 198,016 / 360,624 |
| merge | 0.198 ± 0.004 s | 0.225 ± 0.005 s | 1.14x | 2,688 / 4,480 |
| intersect | 1.002 ± 0.008 s | 2.615 ± 0.011 s | 2.61x | 261,140 / 344,960 |
| subtract | 0.960 ± 0.019 s | 2.992 ± 0.017 s | 3.12x | 289,946 / 344,960 |
| complement | 0.239 ± 0.006 s | 0.297 ± 0.007 s | 1.24x | 5,376 / 4,480 |

All five operations pass the throughput gate on this representative workload.
Sort, merge, intersect, and subtract also use less peak memory. Complement uses
about 20% more peak memory but retains a strict throughput advantage.

The raw JSON is
`/data1/liangjy/rsomics-linux-x86_64-20260730/bed-gate/results/representative-1m.json`,
SHA-256
`277aa73ccd944cc00331f9d7b111467fca5b86ba123ab2a5f5ef187649a6bbbc`.

## Dense-subtract scaling

The Criterion benchmark also contains an adversarial scaling lane. For each
`n`, A contains `n` copies of `[10000,10001)` and B contains the `n` distinct
intervals `[i,20000+i)`. Every A overlaps every raw B record.

The original Apple M2 measurements remain useful for the scaling shape:

| Dense records | rsomics-bed | bedtools 2.31.1 | Speedup |
|---:|---:|---:|---:|
| 500 | 2.800 ms | 8.024 ms | 2.87x |
| 1,000 | 2.887 ms | 21.804 ms | 7.55x |
| 2,000 | 3.312 ms | 75.057 ms | 22.66x |
| 4,000 | 4.190 ms | 294.88 ms | 70.38x |

This lane supports the decision to pre-merge B coverage for subtraction; it is
not used as a representative end-to-end speedup claim.

## Superseded sparse pilot

An earlier 50,000-record Apple M2 pilot reported 72.20x intersect and 68.43x
subtract speedups. Its A and B fixtures did not overlap, so both commands
produced empty output. Those figures describe a sparse no-hit case only and
are superseded by the representative Linux gate above. They are not release
claims.

## Remaining release gates

Evidence revision `76d02dbc9c0fd549782f1e68e2b0ef5e64f13d45`
passed exact-head CI run `30556922040` on native Linux and macOS for both
`x86_64` and `aarch64`.

- native Linux `aarch64` has correctness CI but no representative performance
  host;
- the final publication public API review remains.
