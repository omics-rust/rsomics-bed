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

## 2026-08-20 0.2 candidate Linux gate

The eight-command candidate passes its representative throughput gate. The
exact performance source is
`d7b1507b178053a087862255d84a244e4921f192`; exact-head CI run `32328054861`
passed formatting, strict Clippy, debug and release tests, package verification,
the benchmark smoke, and the pinned live oracle on native Linux and macOS for
both `x86_64` and `aarch64`.

Source and tools:

- release binary SHA-256
  `9a55972f96e11bf087515b61799575505291d12d67902b4c965b8b9e32b51ff8`;
- `rustc 1.95.0 (59807616e 2026-04-14)` and GCC 12.3.0 on the performance
  host; CI independently used the declared minimum Rust 1.91;
- registry `rsomics-common 0.11.0`, `rsomics-help 0.4.0`, and
  `rsomics-intervals 0.3.0`, locked to their checksums;
- bedtools `v2.31.1`, built from the release archive with SHA-256
  `fc7e660c2279b1e008b80aca0165a4a157daf4994d08a533ee925d73ce732b97`;
- bedtools binary SHA-256
  `287fb59cd3f68f43e45df5e807596fdc2f170be38c9b8953aaceb511b20cbfb0`;
- base runner SHA-256
  `ba09fce09fb7152e1a02d983defa282a951d64931b8a9ef4d8de66e64d48f6f2`;
- relation runner SHA-256
  `9d16cbc2b808c93210b5a3b4c637a08fd2a84b6e33545e1a73b6324a9d004e76`.

The host was `dell-Precision-7920-Tower`, Linux 6.8.0-90-generic,
`x86_64`, with two Intel Xeon Gold 6238R CPUs. Commands were bound to physical
cores 48-51 on one socket after those cores and their sibling threads measured
94-100% idle. GNU `time` ran one warmup and ten measurements. The relation
runner alternated implementation order for each pair and recorded wall time,
user plus system CPU time, and maximum RSS; tables report mean and sample
standard deviation for time and median RSS.

The original five-operation one-million-record fixture and every output hash
were unchanged from the 2026-07-31 gate. Re-running the complete recipe after
the shared index refactor produced:

| Operation | rsomics-bed | bedtools 2.31.1 | Speedup | RSS KiB, ours / bedtools | Decision |
|---|---:|---:|---:|---:|---|
| sort | 0.398 ± 0.004 s | 1.097 ± 0.011 s | 2.76x | 209,664 / 360,490 | pass |
| merge | 0.204 ± 0.005 s | 0.226 ± 0.005 s | 1.11x | 2,688 / 4,480 | pass, narrow margin |
| intersect | 0.699 ± 0.006 s | 2.698 ± 0.026 s | 3.86x | 99,678 / 344,960 | pass |
| subtract | 0.602 ± 0.008 s | 3.076 ± 0.022 s | 5.11x | 19,712 / 344,960 | pass |
| complement | 0.243 ± 0.005 s | 0.308 ± 0.004 s | 1.27x | 5,376 / 4,480 | throughput pass |

The relation fixture contains one million A and cluster records across ten
chromosomes. Window B contains one matching record per A plus a duplicate
every fiftieth record. Closest B contains equal-distance left and right records
plus an overlap every fiftieth record. The dense count lane contains 5,000 A
and 5,000 mutually overlapping B records.

| Input | Bytes | SHA-256 |
|---|---:|---|
| cluster BED6 | 33,433,300 | `108c4923808666f7a8942a1699f6f3fffff91c32a211d6e11cd9cceb9b705333` |
| relation A BED6 | 34,877,780 | `dcd79d7e23a118d4d2033d6325a4f9dbc672c87168948f364bdf4fe5d807bc60` |
| closest B BED6 | 72,473,070 | `b56cdd6494a76c497df84bcaf63a93a3617178b2b39730519f452175b37acedc` |
| window B BED6 | 35,575,290 | `0ec499733c02882c82aafaf04a970af6e4e1eb393907b800fd498fff298bd9a5` |
| dense A | 113,890 | `6f238adc556453f0fae3540a06e3f2fee21d38f9d28bfb6fd21ccaca2b60e556` |
| dense B | 107,780 | `726c8bf6160d2a41bf3bea165bf6817b307af0a4b6139fa42ad534c0c97982c5` |

Every complete relation output matched bedtools byte for byte before timing:

| Operation | Output lines | Output bytes | SHA-256 |
|---|---:|---:|---|
| cluster | 1,000,000 | 39,877,775 | `2fdc8c4160e24dcc0349df758ed9a320419383b6b084f9f6d453c6c7748e2101` |
| cluster same-strand | 1,000,000 | 40,322,196 | `b7189c6902ecf6c8a6d76c5ce6027d46025c2411acf2ce12b52228be170872b4` |
| window pairs | 1,020,000 | 71,150,580 | `c8b9ce3929dd859573089a9e198facdffd91bdf49fe25182a81a9c6ef9d58453` |
| closest | 1,980,000 | 140,096,100 | `1f2f82e843e97431220560c6d1a0f4e161203a1569afb27b92144fc46d2fc82f` |
| closest unsigned distance | 1,980,000 | 146,016,100 | `199b21ac7d3abfc22881ab7f8e525ea12e43be5a2e81569534c34d4d96f55429` |
| dense window count | 5,000 | 138,890 | `6cd3bf156bf32c431e95cfdba6bf308ec1228cbf8b92cccaea47d9d11354eb85` |

| Operation | rsomics-bed | bedtools 2.31.1 | Wall / CPU speedup | RSS KiB, ours / bedtools | Decision |
|---|---:|---:|---:|---:|---|
| cluster | 0.229 ± 0.006 s | 0.932 ± 0.010 s | 4.07x / 4.15x | 2,688 / 4,480 | pass |
| cluster same-strand | 0.291 ± 0.003 s | 1.509 ± 0.018 s | 5.19x / 5.33x | 18,816 / 451,372 | pass |
| window pairs | 0.676 ± 0.010 s | 2.820 ± 0.024 s | 4.17x / 4.19x | 158,908 / 470,400 | pass |
| closest | 1.349 ± 0.011 s | 1.662 ± 0.026 s | 1.23x / 1.23x | 316,696 / 28,636 | throughput pass |
| closest unsigned distance | 1.420 ± 0.019 s | 1.927 ± 0.025 s | 1.36x / 1.36x | 316,704 / 13,136 | throughput pass |
| dense window count | 0.157 ± 0.005 s | 2.194 ± 0.020 s | 13.98x / 14.39x | 3,584 / 8,064 | pass |

The first exact benchmark head, `c8e09eeaaba0aa5ce3c20d7f0089cb85c8380264`,
failed the default closest gate at 0.943x and used 577,304 KiB. It was not
accepted. The compact relation-index refactor retains B as one raw byte buffer,
removes duplicate virtual coordinates, and reuses the index's start order. It
reduced closest time from 1.766 to 1.349 seconds and peak RSS from 577,304 to
316,696 KiB on the same fixture.

Closest still uses materially more memory than bedtools because rsomics accepts
arbitrarily ordered B input and retains stable B-file tie order while bedtools
requires coordinated sorted streams. This release claims a measured throughput
advantage for closest, not a memory advantage. Window likewise guarantees
stable B-file hit order instead of reproducing the upstream UCSC-bin traversal
artifact that appears when multiple hits cross internal bin boundaries.

Raw evidence is retained under
`/data3/liangjy/rsomics-linux-x86_64-20260820/bed-gate/results`:

- `base-1m-d7b1507.json`, SHA-256
  `b4fb912d945dd83d1467f40f39ecf534cb8f27641cf8b312b3d12c290c3a0e8e`;
- `relations-1m-d7b1507.json`, SHA-256
  `695c0d7a34020a006251d18e212d8760a1b53a4e59b7c07855ee594c7d8eee14`;
- rejected pre-optimization `relations-1m.json`, SHA-256
  `f66f1248f30e084b1bb78c51c6147ef636097dcf21f41d311cba4800f4486a6a`.

## 2026-07-31 representative Linux gate

Source and tools:

- `rsomics-bed` Git revision
  `e8898dbcb0db7a398fbd84f7627ad2be322f475f`;
- release binary SHA-256
  `e78a4d67091ffa121d12d451c0e10526ca7effccbebd73afb7ba7ffcddbc009f`;
- `rustc 1.91.0 (f8297e351 2025-10-28)`;
- crates.io `rsomics-common 0.7.0`, `rsomics-help 0.4.0`, and
  `rsomics-intervals 0.3.0`, locked to registry checksums;
- bedtools `v2.31.1`, built from the release archive whose SHA-256 is
  `fc7e660c2279b1e008b80aca0165a4a157daf4994d08a533ee925d73ce732b97`;
- bedtools binary SHA-256
  `40c465f999a58dc8ff42959ff8635ede4f54dce51a4309e0fa43afa7c9a38e56`.

The host was `dell-Precision-7920-Tower`, Ubuntu 22.04, Linux
6.8.0-90-generic, `x86_64`, with two Intel Xeon Gold 6238R CPUs. Commands were
bound to physical cores 48-51 on one socket. Source, target output, temporary
files, fixtures, and results remained under
`/data1/liangjy/rsomics-linux-x86_64-20260731/bed-gate`. The Rust toolchain and
Cargo registry cache were reused from the preceding external-disk gate.

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
| sort | 0.390 ± 0.000 s | 1.077 ± 0.016 s | 2.76x | 197,120 / 360,624 |
| merge | 0.197 ± 0.007 s | 0.228 ± 0.004 s | 1.16x | 1,792 / 4,480 |
| intersect | 0.672 ± 0.004 s | 2.683 ± 0.116 s | 3.99x | 114,700 / 344,960 |
| subtract | 0.592 ± 0.009 s | 3.011 ± 0.023 s | 5.09x | 18,816 / 344,960 |
| complement | 0.236 ± 0.007 s | 0.302 ± 0.015 s | 1.28x | 4,480 / 4,480 |

All five operations pass the throughput gate on this representative workload.
Sort, merge, intersect, and subtract also use less peak memory. Complement has
equal measured peak memory and a strict throughput advantage.

The raw JSON is
`/data1/liangjy/rsomics-linux-x86_64-20260731/bed-gate/results/representative-1m.json`,
SHA-256
`2133285863696e1b22f10aadbd69e5f80e16af235d8a2ba640cfb1d567070a46`.

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

Code revision `e8898dbcb0db7a398fbd84f7627ad2be322f475f`
passed exact-head CI run `30598213193` on native Linux and macOS for both
`x86_64` and `aarch64`.

- native Linux `aarch64` has correctness CI but no representative performance
  host;
- the current public API and production hot paths passed release review.
