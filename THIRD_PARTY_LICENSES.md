# Upstream behavior and implementation provenance

## BEDTools

BEDTools 2.31.1 is the external behavior oracle for the eight declared BED
operations and their committed golden fixtures. The pinned source revision is
`705ccfdf2c9a77d71560c8adcece0663c2f5e18e` (`v2.31.1`). The original root
`LICENSE` body is retained unchanged as
[licenses/bedtools-2.31.1-MIT.txt](licenses/bedtools-2.31.1-MIT.txt):

Copyright (c) 2009-2019 Aaron Quinlan.

Acquisition used the exact-revision
[upstream LICENSE](https://raw.githubusercontent.com/arq5x/bedtools2/705ccfdf2c9a77d71560c8adcece0663c2f5e18e/LICENSE),
SHA-256 `0e44360fb7c5f08f49575e8de98606a10a9bc3e785e2c65cbef63e676190c86d`.
The CI executable is built separately from the checksum-pinned release
archive; neither that executable nor its source is redistributed in the
native BED bundle. This root notice is not an assertion that every upstream
component has identical terms or permission to copy future components.

## COITrees and other locked Rust dependencies

`coitrees 0.4.0` supplies the interval-tree backend used by the product-private
checked relation index. Its original MIT terms and copyright are included at
[licenses/registry/coitrees-0.4.0/LICENSE](licenses/registry/coitrees-0.4.0/LICENSE):

Copyright © 2019 Daniel C. Jones (dcjones@cs.washington.edu).

The backend is a dependency, not a new rsomics public interval policy API.
The remaining original locked dependency terms and graph scopes are recorded
in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md) and
[licenses/provenance.json](licenses/provenance.json).

## Historical team-owned assets

The retired operation-sized rsomics repositories are team-owned migration
assets, not current standalone installation targets or third-party packages.
[MIGRATION.md](MIGRATION.md) records the exact source revisions, selected
algorithms, fixtures and benchmark recipes for sort, merge, intersect,
subtract, complement, cluster, window and closest. Their provenance does not
extend compatibility or performance beyond the declared product contracts.

The product's original [MIT](LICENSE-MIT) and
[Apache-2.0](LICENSE-APACHE-2.0) license bodies are unchanged.
