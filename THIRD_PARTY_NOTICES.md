# Enabled Rust dependency notices

Native archives include original terms for the unchanged locked normal/build
dependency graphs of the four supported Linux/macOS targets, including build
dependencies and procedural macros. [licenses/provenance.json](licenses/provenance.json)
records each package, declared expression or license file, original crate
checksum, exact copied body identity, and per-target graph. All 42 packages
in the shipping union were checked against their cached checksum-locked
registry archives.

Original texts are under [licenses/registry](licenses/registry/). Applicable
additional terms, including unicode-ident's Unicode data terms, rustix and
linux-raw-sys copyright/LLVM exceptions, and memchr/same-file alternative
terms, are retained without normalization or replacement by templates.

The provenance also distinguishes 32 packages present only in the broader
resolved source/dev/other-target graph. Their code is not included in the
native shipping graph; the Cargo source package retains the unchanged complete
manifest and lock for rebuilding and testing. This is not an inventory of
linked system libraries or a legal-clearance claim.

Upstream behavior and historical team-owned code provenance are described in
[THIRD_PARTY_LICENSES.md](THIRD_PARTY_LICENSES.md). BEDTools executables,
compiler runtimes, build scratch and compatibility-oracle source trees are
not bundled in the native archive.
