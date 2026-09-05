---
name: docs-maintenance
description: Map of which documentation must be updated when crate functionality changes (public API, cargo features, FFI bindings and build linkage, molecule TOML/JSON parsing, GTO evaluation, env variables, CI matrix, examples, claims). Use while or after making such a change, before committing, so readme.md, crate rustdoc, from_toml_docs.md, gto_notation.md, CHANGELOG, and CI stay consistent.
---

## When to use

Any change that touches what the crate **offers** or **promises**: public API
surface (`CInt` methods, `gto`/`parse` modules, `prelude`), supported
integrals/intors, cargo features, FFI bindings or the linked C library,
molecule input format (TOML/JSON), build/link environment variables, examples,
CI matrix, version/compatibility claims. Internal
refactors with identical behavior need nothing from this map except the
verification commands at the end.

## Document inventory

| File | Audience / role |
| --- | --- |
| `readme.md` (root) | project-facing single source of truth: scope, badges, minimal examples, installation, cargo features, env variables, 50-lines-RHF claim |
| `CHANGELOG.md` (root) | per-release notes; prepended by release-plz from conventional commits (`release-plz.toml`, `.github/workflows/release-plz.yml`) |
| `libcint/src/lib.rs` crate docs | overview, "Quick Links", minimal-example doctest (manual mirror of `readme.md`), badge table (manual mirror), RHF example claim |
| `libcint/src/cint*.rs` | rustdoc for the `CInt` API (`integrate`, `integrate_row_major`, `integrate_cross`, `eval_gto`, properties, builder) |
| `libcint/src/parse/from_toml_docs.md` | normative spec of the TOML/JSON molecule input (feature `bse`); rendered into rustdoc on `CIntMol` via `libcint/src/parse/mole.rs` |
| `libcint/src/gto/gto_notation.md` | mathematical notation for GTO/grid evaluation; rendered as the `gto` module docs (`libcint/src/gto/mod.rs`) |
| `libcint/src/ffi/` | FFI bindings and supported-integrator lists (`cint_ffi`/`cint_wrapper`, `cecp_ffi`/`cecp_wrapper`) |
| `libcint/examples/h2o_rhf.rs` | the RHF example; `include_str!`-ed into crate docs, cited by `readme.md` line counts |
| `libcint/scripts/` | PySCF reference/comparison notebooks and FFI header generators (`generate_libcint_header.py`, `generate_cecp_header.py`) |
| `Cargo.toml`, `libcint/Cargo.toml`, `libcint-src/Cargo.toml` | workspace metadata, feature list + comments |
| `.github/workflows/` | CI legs: clippy, test matrix (build-from-source, custom build, static, qcint), release-plz |
| `katex-header.html` (root) | rustdoc KaTeX header for equations |
| `libcint/tests/`, `libcint/src/test_mol.rs` | behavioral ground truth and test-molecule fixtures |

Do not document in `tmp*` (scratch, gitignored) or agent memory.

## Change-type → update map

**New/changed public method on `CInt` (or `CIntMol`)**:
1. rustdoc at the definition (usually `libcint/src/cint.rs` or the relevant
   `cint_*.rs`/`parse/` module).
2. `libcint/src/lib.rs` "Quick Links to Important Parts" if headline-worthy.
3. `readme.md` minimal example only if the headline usage itself changes.
4. `libcint/tests/` coverage (+ `test_mol.rs` fixture if a new molecule helps).

**New intor/integral or FFI surface change**: bindings in
`libcint/src/ffi/cint_ffi.rs` (regenerate via
`libcint/scripts/generate_libcint_header.py` if the C header changed);
supported-integrator list in `ffi/cint_wrapper.rs`; feature gates
(`with_f12`, `with_4c1e`) where applicable; ECP changes go through
`cecp_ffi.rs`/`generate_cecp_header.py`. Badge-table versions in `readme.md`
**and** `libcint/src/lib.rs` if the bound libcint/qcint/PySCF version moves.

**Cargo feature add/change**: `libcint/Cargo.toml` feature list + comment;
`readme.md` "Cargo features" section; CI leg in
`.github/workflows/test-libcint.yml` if the feature needs one;
feature-gated tests.

**Molecule input (TOML/JSON) field change**: `from_toml_docs.md` (the spec —
update first); rustdoc on `CIntMol::from_toml`/`from_json` in
`libcint/src/parse/`; `readme.md` "Building molecule from TOML/JSON" section
(field list); parse-module tests.

**Build/linkage/env-variable change** (new env var, static/dynamic, qcint
switch, source mirror): `libcint/build.rs` and `libcint-src/build.rs`;
`readme.md` "Installation and Cargo Features" + "Shell environment
variables"; `CLAUDE.md` build section if the local dev workflow changes;
the custom-build CI leg.

**Memory layout / shape convention change** (column-major vs row-major,
component axes): `libcint/src/lib.rs` crate-doc layout statements; shape
comments in the `readme.md` minimal example; tests asserting shapes.

**RHF example change**: `libcint/examples/h2o_rhf.rs` (also rendered in
docs); recount the line numbers in `readme.md` and `libcint/src/lib.rs`
("code with core algorithms is N lines"); keep
`libcint/scripts/pyscf_h2o_rhf.py` as the comparison and update its cited
line count too.

**GTO evaluation change** (grids, `blksize`/`NLANE`, derivatives):
`libcint/src/gto/gto_notation.md` if notation moves; rustdoc in
`gto/mod.rs`/`grid_ao_drv.rs`; `readme.md` only if user-visible.

**Version/release change**: workspace version in root `Cargo.toml`
(crates inherit); `CHANGELOG.md` — normally handled by release-plz from
conventional commits, hand-write only for the initial entry;
`readme.md` version strings in the installation snippets.

## Cross-file invariants

1. **`readme.md` minimal example ↔ `libcint/src/lib.rs` doctest are a manual
   mirror** — deliberately duplicated (not `include_str!`). Any edit to one
   must be copied verbatim into the other; verify with `cargo test --doc`.
   The badge table is likewise duplicated in both files.
2. **Symlinked documents — never replace with real files**:
   `libcint/readme.md` and `libcint-src/readme.md` → `../readme.md`;
   `libcint/CHANGELOG.md` → `../CHANGELOG.md`;
   `libcint/katex-header.html` → `../katex-header.html`. Edit the root file.
   (Agent setup follows the same pattern: `AGENTS.md` → `CLAUDE.md`,
   `.claude` → `.agents`.)
3. **Badge/version claims** (libcint, qcint, PySCF-ECP tags) must match what
   `libcint-src` actually downloads/binds, in both `readme.md` and the
   `lib.rs` badge table.
4. **No personal or absolute paths** in tracked docs (repo-wide AI-agent
   notice in `CLAUDE.md`).

## Verify after doc edits

```sh
cargo doc --no-deps -p libcint   # expect 0 warnings (KaTeX header resolves via symlink)
cargo test --doc -p libcint      # readme-mirrored doctest runs
cargo test -p libcint            # full suite; needs CINT_DIR/LD_LIBRARY_PATH/CINT_DEV (see CLAUDE.md)
cargo clippy --workspace --all-targets
```
