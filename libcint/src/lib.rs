#![allow(clippy::excessive_precision)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::approx_constant)]

//! This project contains libcint (C language) FFI bindings, wrapper and
//! build-from-source.
//!
//! # Quick Links to Important Parts
//!
//! - [`CInt`] basic struct for defining molecule for electron integral.
//! - [`CInt::integrate`]: integrate function (partly PySCF's `mol.intor`
//!   counterpart, column-major which is the same to convention of libcint C
//!   library, and the same memory layout to PySCF's 2/3-center integrals).
//! - [`CInt::eval_gto`]: evaluate GTO values on grids, useful for DFT
//!   computation (PySCF's `mol.eval_ao` or `mol.eval_gto` counterpart, with
//!   same memory layout to PySCF, but note the returned shape is in
//!   column-major).
//! - [`CInt::integrate_row_major`]: integrate function with row-major output
//!   (the same shape convention to PySCF, and same memory layout to PySCF
//!   4-center integrals).
//! - [`CInt::integrate_cross`]: integrate function with multiple [`CInt`]s,
//!   useful for integral with auxiliary basis sets.
//! - [`CInt::integrate_with_args`], [`CInt::integrate_args_builder`]: integrate
//!   function with arguments builder (advanced options).
//! - [`CInt::get_cintopt`], [`CIntOpt`]: build a reusable libcint optimizer and
//!   pass it through [`IntegrateArgs::cintopt`] / [`IntorCrossArgs::cintopt`]
//!   instead of rebuilding it on every call. For cross integrals,
//!   [`CInt::get_cintopt_cross`] builds it from the same `mols` list.
//! - [`CIntMol`] (feature `bse`): molecule builder from TOML/JSON, a weaker
//!   version of PySCF's `Mole` class (handling only atom coordinates and basis
//!   sets).
//! - [cint_wrapper](ffi::cint_wrapper) and [cecp_wrapper](ffi::cecp_wrapper):
//!   supported integrators.
//!
//! # Introduction
//!
//! [libcint](https://github.com/sunqm/libcint) is a C library for GTO
//! (gaussian-type orbital) electronic integral, can be applied in
//! computational chemistry, and has already been applied extensively in
//! [PySCF](https://github.com/pyscf/pyscf).
//!
//! | Resources | Badges |
//! |--|--|
//! | Crate | [![Crate](https://img.shields.io/crates/v/libcint.svg)](https://crates.io/crates/libcint) |
//! | API Document | [![API Documentation](https://docs.rs/libcint/badge.svg)](https://docs.rs/libcint) |
//! | FFI Binding (libcint) | [v6.2.0](https://github.com/sunqm/libcint/tree/v6.2.0) |
//! | FFI Binding (qcint) | [v6.2.0](https://github.com/sunqm/qcint/tree/v6.2.0) |
//! | ECP Source Code (PySCF) | [v2.14.0](https://github.com/pyscf/pyscf/tree/v2.14.0) |
//!
//! This crate is not official bindgen project, neither
//! [libcint](https://github.com/sunqm/libcint),
//! [PySCF](https://github.com/pyscf/pyscf),
//! nor [REST](https://gitee.com/RESTGroup/rest).
//! It is originally intended to be some pioneer work for possible future
//! development of [rest_libcint](https://gitee.com/RESTGroup/rest_libcint) wrapper.
//!
//! # Minimal Example
//!
//! The most important function is [`CInt::integrate`]. It is somehow similar
//! to PySCF's `mol.intor(intor, aosym, shls_slice)`.
//!
//! ```rust
//! use libcint::prelude::*;
//!
//! // This is for testing and development purpose only.
//! // For actual usage, you should initialize `CInt` with your own molecule data.
//! let cint_data: CInt = init_h2o_def2_tzvp();
//!
//! // int1e_ipkin: the integral of kinetic energy operator with derivative on the first orbital
//!
//! // [mu, nu, comp], column-major, same data memory layout to PySCF's 2/3-center intor
//! let (out, shape): (Vec<f64>, Vec<usize>) = cint_data.integrate("int1e_ipkin", None, None).into();
//! assert_eq!(shape, [43, 43, 3]);
//!
//! // [comp, mu, nu], row-major, same shape to PySCF intor
//! let (out, shape): (Vec<f64>, Vec<usize>) = cint_data.integrate_row_major("int1e_ipkin", None, None).into();
//! assert_eq!(shape, [3, 43, 43]);
//! ```
//!
//! # Building Molecule from TOML/JSON (feature `bse`)
//!
//! With the `bse` feature, you can build a [`CIntMol`] object from TOML/JSON:
//!
//! ```ignore
//! use libcint::prelude::*;
//!
//! let toml = r#"
//! atom = "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987"
//! basis = "STO-3G"
//! "#;
//!
//! let mol = CIntMol::from_toml(toml).unwrap();
//! let ovlp = mol.cint.integrate("int1e_ovlp", None, None).out.unwrap();
//! ```
//!
//! See [`parse::mole::CIntMol`] and [`parse::serde_build`] for details.
//! Full TOML specification is documented in
//! `libcint/src/parse/from_toml_docs.md`.
//!
//! # For users from PySCF
//!
//! This library corresponds to some parts of `pyscf.gto` module in PySCF.
//!
//! Similar parts are:
//! - [`CInt::integrate`] corresponds to `mol.intor(intor, aosym, shls_slice)`;
//! - Various get/set/with-clauses methods;
//! - [`CInt`] corresponds to `(mol._atm, mol._bas, mol._env)` in PySCF (adding
//!   `mol._ecpbas` for ECP).
//! - [`CIntMol`] (feature `bse`) corresponds to PySCF's `Mole` class for
//!   molecule input parsing (atom coordinates, basis sets, ECP).
//!
//! Differences are:
//! - [`CInt`] only handles C-FFI data (`_atm`, `_bas`, `_env`, `_ecpbas`),
//!   similar to PySCF's internal arrays.
//! - [`CIntMol`] provides molecule input parsing (atom coordinates, basis/ECP
//!   from BSE, unit conversion, ghost atoms), comparable to PySCF's
//!   `Mole.build()` flow.
//!
//! # Installation and Cargo Features
//!
//! ## Install with pre-compiled `libcint.so` (recommended)
//!
//! If you have already compiled `libcint.so`, then put path of this shared
//! library in `CINT_DIR` or `LD_LIBRARY_PATH` (or `REST_EXT_DIR`). Then you
//! just use this library in your `Cargo.toml` file by
//!
//! ```toml
//! [dependencies]
//! libcint = { version = "0.5" }
//! ```
//!
//! ## Install and also build-from-source
//!
//! If you have not compiled `libcint.so` or `libcint.a`, then you are suggested
//! to use this library by specifying some cargo features:
//!
//! ```toml
//! [dependencies]
//! libcint = { version = "0.5", features = ["build_from_source", "static"] }
//! ```
//!
//! The source code will be automatically downloaded from github, and cargo will
//! handle the building process.
//!
//! If github is not accessible (offline builds, vendored sources), `CINT_SRC`
//! also accepts a local source directory, a source archive file, or the URL of
//! a source archive of [sunqm/libcint](https://github.com/sunqm/libcint) or
//! [sunqm/qcint](https://github.com/sunqm/qcint); pin an archive with
//! `CINT_SRC_SHA256`. A git URL or local repository still works as a source
//! mirror (`CINT_VER` selects the tag). Relative paths are resolved against the
//! `libcint-src` package directory, and a local git repository is cloned at
//! `CINT_VER`, so use a non-git directory to build modified sources.
//!
//! ## Cargo features
//!
//! - Default features: `bse` (molecule builder from TOML/JSON with Basis Set
//!   Exchange support).
//! - `bse`: Enable [`CIntMol::from_toml`] and [`CIntMol::from_json`] for
//!   building molecules from text input. Requires `bse` crate for basis set
//!   data.
//! - `api-v6_2`: Enable integrals added in libcint v6.2.0 (`int1e_ipr_origj`,
//!   `int1e_r_origj_ip`, `int1e_iprr_origj`, `int1e_rr_origj_ip`,
//!   `int3c1e_ipip1`, `int3c1e_ip1ip2`, `int3c1e_ipvip1`, `int3c1e_ipip2`).
//!   This feature is **not** enabled by default, so the default API stays
//!   compatible with libcint/qcint v6.1.2; enable it (e.g. `features =
//!   ["api-v6_2"]`) when using v6.2.0 and needing these integrals. With
//!   `qcint`, it additionally enables the integrals that qcint gained in
//!   v6.2.0. For a v6.1.2 source build, set `CINT_VER=v6.1.2`.
//! - `build_from_source`: Trigger of C language library libcint building. This
//!   performs by CMake; source code will be automatically downloaded from
//!   github (if environment variable `CINT_SRC` not specified).
//! - `static`: Use static library for linking. This will require static link
//!   `libcint.a`, and dynamic link `libquadmath.so` on targets other than macOS
//!   (unless `no_quadmath`).
//! - `qcint`: Use [sunqm/qcint](https://github.com/sunqm/qcint) instead of [sunqm/libcint](https://github.com/sunqm/libcint).
//!   Some integrals will not be available if `qcint` does not supports that.
//!   This will also change URL source if cargo feature `build_from_source`
//!   specified.
//! - `with_f12`: Whether F12 integrals (`int2e_stg`, `int2e_yp`, etc.) are
//!   supported.
//! - `with_4c1e`: Whether 4c1e integrals (`int4c1e`, etc.) are supported.
//! - `no_quadmath`: Do not link `libquadmath` with `static`, and with
//!   `build_from_source` build libcint without it (libcint then uses `long
//!   double` where it would use `__float128`). Needed for targets without
//!   libquadmath, such as aarch64 Linux.
//! - `keep_going`: With `build_from_source`, enable libcint's `KEEP_GOING`
//!   option so a Rys quadrature failure reports an error instead of exiting the
//!   process. Callers must check the returned integrals for numerical failures.
//!   Disabled by default; has no effect on prebuilt libraries.
//!
//! ## Shell environment variables
//!
//! - `CINT_DIR`, `LD_LIBRARY_PATH`, `REST_EXT_DIR`: Your compiled library path
//!   of `libcint.so` and `libcint.a`. This crate will try to find if this
//!   library is in directory root, or `directory_root/lib`. Not searched with
//!   cargo feature `build_from_source`, which links the library it built.
//! - `CINT_SRC`: Source of libcint or qcint: a git repository (URL or local
//!   path), a source archive (local path, or URL ending in `.tar.gz`, `.tgz`,
//!   `.tar.bz2`, `.tar.xz` or `.zip`), or an unpacked source directory. Only
//!   works with cargo feature `build_from_source`.
//! - `CINT_VER`: Version of libcint or qcint (default `v6.2.0`; must start with
//!   `v` prefix), as a git tag. Only works with cargo feature
//!   `build_from_source` and a git `CINT_SRC`.
//! - `CINT_SRC_SHA256`: Expected SHA256 of the `CINT_SRC` archive; the build
//!   fails on a mismatch. Only works with cargo feature `build_from_source` and
//!   an archive `CINT_SRC`; a git source or an unpacked directory with this
//!   variable set is rejected.
//!
//! # 50 lines RHF with Rust
//!
//! This is an example code to compute the RHF energy of H2O molecule using
//! - `libcint` as electronic integral library (corresponding to some parts of
//!   `pyscf.gto`);
//! - `rstsr` as tensor library (corresponding to NumPy and SciPy).
//!
//! You will see that except for import and preparation, the **code with core
//! algorithms is 27 lines** (blank and printing lines included). For
//! comparison, the **same code in Python is 23 lines**.
//!
//! ```ignore
#![doc = include_str!("../examples/h2o_rhf.rs")]
//! ```

pub mod prelude;

pub mod cint;

pub mod cint_initialize;

pub mod cint_prop;

pub mod cint_change;

pub mod cint_crafter;
pub mod cint_fill_col_major;
pub mod cint_fill_grids;
pub mod cint_fill_row_major;
pub mod cint_result;

pub mod ffi;

pub mod test_mol;
pub mod util;

pub mod gto;

#[cfg(feature = "bse")]
pub mod parse;

#[allow(unused_imports)]
use prelude::*;
