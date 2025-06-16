#![allow(clippy::excessive_precision)]
#![allow(clippy::needless_range_loop)]

//! This project contains libcint (C language) FFI bindings, wrapper and
//! build-from-source.
//!
//! [libcint](https://github.com/sunqm/libcint) is a C library for GTO
//! (gaussian-type orbital) electronic integral, can be applied in
//! computational chemistry, and has already been applied extensively in
//! [PySCF](https://github.com/pyscf/pyscf).
//!
//! **This is an alpha version and will soon be promoted to v0.1.** After some
//! documentation update and probably symmetry enhancements, a version v0.1 will
//! be released.
//!
//! | Resources | Badges |
//! |--|--|
//! | Crate | [![Crate](https://img.shields.io/crates/v/libcint.svg)](https://crates.io/crates/libcint) |
//! | API Document | [![API Documentation](https://docs.rs/libcint/badge.svg)](https://docs.rs/libcint) |
//! | FFI Binding (libcint) | [v6.1.2](https://github.com/sunqm/libcint/tree/v6.1.2) |
//! | FFI Binding (qcint) | [v6.1.2](https://github.com/sunqm/qcint/tree/v6.1.2) |
//! | ECP Source Code (PySCF) | [v2.9.0](https://github.com/pyscf/pyscf/tree/v2.9.0) |
//!
//! This crate is not official bindgen project, nither
//! [libcint](https://github.com/sunqm/libcint),
//! [PySCF](https://github.com/pyscf/pyscf),
//! nor [REST](https://gitee.com/RESTGroup/rest).
//! It is originally intended to be some pioneer work for possible future
//! development of [rest_libcint](https://gitee.com/RESTGroup/rest_libcint) wrapper.
//!
//! # Minimal Example
//!
//! The most important function is be [`CInt::integrate`]. It is somehow similar
//! to PySCF's `mol.intor(intor, aosym, shls_slice)`.
//!
//! ```rust
//! use libcint::prelude::*;
//!
//! // This is for testing and development purpose only.
//! // For actual usage, you should initialize `CInt` with your own molecule data.
//! let cint_data: CInt = init_h2o_def2_tzvp();
//!
//! // Obtain the integral of kinetic energy operator with derivative on the first orbital:
//! let (out, shape) = cint_data.integrate("int1e_ipkin", None, None).into();
//! assert_eq!(shape, [43, 43, 3]);
//! ```
//!
//! # For most users
//!
//! You are probably most interested in the following parts:
//! - [`CInt::integrate`] and [`CInt::integrate_spinor`] for obtaining integrals
//!   (PySCF's `mol.intor` counterpart);
//! - [`CInt`] for virtually all methods for integrals and molecule data
//!   handling.
//!
//! For some advanced users, you may also want to use
//! - [`CInt::integrate_cross`] for integral with multiple `CInt` molecule data;
//! - `CInt::with_xxx` for temporarily changing molecule data (PySCF's
//!   `mol.with_xxx` with-clause);
//!
//! For even more advanced users, you may want to use
//! - [`CInt::integrate_args_builder`], [`CInt::integrate_args_builder_spinor`]
//!   for full arguments (currently you can also specify `out` as argument,
//!   meaning that you can use your own memory buffer to store the result);
//!   see also
//!   [derive_builder](https://docs.rs/derive_builder/latest/derive_builder/)
//!   for how to use the builder pattern.
//! - And then use [`CInt::integrate_with_args`] and
//!   [`CInt::integrate_with_args_spinor`] to integrate with the arguments.
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
//!
//! Differences are:
//! - Recall that PySCF's `Mole` class handles three parts: Python input from
//!   user (`mol.atom`, `mol.basis`, etc.), Python intermediates (`mol._basis`,
//!   etc.), data for C-FFI (`mol._atm`, `mol._bas`, `mol._env`). In rust's
//!   [`CInt`], it only handles the last part.
//! - This library does not handle molecule input (coords, spin, charge, etc.)
//!   and basis set parse. Molecule initialization should be performed by user.
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
//! libcint = { version = "0.0.1" }
//! ```
//!
//! ## Install and also build-from-source
//!
//! If you have not compiled `libcint.so` or `libcint.a`, then you are suggested
//! to use this library by specifying some cargo features:
//!
//! ```toml
//! [dependencies]
//! libcint = { version = "0.0.1", features = ["build_from_source", "static"] }
//! ```
//!
//! If access to github is not available, you can use environment variable
//! `CINT_SRC` to specify source mirror of
//! [sunqm/libcint](https://github.com/sunqm/libcint)
//! or [sunqm/qcint](https://github.com/sunqm/qcint).
//!
//! ## Cargo features
//!
//! - Default features: None of any listed below (use library provided by system
//!   or user, using [sunqm/libcint](https://github.com/sunqm/libcint), dynamic
//!   linking, without F12 and 4c1e support).
//! - `build_from_source`: Trigger of C language library libcint building. This
//!   performs by CMake, will source code download from github (if environment
//!   variable `CINT_SRC` not specified).
//! - `static`: Use static library for linking. This will require static link
//!   `libcint.a`, and dynamic link `libquadmath.so`.
//! - `qcint`: Use [sunqm/qcint](https://github.com/sunqm/qcint) instead of [sunqm/libcint](https://github.com/sunqm/libcint).
//!   Some integrals will not be available if `qcint` does not supports that.
//!   This will also change URL source if cargo feature `build_from_source`
//!   specified.
//! - `with_f12`: Whether F12 integrals (`int2e_stg`, `int2e_yp`, etc.) are
//!   supported.
//! - `with_4c1e`: Whether 4c1e integrals (`int4c1e`, etc.) are supported.
//!
//! ## Shell environment variables
//!
//! - `CINT_DIR`, `LD_LIBRARY_PATH`, `REST_EXT_DIR`: Your compiled library path
//!   of `libcint.so` and `libcint.a`. This crate will try to find if this
//!   library is in directory root, or `directory_root/lib`. May override the
//!   library built by cargo feature `build_from_source`.
//! - `CINT_SRC`: Source of libcint or qcint (must be a git repository). Only
//!   works with cargo feature `build_from_source`.
//! - `CINT_VIR`: Version of libcint or qcint (e.g. `v6.1.2`, must starts with
//!   `v` prefix). Only works with cargo feature `build_from_source`.
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
#![doc = include_str!("../assets/h2o_rhf.rs")]
//! ```

pub mod cint;
pub mod cint_change;
pub mod cint_crafter;
pub mod cint_prop;
pub mod ffi;
pub mod prelude;
pub mod test_mol;
pub mod util;

#[allow(unused_imports)]
use prelude::*;
