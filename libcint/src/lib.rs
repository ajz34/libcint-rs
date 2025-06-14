#![allow(clippy::excessive_precision)]
#![allow(clippy::needless_range_loop)]

//! FFI bindings and wrapper for libcint (C language) electronic integral
//! library for rust.
//!
//! # Minimal Example
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
//! # To use this library
//!
//! Firstly, to use this library, you need to link the C library `cint` in your
//! `build.rs`:
//!
//! ```ignore
//! // build.rs
//! println!("cargo:rustc-link-search=native={<your_libcint.so_path>}");
//! println!("cargo:rustc-link-lib=cint");
//! ```
//!
//! Or either use the `libcint-src` crate, which will also link the libcint
//! library with some proper configurations.
//!
//! We refer to github actions for hints of how to setup crate `libcint` with
//! `libcint-src`.
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
