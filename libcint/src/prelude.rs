//! Use `libcint::prelude::*` to import this crate.

pub use crate::cint::{CInt, CIntError, CIntKind, CIntOptimizer, CIntOutput, CIntSymm, CIntType};
pub use crate::cint_change::fakemol_for_charges;
pub use crate::ffi::wrapper_traits::Integrator;
pub use crate::util::ShlsSlice;

// for doc testing
pub use crate::test_mol::{cint_fingerprint, cint_fp, init_c10h22_def2_qzvp, init_h2o_def2_jk, init_h2o_def2_tzvp, init_sb2me4_cc_pvtz};

// for developing
pub(crate) use crate::cint::*;
pub(crate) use crate::ffi::cecp_ffi;
pub(crate) use crate::ffi::cecp_ffi::ECPOpt;
pub(crate) use crate::ffi::cint_ffi;
pub(crate) use crate::ffi::cint_ffi::CINTOpt;
pub(crate) use crate::util::*;
pub(crate) use crate::{cint_raise, cint_trace};
pub(crate) use derive_builder::{Builder, UninitializedFieldError};
pub(crate) use duplicate::duplicate_item;
pub(crate) use itertools::Itertools;
pub(crate) use num::complex::{Complex, ComplexFloat};
pub(crate) use std::any::Any;
pub(crate) use std::error::Error;
pub(crate) use std::ffi::{c_int, c_void};
pub(crate) use std::fmt::{Debug, Display, Write};
pub(crate) use std::ptr::NonNull;
pub(crate) use std::ptr::{null, null_mut};
pub(crate) use CIntType::{Cartesian, Spheric, Spinor};
