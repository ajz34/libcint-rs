//! Use `libcint::prelude::*` to import this crate.

pub use crate::cint::{CInt, CIntKind, CIntOptimizer, CIntOutput, CIntSymm, CIntType};
pub use crate::cint_change::fakemol_for_charges;
pub use crate::cint_result::{CIntError, CIntErrorBase};
pub use crate::ffi::wrapper_traits::Integrator;
pub use crate::util::ShlsSlice;

#[cfg(feature = "bse")]
pub use crate::parse::mole::{CIntMol, CIntMolInput, CIntMolInputBuilder};

// for doc testing
pub use crate::test_mol::{cint_fingerprint, cint_fp, init_c10h22_def2_qzvp, init_h2o_def2_jk, init_h2o_def2_tzvp, init_sb2me4_cc_pvtz};

/// The minimum number for a chunk to be processed in parallel by Rayon.
pub const RAYON_PAR_MIN: usize = 32;

// for developing

#[allow(unused_imports)]
pub mod prelude_dev {
    pub use crate::cint::*;
    pub use crate::cint_result::CIntResultAPI;
    pub use crate::ffi::cecp_ffi;
    pub use crate::ffi::cecp_ffi::ECPOpt;
    pub use crate::ffi::cint_ffi;
    pub use crate::ffi::cint_ffi::CINTOpt;
    pub use crate::gto::prelude_dev::*;
    pub use crate::util::*;
    pub use crate::{cint_error, cint_raise, cint_trace};
    pub use derive_builder::{Builder, UninitializedFieldError};
    pub use duplicate::duplicate_item;
    pub use indexmap::IndexMap;
    pub use itertools::Itertools;
    pub use num::complex::ComplexFloat;
    pub use num::{Complex, ToPrimitive};
    pub use rayon::prelude::*;
    pub use rstsr_common::layout::exports::IndexedIterLayout;
    pub use rstsr_common::prelude::*;
    pub use serde::{Deserialize, Deserializer, Serialize, Serializer};
    pub use std::any::Any;
    pub use std::backtrace::Backtrace;
    pub use std::collections::{BTreeMap, HashMap};
    pub use std::error::Error;
    pub use std::ffi::{c_int, c_void};
    pub use std::fmt::{Debug, Display, Write};
    pub use std::ptr::NonNull;
    pub use std::ptr::{null, null_mut};
    pub use CIntType::{Cartesian, Spheric, Spinor};
}

pub(crate) use prelude_dev::*;

#[cfg(feature = "bse")]
pub(crate) use bse::prelude::*;
