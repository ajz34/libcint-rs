//! GTO (gaussian-type atomic orbitals) values on coordinate grids.
#![doc = include_str!("gto_notation.md")]

#[allow(unused_imports)]
use crate::gto::prelude_dev::*;

pub mod deriv_impl;
pub mod deriv_util;
pub mod grid_ao_drv;
pub mod gto_crafter;
pub mod prelude_dev;
