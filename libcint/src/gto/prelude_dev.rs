pub use crate::gto::deriv_0::*;
pub use crate::gto::deriv_1::*;
pub use crate::gto::deriv_2::*;
pub use crate::gto::deriv_3::*;
pub use crate::gto::deriv_4::*;
pub use crate::gto::deriv_ig::*;
pub use crate::gto::deriv_ip::*;
pub use crate::gto::deriv_ipig::*;
pub use crate::gto::deriv_ipr::*;
pub use crate::gto::deriv_iprc::*;

pub use crate::gto::deriv_util::*;
pub use crate::gto::grid_ao_drv::*;

pub use crate::cint_ffi::*;
pub use crate::prelude::*;
pub use core::ffi::c_int;
pub use core::mem::transmute;
pub use core::mem::MaybeUninit;

pub const BLKSIZE: usize = 56;
pub const SIMDD: usize = 8;
pub const BLKSIMDD: usize = BLKSIZE / SIMDD;
pub const NBINS: u8 = 100;
pub const CUTOFF: f64 = 1e-15;
pub const GTOZERO: f64 = 1e-30;

pub const X: usize = 0;
pub const Y: usize = 1;
pub const Z: usize = 2;

pub const XX: usize = 0;
pub const XY: usize = 1;
pub const XZ: usize = 2;
pub const YY: usize = 3;
pub const YZ: usize = 4;
pub const ZZ: usize = 5;

pub const XXX: usize = 0;
pub const XXY: usize = 1;
pub const XXZ: usize = 2;
pub const XYY: usize = 3;
pub const XYZ: usize = 4;
pub const XZZ: usize = 5;
pub const YYY: usize = 6;
pub const YYZ: usize = 7;
pub const YZZ: usize = 8;
pub const ZZZ: usize = 9;

/// Common factor for spherical harmonics, by libcint's convention.
///
/// For efficiency, libcint will pre-scale normalization factor to the integral,
/// for s and p functions. For other angular momenta, the scaling will be
/// performed along with cartesian-spherical transformation.
///
/// - $l = 0$ (s): $\sqrt{\frac{1}{4\pi}}$
/// - $l = 1$ (p): $\sqrt{\frac{3}{4\pi}}$
/// - otherwise: 1
///
/// libcint C function: `CINTcommon_fac_sp`.
pub fn cint_common_fac_sp(l: c_int) -> f64 {
    match l {
        0 => 0.282094791773878143, // sqrt(1 / (4 * pi))
        1 => 0.488602511902919921, // sqrt(3 / (4 * pi))
        _ => 1.0,
    }
}

unsafe extern "C" {
    pub fn CINTcommon_fac_sp(l: c_int) -> f64;
}
