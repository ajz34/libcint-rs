pub use crate::gto::deriv_impl::deriv_0::*;
pub use crate::gto::deriv_impl::deriv_1::*;
pub use crate::gto::deriv_impl::deriv_2::*;
pub use crate::gto::deriv_impl::deriv_3::*;
pub use crate::gto::deriv_impl::deriv_4::*;
pub use crate::gto::deriv_impl::deriv_ig::*;
pub use crate::gto::deriv_impl::deriv_ip::*;
pub use crate::gto::deriv_impl::deriv_ipig::*;
pub use crate::gto::deriv_impl::deriv_ipipsp::*;
pub use crate::gto::deriv_impl::deriv_ipr::*;
pub use crate::gto::deriv_impl::deriv_iprc::*;
pub use crate::gto::deriv_impl::deriv_ipsp::*;
pub use crate::gto::deriv_impl::deriv_sp::*;

pub use crate::gto::deriv_util::*;
pub use crate::gto::grid_ao_drv::*;

pub(crate) use crate::cint_ffi::*;
pub(crate) use crate::prelude::*;
pub(crate) use core::ffi::c_int;
pub(crate) use core::mem::transmute;
pub(crate) use core::mem::MaybeUninit;

/// Number of SIMD lanes for double precision (SIMD for double).
///
/// It should be 8 for AVX-512.
///
/// For other architectures, it is better to be 4 (AVX2) or 2 (SSE2), but
/// currently we fix it to 8. So other architectures will still be optimized by
/// SIMD, it can be somehow sub-optimal due to insufficient instructions in
/// registers.
pub const SIMDD: usize = 8;

/// Block size for grid AO evaluation.
///
/// This [`BLKSIZE`] value (48) should be multiple of [`SIMDD`] (8). It should
/// be better to be multiple of 16 because of microkernel in matmul usually
/// requires 2 lanes of SIMD (for AVX-512, it is 16 f64). Currently this value
/// is fixed to 48, in that for the most-used deriv1 case, the processed grids
/// per function call is $(n_\mathrm{comp}, n_\mathrm{ctr} \times
/// n_\mathrm{cart}(l), n_\mathrm{grids})$ will usually be up to (4, 15, 48) or
/// 22.5 KB, which fits in L1d cache (32 KB in most micro-architectures)
/// together with other data.
pub const BLKSIZE: usize = 48;

/// Number of SIMD blocks in a block.
///
/// This value will be 7 for the current settings, and will be used in loops
/// over blocks.
pub const BLKSIMDD: usize = BLKSIZE / SIMDD;

/// Number of histogram bins for GTO screening [`gto_screen_index`].
pub const NBINS: u8 = 100;

/// Cutoff value for GTO screening [`gto_screen_index`].
pub const CUTOFF: f64 = 1e-22;

/// Cutoff value to consider a GTO exponents as zero in implementations of trait
/// function [`GtoEvalAPI::gto_shell_eval`].
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
/// libcint C function: `CINTcommon_fac_sp`. This function is not included in
/// the public C API of libcint, so we re-implement it in Rust here.
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
