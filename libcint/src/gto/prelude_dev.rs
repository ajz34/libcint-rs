pub use crate::cint_ffi::*;
pub use crate::prelude::*;
pub use core::ffi::c_int;
pub use core::mem::transmute;

pub const BLKSIZE: usize = 56;
pub const BLKSIMD: usize = 7;
// maximum number of cartesian functions for a shell
// 128s, 42p, 21d, 12f, 8g, 6h, 4i, 3j
pub const NCTR_CART: usize = 128;
// maximum components to be evaluated
// derivatives: 1, 4, 10, 20, 35 (to 4-th order)
// ipipsp: 4 * 9
pub const NCOMP_MAX: usize = 36;

// 2 slots of int param[]
pub const POS_E1: usize = 0;
pub const TENSOR: usize = 1;

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

use core::ops::*;

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

#[repr(align(64))]
#[derive(Clone, Debug, Copy)]
pub struct Blk<T: Copy>(pub [T; BLKSIZE]);

#[allow(non_camel_case_types)]
pub type f64blk = Blk<f64>;

impl<T: Copy> Index<usize> for Blk<T> {
    type Output = T;
    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T: Copy> IndexMut<usize> for Blk<T> {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<T: Copy> Blk<T> {
    /// # Safety
    ///
    /// This function returns an uninitialized `BlkF64`. The caller must ensure
    /// that the returned value is properly initialized before use.
    #[inline(always)]
    #[allow(clippy::uninit_assumed_init)]
    #[allow(invalid_value)]
    pub const unsafe fn uninit() -> Self {
        Blk(core::mem::MaybeUninit::uninit().assume_init())
    }

    pub const fn splat(val: T) -> Self {
        Blk([val; BLKSIZE])
    }

    pub const fn fill(&mut self, val: T) {
        let mut i = 0;
        while i < BLKSIZE {
            self.0[i] = val;
            i += 1;
        }
    }

    /// # Safety
    ///
    /// The caller must ensure that `src` has at least `BLKSIZE` elements.
    #[inline(always)]
    pub const unsafe fn read_ensure(&mut self, src: &[T]) {
        let mut i = 0;
        while i < BLKSIZE {
            self.0[i] = src[i];
            i += 1;
        }
    }

    #[inline]
    pub fn read(&mut self, src: &[T]) {
        let len_slc = if src.len() < BLKSIZE { src.len() } else { BLKSIZE };
        if len_slc >= BLKSIZE {
            unsafe { self.read_ensure(src) }
        } else {
            self.0.copy_from_slice(src);
        }
    }

    /// # Safety
    ///
    /// The caller must ensure that `dst` has at least `BLKSIZE` elements.
    #[inline(always)]
    pub const unsafe fn write_ensure(&self, dst: &mut [T]) {
        let mut i = 0;
        while i < BLKSIZE {
            dst[i] = self.0[i];
            i += 1;
        }
    }

    #[inline]
    pub fn write(&self, dst: &mut [T]) {
        let len_slc = if dst.len() < BLKSIZE { dst.len() } else { BLKSIZE };
        if len_slc >= BLKSIZE {
            unsafe { self.write_ensure(dst) }
        } else {
            dst.copy_from_slice(&self.0[..len_slc]);
        }
    }
}

#[test]
fn test_blkf64_uninit() {
    #[allow(clippy::uninit_assumed_init)]
    #[allow(invalid_value)]
    let v: [[f64blk; 10]; 10] = unsafe { [[f64blk::uninit(); 10]; 10] };
    println!("{:?}", v);
}
