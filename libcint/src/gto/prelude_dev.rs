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

/* #region simple f64x8 */

// TODO: use fearless_simd or pulp, they are more heavy but more complete SIMD
// abstraction crates

#[repr(align(64))]
#[derive(Clone, Debug, Copy)]
#[allow(non_camel_case_types)]
pub struct f64x8(pub [f64; 8]);

impl Index<usize> for f64x8 {
    type Output = f64;
    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl IndexMut<usize> for f64x8 {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl f64x8 {
    #[inline(always)]
    pub const fn zero() -> Self {
        f64x8([0.0f64; 8])
    }

    /// # Safety
    ///
    /// This function returns an uninitialized object. The caller must ensure
    /// that the returned value is properly initialized before use.
    #[inline(always)]
    #[allow(clippy::uninit_assumed_init)]
    #[allow(invalid_value)]
    pub const unsafe fn uninit() -> Self {
        core::mem::MaybeUninit::uninit().assume_init()
    }

    #[inline(always)]
    pub const fn splat(val: f64) -> Self {
        f64x8([val; 8])
    }

    #[inline(always)]
    pub const fn fill(&mut self, val: f64) {
        self.0 = [val; 8];
    }

    pub fn map<F>(&self, f: F) -> Self
    where
        F: Fn(f64) -> f64,
    {
        f64x8(self.0.map(f))
    }
}

#[duplicate_item(
    Trait trait_fn;
    [Add] [add];
    [Sub] [sub];
    [Mul] [mul];
    [Div] [div];
)]
mod impl_trait_for_f64x8 {
    use super::*;

    // simd * simd
    impl Trait for f64x8 {
        type Output = Self;
        #[inline(always)]
        fn trait_fn(self, rhs: Self) -> Self::Output {
            f64x8([
                f64::trait_fn(self[0], rhs[0]),
                f64::trait_fn(self[1], rhs[1]),
                f64::trait_fn(self[2], rhs[2]),
                f64::trait_fn(self[3], rhs[3]),
                f64::trait_fn(self[4], rhs[4]),
                f64::trait_fn(self[5], rhs[5]),
                f64::trait_fn(self[6], rhs[6]),
                f64::trait_fn(self[7], rhs[7]),
            ])
        }
    }

    // simd * scalar
    impl Trait<f64> for f64x8 {
        type Output = Self;
        #[inline(always)]
        fn trait_fn(self, rhs: f64) -> Self::Output {
            f64x8([
                f64::trait_fn(self[0], rhs),
                f64::trait_fn(self[1], rhs),
                f64::trait_fn(self[2], rhs),
                f64::trait_fn(self[3], rhs),
                f64::trait_fn(self[4], rhs),
                f64::trait_fn(self[5], rhs),
                f64::trait_fn(self[6], rhs),
                f64::trait_fn(self[7], rhs),
            ])
        }
    }

    // scalar * simd
    impl Trait<f64x8> for f64 {
        type Output = f64x8;
        #[inline(always)]
        fn trait_fn(self, rhs: f64x8) -> Self::Output {
            f64x8([
                f64::trait_fn(self, rhs[0]),
                f64::trait_fn(self, rhs[1]),
                f64::trait_fn(self, rhs[2]),
                f64::trait_fn(self, rhs[3]),
                f64::trait_fn(self, rhs[4]),
                f64::trait_fn(self, rhs[5]),
                f64::trait_fn(self, rhs[6]),
                f64::trait_fn(self, rhs[7]),
            ])
        }
    }
}

#[duplicate_item(
    Trait trait_fn;
    [AddAssign] [add_assign];
    [SubAssign] [sub_assign];
    [MulAssign] [mul_assign];
    [DivAssign] [div_assign];
)]
impl Trait for f64x8 {
    #[inline(always)]
    fn trait_fn(&mut self, rhs: Self) {
        self[0].trait_fn(rhs[0]);
        self[1].trait_fn(rhs[1]);
        self[2].trait_fn(rhs[2]);
        self[3].trait_fn(rhs[3]);
        self[4].trait_fn(rhs[4]);
        self[5].trait_fn(rhs[5]);
        self[6].trait_fn(rhs[6]);
        self[7].trait_fn(rhs[7]);
    }
}

impl f64x8 {
    #[inline]
    pub fn add_mul(self, b: f64x8, c: f64x8) -> f64x8 {
        f64x8([
            self[0].mul_add(b[0], c[0]),
            self[1].mul_add(b[1], c[1]),
            self[2].mul_add(b[2], c[2]),
            self[3].mul_add(b[3], c[3]),
            self[4].mul_add(b[4], c[4]),
            self[5].mul_add(b[5], c[5]),
            self[6].mul_add(b[6], c[6]),
            self[7].mul_add(b[7], c[7]),
        ])
    }
}

/* #endregion */

#[repr(align(64))]
#[derive(Clone, Debug, Copy)]
#[allow(non_camel_case_types)]
pub struct f64blk(pub [f64; BLKSIZE]);

impl Index<usize> for f64blk {
    type Output = f64;
    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl IndexMut<usize> for f64blk {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl f64blk {
    /// # Safety
    ///
    /// This function returns an uninitialized `BlkF64`. The caller must ensure
    /// that the returned value is properly initialized before use.
    #[inline(always)]
    #[allow(clippy::uninit_assumed_init)]
    #[allow(invalid_value)]
    pub const unsafe fn uninit() -> Self {
        f64blk(core::mem::MaybeUninit::uninit().assume_init())
    }

    pub const fn zero() -> Self {
        f64blk([0.0f64; BLKSIZE])
    }

    pub const fn splat(val: f64) -> Self {
        f64blk([val; BLKSIZE])
    }

    pub const fn fill(&mut self, val: f64) {
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
    pub const unsafe fn read_ensure(&mut self, src: &[f64]) {
        let mut i = 0;
        while i < BLKSIZE {
            self.0[i] = src[i];
            i += 1;
        }
    }

    #[inline]
    pub fn read(&mut self, src: &[f64]) {
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
    pub const unsafe fn write_ensure(&self, dst: &mut [f64]) {
        let mut i = 0;
        while i < BLKSIZE {
            dst[i] = self.0[i];
            i += 1;
        }
    }

    #[inline]
    pub fn write(&self, dst: &mut [f64]) {
        let len_slc = if dst.len() < BLKSIZE { dst.len() } else { BLKSIZE };
        if len_slc >= BLKSIZE {
            unsafe { self.write_ensure(dst) }
        } else {
            dst.copy_from_slice(&self.0[..len_slc]);
        }
    }

    pub const fn as_f64x8_slice(&self) -> &[f64x8; BLKSIMD] {
        unsafe { transmute(&self.0) }
    }

    pub const fn as_f64x8_slice_mut(&mut self) -> &mut [f64x8; BLKSIMD] {
        unsafe { transmute(&mut self.0) }
    }
}

#[test]
fn test_blkf64_uninit() {
    #[allow(clippy::uninit_assumed_init)]
    #[allow(invalid_value)]
    let v: [[f64blk; 10]; 10] = unsafe { [[f64blk::uninit(); 10]; 10] };
    println!("{:?}", v);
}
