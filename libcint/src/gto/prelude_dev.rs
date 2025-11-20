pub use crate::cint_ffi::*;
pub use crate::prelude::*;
pub use core::ffi::c_int;
pub use core::mem::transmute;

pub const BLKSIZE: usize = 56;
pub const BLKSIMD: usize = 7;
// 128s42p21d12f8g6h4i3j
pub const NCTR_CART: usize = 128;

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

use core::ops::*;

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
                f64::trait_fn(self.0[0], rhs.0[0]),
                f64::trait_fn(self.0[1], rhs.0[1]),
                f64::trait_fn(self.0[2], rhs.0[2]),
                f64::trait_fn(self.0[3], rhs.0[3]),
                f64::trait_fn(self.0[4], rhs.0[4]),
                f64::trait_fn(self.0[5], rhs.0[5]),
                f64::trait_fn(self.0[6], rhs.0[6]),
                f64::trait_fn(self.0[7], rhs.0[7]),
            ])
        }
    }

    // simd * scalar
    impl Trait<f64> for f64x8 {
        type Output = Self;
        #[inline(always)]
        fn trait_fn(self, rhs: f64) -> Self::Output {
            f64x8([
                f64::trait_fn(self.0[0], rhs),
                f64::trait_fn(self.0[1], rhs),
                f64::trait_fn(self.0[2], rhs),
                f64::trait_fn(self.0[3], rhs),
                f64::trait_fn(self.0[4], rhs),
                f64::trait_fn(self.0[5], rhs),
                f64::trait_fn(self.0[6], rhs),
                f64::trait_fn(self.0[7], rhs),
            ])
        }
    }

    // scalar * simd
    impl Trait<f64x8> for f64 {
        type Output = f64x8;
        #[inline(always)]
        fn trait_fn(self, rhs: f64x8) -> Self::Output {
            f64x8([
                f64::trait_fn(self, rhs.0[0]),
                f64::trait_fn(self, rhs.0[1]),
                f64::trait_fn(self, rhs.0[2]),
                f64::trait_fn(self, rhs.0[3]),
                f64::trait_fn(self, rhs.0[4]),
                f64::trait_fn(self, rhs.0[5]),
                f64::trait_fn(self, rhs.0[6]),
                f64::trait_fn(self, rhs.0[7]),
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
        self.0[0].trait_fn(rhs.0[0]);
        self.0[1].trait_fn(rhs.0[1]);
        self.0[2].trait_fn(rhs.0[2]);
        self.0[3].trait_fn(rhs.0[3]);
        self.0[4].trait_fn(rhs.0[4]);
        self.0[5].trait_fn(rhs.0[5]);
        self.0[6].trait_fn(rhs.0[6]);
        self.0[7].trait_fn(rhs.0[7]);
    }
}

impl f64x8 {
    #[inline(always)]
    pub fn fma(&mut self, a: &f64x8, b: &f64x8) {
        self.0[0] = self.0[0].mul_add(a.0[0], b.0[0]);
        self.0[1] = self.0[1].mul_add(a.0[1], b.0[1]);
        self.0[2] = self.0[2].mul_add(a.0[2], b.0[2]);
        self.0[3] = self.0[3].mul_add(a.0[3], b.0[3]);
        self.0[4] = self.0[4].mul_add(a.0[4], b.0[4]);
        self.0[5] = self.0[5].mul_add(a.0[5], b.0[5]);
        self.0[6] = self.0[6].mul_add(a.0[6], b.0[6]);
        self.0[7] = self.0[7].mul_add(a.0[7], b.0[7]);
    }
}

/* #endregion */

#[repr(align(64))]
#[derive(Clone, Debug, Copy)]
pub struct BlkF64(pub [f64; BLKSIZE]);

impl Index<usize> for BlkF64 {
    type Output = f64;
    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl IndexMut<usize> for BlkF64 {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl BlkF64 {
    /// # Safety
    ///
    /// This function returns an uninitialized `BlkF64`. The caller must ensure
    /// that the returned value is properly initialized before use.
    #[inline(always)]
    #[allow(clippy::uninit_assumed_init)]
    #[allow(invalid_value)]
    pub const unsafe fn uninit() -> Self {
        BlkF64(core::mem::MaybeUninit::uninit().assume_init())
    }

    pub const fn zero() -> Self {
        BlkF64([0.0f64; BLKSIZE])
    }

    pub const fn splat(val: f64) -> Self {
        BlkF64([val; BLKSIZE])
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
}

#[test]
fn test_blkf64_uninit() {
    #[allow(clippy::uninit_assumed_init)]
    #[allow(invalid_value)]
    let v: [[BlkF64; 10]; 10] = unsafe { [[BlkF64::uninit(); 10]; 10] };
    println!("{:?}", v);
}
