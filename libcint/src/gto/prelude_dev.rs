pub use crate::gto::deriv_0::*;
pub use crate::gto::deriv_1::*;
pub use crate::gto::grid_ao_drv::*;
use num::traits::MulAdd;
use num::traits::NumAssignOps;
use num::Num;

pub use crate::cint_ffi::*;
pub use crate::prelude::*;
pub use core::ffi::c_int;
pub use core::mem::transmute;
pub use core::mem::MaybeUninit;

pub const BLKSIZE: usize = 56;
pub const SIMDD: usize = 8;
pub const BLKSIMD: usize = BLKSIZE / SIMDD;
pub const NBINS: u8 = 100;
pub const CUTOFF: f64 = 1e-15;

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
pub struct FpSimd<T: Copy>(pub [T; SIMDD]);

impl<T: Copy> Index<usize> for FpSimd<T> {
    type Output = T;
    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T: Copy> IndexMut<usize> for FpSimd<T> {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<T: Num + Copy> FpSimd<T> {
    #[inline(always)]
    pub fn zero() -> Self {
        FpSimd([T::zero(); 8])
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
    pub const fn splat(val: T) -> Self {
        FpSimd([val; 8])
    }

    #[inline(always)]
    pub const fn fill(&mut self, val: T) {
        self.0 = [val; 8];
    }

    pub fn map<F>(&self, f: F) -> Self
    where
        F: Fn(T) -> T,
    {
        FpSimd(self.0.map(f))
    }
}

#[duplicate_item(
    Trait trait_fn;
    [Add] [add];
    [Sub] [sub];
    [Mul] [mul];
    [Div] [div];
)]
mod impl_trait_for_f64simd {
    use super::*;

    // simd * simd
    impl<T: Num + Copy> Trait for FpSimd<T> {
        type Output = Self;
        #[inline(always)]
        fn trait_fn(self, rhs: Self) -> Self::Output {
            FpSimd([
                T::trait_fn(self[0], rhs[0]),
                T::trait_fn(self[1], rhs[1]),
                T::trait_fn(self[2], rhs[2]),
                T::trait_fn(self[3], rhs[3]),
                T::trait_fn(self[4], rhs[4]),
                T::trait_fn(self[5], rhs[5]),
                T::trait_fn(self[6], rhs[6]),
                T::trait_fn(self[7], rhs[7]),
            ])
        }
    }

    // simd * scalar
    impl<T: Num + Copy> Trait<T> for FpSimd<T> {
        type Output = Self;
        #[inline(always)]
        fn trait_fn(self, rhs: T) -> Self::Output {
            FpSimd([
                T::trait_fn(self[0], rhs),
                T::trait_fn(self[1], rhs),
                T::trait_fn(self[2], rhs),
                T::trait_fn(self[3], rhs),
                T::trait_fn(self[4], rhs),
                T::trait_fn(self[5], rhs),
                T::trait_fn(self[6], rhs),
                T::trait_fn(self[7], rhs),
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
impl<T: NumAssignOps + Copy> Trait for FpSimd<T> {
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

impl<T> FpSimd<T>
where
    T: MulAdd<Output = T> + Copy,
{
    #[inline(always)]
    pub fn add_mul(self, b: FpSimd<T>, c: FpSimd<T>) -> FpSimd<T> {
        FpSimd([
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

/* #region Blk<T> */

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
        Blk(MaybeUninit::uninit().assume_init())
    }

    #[inline(always)]
    pub const fn splat(val: T) -> Self {
        Blk([val; BLKSIZE])
    }

    #[inline(always)]
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

impl<T: Copy> Blk<T> {
    #[inline(always)]
    pub const fn as_simd_slice(&self) -> &[FpSimd<T>; BLKSIMD] {
        unsafe { transmute(&self.0) }
    }

    #[inline(always)]
    pub const fn as_simd_slice_mut(&mut self) -> &mut [FpSimd<T>; BLKSIMD] {
        unsafe { transmute(&mut self.0) }
    }

    #[inline(always)]
    pub const fn get_simd(&self, index: usize) -> FpSimd<T> {
        let slice = self.as_simd_slice();
        slice[index]
    }

    #[inline(always)]
    pub const fn get_simd_mut(&mut self, index: usize) -> &mut FpSimd<T> {
        let slice = self.as_simd_slice_mut();
        &mut slice[index]
    }
}

/* #endregion */

#[test]
fn test_blkf64_uninit() {
    #[allow(clippy::uninit_assumed_init)]
    #[allow(invalid_value)]
    let v: [[f64blk; 10]; 10] = unsafe { [[f64blk::uninit(); 10]; 10] };
    println!("{:?}", v);
}
