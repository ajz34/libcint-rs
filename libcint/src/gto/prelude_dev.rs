pub use crate::cint_ffi::*;
pub use crate::prelude::*;
pub use core::ffi::c_int;
pub use core::mem::transmute;

pub const BLKSIZE: usize = 56;
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

unsafe extern "C" {
    pub fn CINTcommon_fac_sp(l: c_int) -> f64;
}

#[repr(align(64))]
#[derive(Clone, Debug, Copy)]
pub struct BlkF64(pub [f64; BLKSIZE]);

impl std::ops::Index<usize> for BlkF64 {
    type Output = f64;
    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl std::ops::IndexMut<usize> for BlkF64 {
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
