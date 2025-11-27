use crate::gto::prelude_dev::*;
use num::traits::{MulAdd, NumAssignOps};
use num::Num;

/* #region simple f64x8 */

// TODO: use fearless_simd or pulp, they are more heavy but more complete SIMD
// abstraction crates

#[repr(align(64))]
#[derive(Clone, Debug, Copy)]
pub struct FpSimd<T: Copy, const N: usize = SIMDD>(pub [T; N]);

#[allow(non_camel_case_types)]
pub type f64simd = FpSimd<f64, SIMDD>;

impl<T: Copy, const N: usize> Index<usize> for FpSimd<T, N> {
    type Output = T;
    #[inline(always)]
    fn index(&self, index: usize) -> &Self::Output {
        &self.0[index]
    }
}

impl<T: Copy, const N: usize> IndexMut<usize> for FpSimd<T, N> {
    #[inline(always)]
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.0[index]
    }
}

impl<T: Num + Copy, const N: usize> FpSimd<T, N> {
    #[inline(always)]
    pub fn zero() -> Self {
        FpSimd([T::zero(); N])
    }

    /// # Safety
    ///
    /// This function returns an uninitialized object. The caller must ensure
    /// that the returned value is properly initialized before use.
    #[inline(always)]
    #[allow(clippy::uninit_assumed_init)]
    #[allow(invalid_value)]
    pub unsafe fn uninit() -> Self {
        core::mem::MaybeUninit::uninit().assume_init()
    }

    #[inline(always)]
    pub const fn splat(val: T) -> Self {
        FpSimd([val; N])
    }

    #[inline(always)]
    pub fn fill(&mut self, val: T) {
        self.0 = [val; N];
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
    impl<T: Num + Copy> Trait for FpSimd<T, SIMDD> {
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
    impl<T: Num + Copy> Trait<T> for FpSimd<T, SIMDD> {
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
impl<T: NumAssignOps + Copy> Trait for FpSimd<T, SIMDD> {
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

impl<T: Neg<Output = T> + Copy> Neg for FpSimd<T, SIMDD> {
    type Output = Self;
    #[inline(always)]
    fn neg(self) -> Self::Output {
        FpSimd([-self[0], -self[1], -self[2], -self[3], -self[4], -self[5], -self[6], -self[7]])
    }
}

impl<T> FpSimd<T, SIMDD>
where
    T: MulAdd<Output = T> + Copy,
{
    #[inline(always)]
    pub fn add_mul(self, b: FpSimd<T, SIMDD>, c: FpSimd<T, SIMDD>) -> FpSimd<T, SIMDD> {
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

    pub fn fma_from(&mut self, b: FpSimd<T, SIMDD>, c: FpSimd<T, SIMDD>) {
        *self = FpSimd([
            b[0].mul_add(c[0], self[0]),
            b[1].mul_add(c[1], self[1]),
            b[2].mul_add(c[2], self[2]),
            b[3].mul_add(c[3], self[3]),
            b[4].mul_add(c[4], self[4]),
            b[5].mul_add(c[5], self[5]),
            b[6].mul_add(c[6], self[6]),
            b[7].mul_add(c[7], self[7]),
        ])
    }
}

impl f64simd {
    #[inline(always)]
    pub fn is_gto_zero(&self) -> bool {
        for i in 0..SIMDD {
            if self[i].abs() > GTOZERO {
                return false;
            }
        }
        true
    }
}

/* #endregion */

/* #region Blk<T> */

#[repr(align(64))]
#[derive(Clone, Debug, Copy)]
pub struct Blk<T: Copy>(pub [T; BLKSIZE]);

#[allow(non_camel_case_types)]
pub type f64blk = Blk<f64>;
#[allow(non_camel_case_types)]
pub type c64blk = Blk<Complex<f64>>;

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
    pub unsafe fn uninit() -> Self {
        Blk(MaybeUninit::uninit().assume_init())
    }

    #[inline(always)]
    pub const fn splat(val: T) -> Self {
        Blk([val; BLKSIZE])
    }

    #[inline(always)]
    pub fn fill(&mut self, val: T) {
        for i in 0..BLKSIZE {
            self.0[i] = val;
        }
    }

    /// # Safety
    ///
    /// The caller must ensure that `src` has at least `BLKSIZE` elements.
    #[inline(always)]
    pub unsafe fn read_ensure(&mut self, src: &[T]) {
        self.0.copy_from_slice(&src[0..BLKSIZE]);
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
    pub unsafe fn write_ensure(&self, dst: &mut [T]) {
        dst.copy_from_slice(&self.0);
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
    pub const fn as_simdd_slice(&self) -> &[FpSimd<T, SIMDD>; BLKSIMDD] {
        unsafe { transmute(&self.0) }
    }

    #[inline(always)]
    pub fn as_simdd_slice_mut(&mut self) -> &mut [FpSimd<T, SIMDD>; BLKSIMDD] {
        unsafe { transmute(&mut self.0) }
    }

    #[inline(always)]
    pub const fn get_simdd(&self, index: usize) -> FpSimd<T, SIMDD> {
        let slice = self.as_simdd_slice();
        slice[index]
    }

    #[inline(always)]
    pub fn get_simdd_mut(&mut self, index: usize) -> &mut FpSimd<T, SIMDD> {
        let slice = self.as_simdd_slice_mut();
        &mut slice[index]
    }
}

/* #endregion */

pub fn gto_l_iter(l: usize) -> Box<dyn Iterator<Item = (usize, usize, usize)>> {
    match l {
        0 => Box::new(std::iter::once((0, 0, 0))),
        1 => Box::new([(1, 0, 0), (0, 1, 0), (0, 0, 1)].into_iter()),
        2 => Box::new([(2, 0, 0), (1, 1, 0), (1, 0, 1), (0, 2, 0), (0, 1, 1), (0, 0, 2)].into_iter()),
        3 => Box::new(
            [(3, 0, 0), (2, 1, 0), (2, 0, 1), (1, 2, 0), (1, 1, 1), (1, 0, 2), (0, 3, 0), (0, 2, 1), (0, 1, 2), (0, 0, 3)].into_iter(),
        ),
        4 => Box::new(
            [
                (4, 0, 0),
                (3, 1, 0),
                (3, 0, 1),
                (2, 2, 0),
                (2, 1, 1),
                (2, 0, 2),
                (1, 3, 0),
                (1, 2, 1),
                (1, 1, 2),
                (1, 0, 3),
                (0, 4, 0),
                (0, 3, 1),
                (0, 2, 2),
                (0, 1, 3),
                (0, 0, 4),
            ]
            .into_iter(),
        ),
        5 => Box::new(
            [
                (5, 0, 0),
                (4, 1, 0),
                (4, 0, 1),
                (3, 2, 0),
                (3, 1, 1),
                (3, 0, 2),
                (2, 3, 0),
                (2, 2, 1),
                (2, 1, 2),
                (2, 0, 3),
                (1, 4, 0),
                (1, 3, 1),
                (1, 2, 2),
                (1, 1, 3),
                (1, 0, 4),
                (0, 5, 0),
                (0, 4, 1),
                (0, 3, 2),
                (0, 2, 3),
                (0, 1, 4),
                (0, 0, 5),
            ]
            .into_iter(),
        ),
        _ => Box::new((0..=l).flat_map(move |lx| (0..=(l - lx)).map(move |ly| (lx, ly, l - lx - ly)))),
    }
}

pub fn gto_nabla1_simdd(f1: &mut [[f64simd; 3]], f0: &[[f64simd; 3]], l: usize, alpha: f64) {
    let a2 = FpSimd::<f64>::splat(-2.0 * alpha);
    // first derivative
    f1[0][X] = a2 * f0[1][X];
    f1[0][Y] = a2 * f0[1][Y];
    f1[0][Z] = a2 * f0[1][Z];
    // recursive derivatives
    for i in 1..=l {
        let i_f64 = FpSimd::<f64>::splat(i as f64);
        f1[i][X] = i_f64 * f0[i - 1][X] + a2 * f0[i + 1][X];
        f1[i][Y] = i_f64 * f0[i - 1][Y] + a2 * f0[i + 1][Y];
        f1[i][Z] = i_f64 * f0[i - 1][Z] + a2 * f0[i + 1][Z];
    }
}

pub fn gto_x1_simdd(f1: &mut [[f64simd; 3]], f0: &[[f64simd; 3]], l: usize, ri: [f64; 3]) {
    let ri_x = f64simd::splat(ri[X]);
    let ri_y = f64simd::splat(ri[Y]);
    let ri_z = f64simd::splat(ri[Z]);
    for i in 0..=l {
        f1[i][X] = f0[i][X].add_mul(ri_x, f0[i + 1][X]);
        f1[i][Y] = f0[i][Y].add_mul(ri_y, f0[i + 1][Y]);
        f1[i][Z] = f0[i][Z].add_mul(ri_z, f0[i + 1][Z]);
    }
}

pub fn gto_r_simdd(f1: &mut [[f64simd; 3]], f0: &[[f64simd; 3]], l: usize) {
    for i in 0..=l {
        f1[i][X] = f0[i + 1][X];
        f1[i][Y] = f0[i + 1][Y];
        f1[i][Z] = f0[i + 1][Z];
    }
}

pub fn gto_prim_exp(
    // arguments
    eprim: &mut [f64blk],
    coord: &[f64blk; 3],
    alpha: &[f64],
    fac: f64,
    // dimensions
    nprim: usize,
) {
    let mut rr = unsafe { f64blk::uninit() };
    for i in 0..BLKSIZE {
        rr[i] = coord[0][i] * coord[0][i] + coord[1][i] * coord[1][i] + coord[2][i] * coord[2][i];
    }

    for p in 0..nprim {
        for g in 0..BLKSIZE {
            let arr = alpha[p] * rr[g];
            eprim[p][g] = fac * (-arr).exp();
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
