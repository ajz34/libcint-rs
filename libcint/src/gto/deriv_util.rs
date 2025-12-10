use crate::gto::prelude_dev::*;
use num::traits::{MulAdd, NumAssignOps};
use num::Num;

/* #region simple f64x8 */

// TODO: use fearless_simd or pulp, they are more heavy but more complete SIMD
// abstraction crates

/// (dev) GTO internal SIMD type.
///
/// In most cases, we use f64x8 as the SIMD type, which corresponds to AVX-512.
///
/// This type implements basic arithmetic operations and some utility functions,
/// such as arithmetics, fmadd, map, splat, etc.
///
/// Please note this is only a simple implementation. This does not cover any
/// target of protable SIMD propose ([#86656](https://github.com/rust-lang/rust/issues/86656)).
///
/// To fully utilize SIMD capabilities, you need to compile by `RUSTFLAGS="-C
/// target-cpu=native"` or similar flags.
#[repr(align(64))]
#[derive(Clone, Debug, Copy)]
pub struct FpSimd<T: Copy, const N: usize = SIMDD>(pub [T; N]);

/// Type alias for f64x8 SIMD type [`FpSimd<f64, 8>`].
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
    /// Returns a SIMD object with all lanes set to zero.
    #[inline(always)]
    pub fn zero() -> Self {
        FpSimd([T::zero(); N])
    }

    /// Returns an uninitialized SIMD object.
    ///
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

    /// Returns a SIMD object with all lanes set to `val`.
    #[inline(always)]
    pub const fn splat(val: T) -> Self {
        FpSimd([val; N])
    }

    /// Sets all lanes to `val`.
    #[inline(always)]
    pub fn fill(&mut self, val: T) {
        self.0 = [val; N];
    }

    /// Applies function `f` to each lane and returns a new SIMD object.
    #[inline]
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
    /// Performs fused multiply-add: `self * b + c`.
    ///
    /// This is similar function to [`MulAdd::mul_add`].
    #[inline(always)]
    pub fn mul_add(self, b: FpSimd<T, SIMDD>, c: FpSimd<T, SIMDD>) -> FpSimd<T, SIMDD> {
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

    /// Performs fused multiply-add: `self = self + b * c`.
    ///
    /// Note that the order of multiplication and addition is different from
    /// [`FpSimd::mul_add`].
    #[inline(always)]
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
    /// Checks if all lanes are smaller than the GTO zero cutoff [`GTOZERO`].
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

/// (dev) GTO internal block type.
///
/// This type represents a block of [`BLKSIZE`] elements of type `T`, aligned to
/// 64 bytes (AVX-512 alignment).
///
/// This type is a basic building block for GTO grid evaluations. Grids will be
/// always batched by [`BLKSIZE`].
///
/// This [`BLKSIZE`] value (48) should be multiple of [`SIMDD`] (8). It should
/// be better to be multiple of 16 because of microkernel in matmul usually
/// requires 2 lanes of SIMD (for AVX-512, it is 16 f64). Currently this value
/// is fixed to 48, in that for the most-used deriv1 case, the processed grids
/// per function call is $(n_\mathrm{comp}, n_\mathrm{ctr} \times
/// n_\mathrm{cart}(l), n_\mathrm{grids})$ will usually be up to (4, 15, 48) or
/// 22.5 KB, which fits in L1d cache (32 KB in most micro-architectures)
/// together with other data.
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
    /// Returns an uninitialized bulk.
    ///
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

    /// Returns a block with all elements set to `val`.
    #[inline(always)]
    pub const fn splat(val: T) -> Self {
        Blk([val; BLKSIZE])
    }

    /// Sets all elements to `val`.
    #[inline(always)]
    pub fn fill(&mut self, val: T) {
        for i in 0..BLKSIZE {
            self.0[i] = val;
        }
    }

    /// Reads data from `src` slice into the block, ensuring [`BLKSIZE`]
    /// elements.
    ///
    /// # Safety
    ///
    /// The caller must ensure that `src` has at least `BLKSIZE` elements.
    #[inline(always)]
    pub unsafe fn read_ensure(&mut self, src: &[T]) {
        self.0.copy_from_slice(&src[0..BLKSIZE]);
    }

    /// Reads data from `src` slice into the block, at most [`BLKSIZE`]
    /// elements.
    #[inline]
    pub fn read(&mut self, src: &[T]) {
        let len_slc = if src.len() < BLKSIZE { src.len() } else { BLKSIZE };
        if len_slc >= BLKSIZE {
            unsafe { self.read_ensure(src) }
        } else {
            self.0.copy_from_slice(src);
        }
    }

    /// Writes data from the block into `dst` slice, ensuring [`BLKSIZE`]
    /// elements.
    ///
    /// # Safety
    ///
    /// The caller must ensure that `dst` has at least `BLKSIZE` elements.
    #[inline(always)]
    pub unsafe fn write_ensure(&self, dst: &mut [T]) {
        dst.copy_from_slice(&self.0);
    }

    /// Writes data from the block into `dst` slice, at most [`BLKSIZE`]
    /// elements.
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
    /// Returns a slice of SIMD (double float) blocks view of this block.
    #[inline(always)]
    pub const fn as_simdd_slice(&self) -> &[FpSimd<T, SIMDD>; BLKSIMDD] {
        unsafe { transmute(&self.0) }
    }

    /// Returns a mutable slice of SIMD (double float) blocks view of this
    /// block.
    #[inline(always)]
    pub fn as_simdd_slice_mut(&mut self) -> &mut [FpSimd<T, SIMDD>; BLKSIMDD] {
        unsafe { transmute(&mut self.0) }
    }

    /// Gets the SIMD (double float) block at `index`.
    #[inline(always)]
    pub const fn get_simdd(&self, index: usize) -> FpSimd<T, SIMDD> {
        let slice = self.as_simdd_slice();
        slice[index]
    }

    /// Gets a mutable reference to the SIMD (double float) block at `index`.
    #[inline(always)]
    pub fn get_simdd_mut(&mut self, index: usize) -> &mut FpSimd<T, SIMDD> {
        let slice = self.as_simdd_slice_mut();
        &mut slice[index]
    }
}

/* #endregion */

/// Returns an iterator over ($l_x$, $l_y$, $l_z$) tuples for given angular
/// momentum `l`.
///
/// Iteration order: $l_x$ varies fastest, then $l_y$, then $l_z$.
///
/// This function supports up to `l = 5` (g) with pre-defined tuples for
/// efficiency.
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

/// Computes first derivatives of polynomial part of GTO in SIMD blocks.
///
/// # Formula of Gaussian
///
/// As an example, cartesian gaussian's derivative to $x$ is computed by:
///
/// $$
/// \frac{\partial}{\partial x} | l_x, l_y, l_z, \alpha, \bm r \rangle = l_x |
/// l_x, l_y, l_z, \alpha, \bm r \rangle - 2 \alpha | l_x + 1, l_y, l_z, \alpha,
/// \bm r \rangle
/// $$
///
/// where the recursion starts from:
///
/// $$
/// \frac{\partial}{\partial x} | 0, l_y, l_z, \alpha, \bm r \rangle = - 2
/// \alpha | 1, l_y, l_z, \alpha, \bm r \rangle
/// $$
///
/// The same applies to $y$ and $z$ directions.
///
/// See [`gto`](crate::gto) module documentation for notation details.
///
/// # Formula of recursion relations of polynomial part
///
/// This function evaluates the polynomial part of cartesian GTO (the whole $|
/// \bm l, \alpha, \bm r \rangle$ without $\exp(- \alpha r^2)$). The recursion
/// correlation of this polynomial (along each coordinate component) is given
/// by:
///
/// $$
/// \begin{alignat*}{5}
/// F_0^{(1)} &= -2 \alpha F_1^{(0)} &\quad& (\text{initial}) \\\\
/// F_l^{(1)} &= l F_{l-1}^{(0)} - 2 \alpha F_{l+1}^{(0)} &\quad& (l \geq 1)
/// \end{alignat*}
/// $$
///
/// Please note that we are evaluating the polynomial part at each grid point.
/// So in the program, we are actually handling $F_{l, g}^{(0)}$ and $F_{l,
/// g}^{(1)}$ as two-dimension tensors.
///
/// # Argument Table
///
/// | variable | formula | dimension | shape | notes |
/// |--|--|--|--|--|
/// | `f1` | $[F_{l_x, g}^{(1)}, F_{l_y, g}^{(1)}, F_{l_z, g}^{(1)}]$ | $(l + 1, t, g)$ | `(l + 1, 3, BLKSIMDD)` | output buffer for first derivatives |
/// | `f0` | $[F_{l_x, g}^{(0)}, F_{l_y, g}^{(0)}, F_{l_z, g}^{(0)}]$ | $(l + 2, t, g)$ | `(l + 2, 3, BLKSIMDD)` | input buffer for function values |
/// | `l` | $l$ | | scalar | angular momentum of current GTO shell |
/// | `alpha` | $\alpha$ | | scalar | exponent value of current GTO shell |
///
/// # PySCF equivalent
///
/// `libcgto.so`: `void GTOnabla1`
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

/// Computes elementwise product of polynomial part of GTO with 3-component
/// vector $\bm R$ in SIMD blocks.
///
/// # Formula of Gaussian
///
/// This usually evaluates the grid-point coordinate, which is composition of
/// electronic coordinate (centered by atom) $\bm r$ and the atomic center
/// position $\bm R$:
///
/// $$
/// (\bm r + \bm R) | l_x, l_y, l_z, \alpha, \bm r \rangle
/// $$
///
/// Take the example of $x$ component:
///
/// $$
/// \begin{align*}
/// (x + R_x) | l_x, l_y, l_z, \alpha, \bm r \rangle
/// &= x | l_x, l_y, l_z, \alpha, \bm r \rangle + R_x | l_x, l_y,
/// l_z, \alpha, \bm r \rangle
/// \\\\
/// &= R_x | l_x, l_y, l_z, \alpha, \bm r \rangle + | l_x + 1, l_y,
/// l_z, \alpha, \bm r \rangle \end{align*}
/// $$
///
/// The same applies to $y$ and $z$ directions.
///
/// See [`gto`](crate::gto) module documentation for notation details.
///
/// # Formula of recursion relations of polynomial part
///
/// This function evaluates the polynomial part of cartesian GTO (the whole $|
/// \bm l, \alpha, \bm r \rangle$ without $\exp(- \alpha r^2)$). The recursion
/// correlation of this polynomial (along each coordinate component) is given
/// by:
///
/// $$
/// F_{l_t}^\texttt{1} = R_t F_{l_t}^\texttt{0} +
/// F_{l_t+1}^\texttt{0} $$
///
/// # Argument Table
///
/// | variable | formula | dimension | shape | notes |
/// |--|--|--|--|--|
/// | `f1` | $[F_{l_x, g}^\texttt{1}, F_{l_y, g}^\texttt{1}, F_{l_z, g}^\texttt{1}]$ | $(l + 1, t, g)$ | `(l + 1, 3, BLKSIMDD)` | output buffer for first derivatives |
/// | `f0` | $[F_{l_x, g}^\texttt{0}, F_{l_y, g}^\texttt{0}, F_{l_z, g}^\texttt{0}]$ | $(l + 2, t, g)$ | `(l + 2, 3, BLKSIMDD)` | input buffer for function values |
/// | `l` | $l$ | | scalar | angular momentum of current GTO shell |
/// | `ri` | $\bm R$ | $(t)$ | `(3)` | center position of current GTO shell |
///
/// # PySCF equivalent
///
/// `libcgto.so`: `void GTOx1`
pub fn gto_x1_simdd(f1: &mut [[f64simd; 3]], f0: &[[f64simd; 3]], l: usize, ri: [f64; 3]) {
    let ri_x = f64simd::splat(ri[X]);
    let ri_y = f64simd::splat(ri[Y]);
    let ri_z = f64simd::splat(ri[Z]);
    for i in 0..=l {
        f1[i][X] = f0[i][X].mul_add(ri_x, f0[i + 1][X]);
        f1[i][Y] = f0[i][Y].mul_add(ri_y, f0[i + 1][Y]);
        f1[i][Z] = f0[i][Z].mul_add(ri_z, f0[i + 1][Z]);
    }
}

/// Computes polynomial part of GTO multiplied by electronic coordinate $\bm r$
/// in SIMD blocks.
///
/// # Formula of Gaussian
///
/// As an example, cartesian gaussian's multiplication by $x$ is computed by:
///
/// $$
/// x | l_x, l_y, l_z, \alpha, \bm r \rangle = | l_x + 1, l_y, l_z, \alpha,
/// \bm r \rangle
/// $$
///
/// The same applies to $y$ and $z$ directions.
///
/// See [`gto`](crate::gto) module documentation for notation details.
///
/// # Formula of recursion relations of polynomial part
///
/// This function evaluates the polynomial part of cartesian GTO (the whole $|
/// \bm l, \alpha, \bm r \rangle$ without $\exp(- \alpha r^2)$). The recursion
/// correlation of this polynomial (along each coordinate component) is given
/// by:
///
/// $$
/// F_l^\texttt{1} = F_{l+1}^\texttt{0}
/// $$
///
/// # Argument Table
///
/// | variable | formula | dimension | shape | notes |
/// |--|--|--|--|--|
/// | `f1` | $[F_{l_x, g}^\texttt{1}, F_{l_y, g}^\texttt{1}, F_{l_z, g}^\texttt{1}]$ | $(l + 1, t, g)$ | `(l + 1, 3, BLKSIMDD)` | output buffer for first derivatives |
/// | `f0` | $[F_{l_x, g}^\texttt{0}, F_{l_y, g}^\texttt{0}, F_{l_z, g}^\texttt{0}]$ | $(l + 2, t, g)$ | `(l + 2, 3, BLKSIMDD)` | input buffer for function values |
/// | `l` | $l$ | | scalar | angular momentum of current GTO shell |
///
/// # PySCF equivalent
///
/// `grid_ao_drv.h`: macro `GTO_R_I`
pub fn gto_r_simdd(f1: &mut [[f64simd; 3]], f0: &[[f64simd; 3]], l: usize) {
    for i in 0..=l {
        f1[i][X] = f0[i + 1][X];
        f1[i][Y] = f0[i + 1][Y];
        f1[i][Z] = f0[i + 1][Z];
    }
}

/// Computes GTO exponents without angular momentum and normalization in SIMD
/// blocks of a shell.
///
/// # Formula
///
/// This function evaluates only the exponent part (without polynomial part) of
/// GTO.
///
/// $$
/// | \bm 0, \alpha, \bm r \rangle = \exp(- \alpha r^2)
/// $$
///
/// where $r^2 = x^2 + y^2 + z^2$.
///
/// # Argument Table
///
/// | variable | formula | dimension | shape | notes |
/// |--|--|--|--|--|
/// | `eprim` | $\vert \bm 0, \alpha_p, \bm r_g \rangle$ | $(p, g)$ | `(nprim, BLKSIZE)` | output buffer for GTO exponents |
/// | `coord` | $\bm r_g$ | $(t, g)$ | `(3, BLKSIZE)` | grid coordinates |
/// | `alpha` | $\alpha_p$ | $(p)$ | `(nprim)` | exponent values for primitives |
/// | `fac` | | | scalar | prefactor for GTO exponents |
/// | `nprim` | | | dim-size | number of primitives (index of $p$) |
///
/// # See also
///
/// This function is a low-level function to implement the
/// [`GtoEvalAPI::gto_exp`] trait method.
///
/// # PySCF equivalent
///
/// `libcgto.so`: `int GTOprim_exp`
pub fn gto_prim_exp(eprim: &mut [f64blk], coord: &[f64blk; 3], alpha: &[f64], fac: f64, nprim: usize) {
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
