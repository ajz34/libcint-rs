//! Integral crafter for `CInt` instance (the main integral functions).

#![allow(dead_code)]

use crate::prelude::*;
use rayon::prelude::*;

/// Implementation of integral at lower API level.
impl CInt {
    /// Check if the integral can be integrated with the given integrator.
    pub fn check_integratable(
        &self,
        integrator: Box<dyn Integrator>,
        shls_slice: &[[c_int; 2]],
    ) -> Result<(), CIntError> {
        // length of shls_slice must be either 0, or the same value to number of
        // centers of integrator.
        if !shls_slice.is_empty() && shls_slice.len() != integrator.n_center() {
            return Err(CIntError::IntegratorNotAvailable(format!(
                "Integrator {} requires {} centers, but got {}",
                integrator.name(),
                integrator.n_center(),
                shls_slice.len()
            )));
        }

        // Some integrators may not support certain shell types.
        match self.cint_type {
            Spheric => {
                if !integrator.is_sph_available() {
                    return Err(CIntError::IntegratorNotAvailable(format!(
                        "Integrator {} does not support spherical integrals",
                        integrator.name()
                    )));
                }
            },
            Cartesian => {
                if !integrator.is_cart_available() {
                    return Err(CIntError::IntegratorNotAvailable(format!(
                        "Integrator {} does not support Cartesian integrals",
                        integrator.name()
                    )));
                }
            },
            Spinor => {
                if !integrator.is_spinor_available() {
                    return Err(CIntError::IntegratorNotAvailable(format!(
                        "Integrator {} does not support spinor integrals",
                        integrator.name()
                    )));
                }
            },
        }

        // Every shell in shls_slice must be within the range of number of shells.
        let nbas = self.get_nbas() as c_int;
        for (i, shl) in shls_slice.iter().enumerate() {
            if !(0 <= shl[0] && shl[0] <= shl[1] && shl[1] <= nbas) {
                return Err(CIntError::IntegratorNotAvailable(format!(
                    "Shell slice {:?} at index {} is not proper within [{}, {}]",
                    shl, i, 0, nbas
                )));
            }
        }

        Ok(())
    }

    /// Obtain cache size for integral.
    ///
    /// If the shell slice is not known to you currently, just pass empty
    /// `shls_slice = vec![]`, then it should give the maximum cache size
    /// for this molecule/intor.
    ///
    /// # Safety
    ///
    /// - Does not check the validity of the `integrator`.
    /// - Panic if `shls_slice` exceeds the number of shells in the molecule.
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// use libcint::ffi::cint_wrapper::int2e_ip1;
    /// let cint_data = init_h2o_def2_tzvp();
    /// let cache_size = unsafe { cint_data.size_of_cache::<int2e_ip1>(&[]) };
    /// ```
    pub unsafe fn size_of_cache<T>(&self, shls_slice: &[[c_int; 2]]) -> usize
    where
        T: Integrator + Default,
    {
        let integrator = T::default();
        unsafe { self.size_of_cache_dyn(&integrator, shls_slice) }
    }

    /// Obtain cache size for integral.
    ///
    /// If the shell slice is not known to you currently, just pass empty
    /// `shls_slice = vec![]`, then it should give the maximum cache size
    /// for this molecule/intor.
    ///
    /// # Safety
    ///
    /// - Does not check the validity of the `integrator`.
    /// - Panic if `shls_slice` exceeds the number of shells in the molecule.
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// let integrator = CInt::get_integrator("int2e_ip1");
    /// let cache_size = unsafe { cint_data.size_of_cache_dyn(&*integrator, &[]) };
    /// ```
    pub unsafe fn size_of_cache_dyn(
        &self,
        integrator: &dyn Integrator,
        shls_slice: &[[c_int; 2]],
    ) -> usize {
        let natm = self.get_natm() as c_int;
        let nbas = self.get_nbas() as c_int;
        let shls_min = shls_slice.iter().map(|x| x[0]).min().unwrap_or(0);
        let shls_max = shls_slice.iter().map(|x| x[1]).max().unwrap_or(nbas);

        (shls_min..shls_max)
            .into_par_iter()
            .map(|shl| unsafe {
                let shls = [shl; 4];
                match self.cint_type {
                    Spheric => integrator.integral_sph(
                        null_mut(),
                        null(),
                        shls.as_ptr(),
                        self.get_atm_ptr(),
                        natm,
                        self.get_bas_ptr(),
                        nbas,
                        self.get_env_ptr(),
                        null(),
                        null_mut(),
                    ),
                    Cartesian => integrator.integral_cart(
                        null_mut(),
                        null(),
                        shls.as_ptr(),
                        self.get_atm_ptr(),
                        natm,
                        self.get_bas_ptr(),
                        nbas,
                        self.get_env_ptr(),
                        null(),
                        null_mut(),
                    ),
                    Spinor => integrator.integral_spinor(
                        null_mut(),
                        null(),
                        shls.as_ptr(),
                        self.get_atm_ptr(),
                        natm,
                        self.get_bas_ptr(),
                        nbas,
                        self.get_env_ptr(),
                        null(),
                        null_mut(),
                    ),
                }
            })
            .max()
            .unwrap() as usize
    }
}

/* #region util */

#[allow(clippy::mut_from_ref)]
pub(crate) unsafe fn cast_mut_slice<T>(slc: &[T]) -> &mut [T] {
    let len = slc.len();
    let ptr = slc.as_ptr() as *mut T;
    unsafe { std::slice::from_raw_parts_mut(ptr, len) }
}

#[inline(always)]
pub(crate) fn get_f_index_3d(indices: &[usize; 3], shape: &[usize; 3]) -> usize {
    indices[0] + shape[0] * (indices[1] + shape[1] * (indices[2]))
}

#[inline(always)]
pub(crate) fn get_f_index_4d(indices: &[usize; 4], shape: &[usize; 4]) -> usize {
    indices[0] + shape[0] * (indices[1] + shape[1] * (indices[2] + shape[2] * (indices[3])))
}

#[inline(always)]
pub(crate) fn get_f_index_5d(indices: &[usize; 5], shape: &[usize; 5]) -> usize {
    indices[0]
        + shape[0]
            * (indices[1]
                + shape[1] * (indices[2] + shape[2] * (indices[3] + shape[3] * (indices[4]))))
}

#[inline(always)]
pub(crate) fn get_f_index_3d_s2ij(indices: &[usize; 3], shape: &[usize; 2]) -> usize {
    indices[0] + indices[1] * (indices[1] + 1) / 2 + shape[0] * (indices[2])
}

#[inline(always)]
pub(crate) fn get_f_index_4d_s2ij(indices: &[usize; 4], shape: &[usize; 3]) -> usize {
    indices[0]
        + indices[1] * (indices[1] + 1) / 2
        + shape[0] * (indices[2] + shape[1] * (indices[3]))
}

#[inline(always)]
pub(crate) fn get_f_index_5d_s2ij(indices: &[usize; 5], shape: &[usize; 4]) -> usize {
    indices[0]
        + indices[1] * (indices[1] + 1) / 2
        + shape[0] * (indices[2] + shape[1] * (indices[3] + shape[2] * (indices[4])))
}

#[inline(always)]
pub(crate) fn copy_3d_s2ij_offdiag<T>(
    out: &mut [T],
    out_offsets: &[usize; 3],
    out_s2ij_shape: &[usize; 2],
    buf: &[T],
    buf_shape: &[usize; 3],
) where
    T: Copy,
{
    for c in 0..buf_shape[2] {
        for j in 0..buf_shape[1] {
            for i in 0..buf_shape[0] {
                let out_indices = [out_offsets[0] + i, out_offsets[1] + j, out_offsets[2] + c];
                let buf_indices = [i, j, c];
                let out_index = get_f_index_3d_s2ij(&out_indices, out_s2ij_shape);
                let buf_index = get_f_index_3d(&buf_indices, buf_shape);
                out[out_index] = buf[buf_index];
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_3d_s2ij_diag<T>(
    out: &mut [T],
    out_offsets: &[usize; 3],
    out_s2ij_shape: &[usize; 2],
    buf: &[T],
    buf_shape: &[usize; 3],
) where
    T: Copy,
{
    for c in 0..buf_shape[2] {
        for j in 0..buf_shape[1] {
            for i in 0..(j + 1) {
                let out_indices = [out_offsets[0] + i, out_offsets[1] + j, out_offsets[2] + c];
                let buf_indices = [i, j, c];
                let out_index = get_f_index_3d_s2ij(&out_indices, out_s2ij_shape);
                let buf_index = get_f_index_3d(&buf_indices, buf_shape);
                out[out_index] = buf[buf_index];
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_4d_s2ij_offdiag<T>(
    out: &mut [T],
    out_offsets: &[usize; 4],
    out_s2ij_shape: &[usize; 3],
    buf: &[T],
    buf_shape: &[usize; 4],
) where
    T: Copy,
{
    for c in 0..buf_shape[3] {
        for k in 0..buf_shape[2] {
            for j in 0..buf_shape[1] {
                for i in 0..buf_shape[0] {
                    let out_indices = [
                        out_offsets[0] + i,
                        out_offsets[1] + j,
                        out_offsets[2] + k,
                        out_offsets[3] + c,
                    ];
                    let buf_indices = [i, j, k, c];
                    let out_index = get_f_index_4d_s2ij(&out_indices, out_s2ij_shape);
                    let buf_index = get_f_index_4d(&buf_indices, buf_shape);
                    out[out_index] = buf[buf_index];
                }
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_4d_s2ij_diag<T>(
    out: &mut [T],
    out_offsets: &[usize; 4],
    out_s2ij_shape: &[usize; 3],
    buf: &[T],
    buf_shape: &[usize; 4],
) where
    T: Copy,
{
    for c in 0..buf_shape[3] {
        for k in 0..buf_shape[2] {
            for j in 0..buf_shape[1] {
                for i in 0..(j + 1) {
                    let out_indices = [
                        out_offsets[0] + i,
                        out_offsets[1] + j,
                        out_offsets[2] + k,
                        out_offsets[3] + c,
                    ];
                    let buf_indices = [i, j, k, c];
                    let out_index = get_f_index_4d_s2ij(&out_indices, out_s2ij_shape);
                    let buf_index = get_f_index_4d(&buf_indices, buf_shape);
                    out[out_index] = buf[buf_index];
                }
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_5d_s2ij_offdiag<T>(
    out: &mut [T],
    out_offsets: &[usize; 5],
    out_s2ij_shape: &[usize; 4],
    buf: &[T],
    buf_shape: &[usize; 5],
) where
    T: Copy,
{
    for c in 0..buf_shape[4] {
        for l in 0..buf_shape[3] {
            for k in 0..buf_shape[2] {
                for j in 0..buf_shape[1] {
                    for i in 0..buf_shape[0] {
                        let out_indices = [
                            out_offsets[0] + i,
                            out_offsets[1] + j,
                            out_offsets[2] + k,
                            out_offsets[3] + l,
                            out_offsets[4] + c,
                        ];
                        let buf_indices = [i, j, k, l, c];
                        let out_index = get_f_index_5d_s2ij(&out_indices, out_s2ij_shape);
                        let buf_index = get_f_index_5d(&buf_indices, buf_shape);
                        out[out_index] = buf[buf_index];
                    }
                }
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_5d_s2ij_diag<T>(
    out: &mut [T],
    out_offsets: &[usize; 5],
    out_s2ij_shape: &[usize; 4],
    buf: &[T],
    buf_shape: &[usize; 5],
) where
    T: Copy,
{
    for c in 0..buf_shape[4] {
        for l in 0..buf_shape[3] {
            for k in 0..buf_shape[2] {
                for j in 0..buf_shape[1] {
                    for i in 0..(j + 1) {
                        let out_indices = [
                            out_offsets[0] + i,
                            out_offsets[1] + j,
                            out_offsets[2] + k,
                            out_offsets[3] + l,
                            out_offsets[4] + c,
                        ];
                        let buf_indices = [i, j, k, l, c];
                        let out_index = get_f_index_5d_s2ij(&out_indices, out_s2ij_shape);
                        let buf_index = get_f_index_5d(&buf_indices, buf_shape);
                        out[out_index] = buf[buf_index];
                    }
                }
            }
        }
    }
}

/// Create an unaligned uninitialized vector with the given size.
///
/// # Safety
///
/// Caller must ensure that the vector is properly initialized before using it.
///
/// This is not a very good function, since `set_len` on uninitialized memory is
/// undefined-behavior (UB).
/// Nevertheless, if `T` is some type of `MaybeUninit`, then this will not UB.
#[allow(clippy::uninit_vec)]
#[inline]
pub unsafe fn unaligned_uninitialized_vec<T>(size: usize) -> Vec<T> {
    let mut v: Vec<T> = vec![];
    v.try_reserve_exact(size).unwrap();
    unsafe { v.set_len(size) };
    v
}

/// Create an uninitialized vector with the given size and alignment.
///
/// - None: if the size is 0 or allocation fails.
/// - Some: pointer to the allocated memory.
///
/// https://users.rust-lang.org/t/how-can-i-allocate-aligned-memory-in-rust/33293
#[inline]
pub fn aligned_alloc(numbytes: usize, alignment: usize) -> Option<NonNull<()>> {
    extern crate alloc;

    if numbytes == 0 {
        return None;
    }
    let layout = alloc::alloc::Layout::from_size_align(numbytes, alignment).unwrap();
    NonNull::new(unsafe { alloc::alloc::alloc(layout) }).map(|p| p.cast::<()>())
}

/// Create an conditionally aligned uninitialized vector with the given size.
///
/// - `N`: condition for alignment; if `N < size`, then this function will not
///   allocate aligned vector.
///
/// # Safety
///
/// Caller must ensure that the vector is properly initialized before using it.
///
/// This is not a very good function, since `set_len` on uninitialized memory is
/// undefined-behavior (UB).
/// Nevertheless, if `T` is some type of `MaybeUninit`, then this will not UB.
#[inline]
pub unsafe fn aligned_uninitialized_vec<T, const N: usize>(
    size: usize,
    alignment: usize,
) -> Vec<T> {
    const MIN_ALIGN: usize = 64;

    if size == 0 {
        vec![]
    } else if size < N {
        unsafe { unaligned_uninitialized_vec(size) }
    } else {
        let sizeof = core::mem::size_of::<T>();
        let pointer = aligned_alloc(size * sizeof, alignment);
        if let Some(pointer) = pointer {
            let mut v = unsafe { Vec::from_raw_parts(pointer.as_ptr() as *mut T, size, size) };
            unsafe { v.set_len(size) };
            return v;
        } else {
            panic!("Allocation failed (probably due to out-of-memory)")
        }
    }
}

/* #endregion */

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn playground() {
        use crate::ffi::cint_wrapper::int2e_ip1;
        let cint_data = init_h2o_def2_tzvp();
        let integrator = CInt::get_integrator("int2e_ip1");
        let cache_size = unsafe { cint_data.size_of_cache_dyn(&*integrator, &[]) };
        println!("Cache size: {cache_size}");
        let cache_size = unsafe { cint_data.size_of_cache::<int2e_ip1>(&[]) };
        println!("Cache size: {cache_size}");
    }
}
