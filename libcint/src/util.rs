use duplicate::duplicate_item;

use crate::prelude::*;

/* #region shls_slice parsing */

/// Convenient structure for parsing shells slice.
pub struct ShlsSlice {
    pub shls_slice: Vec<[usize; 2]>,
}

#[duplicate_item(
    T                   N;
    [Vec<[usize; 2]>]   [];
    [&Vec<[usize; 2]>]  [];
    [&[[usize; 2]]]     [];
    [[[usize; 2]; N]]   [const N: usize];
    [&[[usize; 2]; N]]  [const N: usize];
)]
impl<N> From<T> for ShlsSlice {
    fn from(shls_slice: T) -> Self {
        let shls_slice = shls_slice.to_vec();
        Self { shls_slice }
    }
}

#[duplicate_item(
    T                   N;
    [Vec<[isize; 2]>]   [];
    [&Vec<[isize; 2]>]  [];
    [&[[isize; 2]]]     [];
    [[[isize; 2]; N]]   [const N: usize];
    [&[[isize; 2]; N]]  [const N: usize];
    [Vec<[i32; 2]>]   [];
    [&Vec<[i32; 2]>]  [];
    [&[[i32; 2]]]     [];
    [[[i32; 2]; N]]   [const N: usize];
    [&[[i32; 2]; N]]  [const N: usize];
    [Vec<[i64; 2]>]   [];
    [&Vec<[i64; 2]>]  [];
    [&[[i64; 2]]]     [];
    [[[i64; 2]; N]]   [const N: usize];
    [&[[i64; 2]; N]]  [const N: usize];
)]
impl<N> From<T> for ShlsSlice {
    fn from(shls_slice: T) -> Self {
        let shls_slice = shls_slice.iter().map(|arr| [arr[0] as usize, arr[1] as usize]).collect();
        Self { shls_slice }
    }
}

// This implementation will work when using `None` as input.
impl From<Option<&[[usize; 2]]>> for ShlsSlice {
    fn from(shls_slice: Option<&[[usize; 2]]>) -> Self {
        match shls_slice {
            Some(slice) => slice.to_vec().into(),
            None => ShlsSlice { shls_slice: vec![] },
        }
    }
}

#[duplicate_item(
    T                   N;
    [Vec<Vec<usize>>]   [];
    [&Vec<Vec<usize>>]  [];
    [&[Vec<usize>]]     [];
    [[Vec<usize>; N]]   [const N: usize];
    [&[Vec<usize>; N]]  [const N: usize];
)]
impl<N> From<T> for ShlsSlice {
    fn from(shls_slice: T) -> Self {
        // Convert Vec<Vec<usize>> to Vec<[usize; 2]>
        let shls_slice: Vec<[usize; 2]> = shls_slice
            .iter()
            .map(|vec| {
                if vec.len() == 2 {
                    [vec[0], vec[1]]
                } else {
                    panic!("Each inner vector must have exactly two elements.");
                }
            })
            .collect();
        Self { shls_slice }
    }
}

impl AsRef<[[usize; 2]]> for ShlsSlice {
    fn as_ref(&self) -> &[[usize; 2]] {
        self.shls_slice.as_slice()
    }
}

/* #endregion */

/* #region mathematical computation */

/// Computes the Gaussian integral:
///
/// $$
/// \int_0^{+\inf} x^n e^{-\alpha x^2} dx = \frac{n!}{2\alpha^{(n+1)/2}}
/// $$
pub fn gaussian_int(n: f64, alpha: f64) -> f64 {
    let n1 = (n + 1.0) / 2.0;
    libm::tgamma(n1) / (2.0 * alpha.powf(n1))
}

/* #endregion */

/* #region unsafe cast */

#[allow(clippy::mut_from_ref)]
pub(crate) unsafe fn cast_mut_slice<T>(slc: &[T]) -> &mut [T] {
    let len = slc.len();
    let ptr = slc.as_ptr() as *mut T;
    unsafe { std::slice::from_raw_parts_mut(ptr, len) }
}

/* #endregion */

/* #region indices computation */

#[inline(always)]
pub(crate) fn unravel_s2_indices(x: usize) -> [usize; 2] {
    let j = (((x * 2 + 1) as f64).sqrt() - 0.5) as usize;
    let i = x - j * (j + 1) / 2;
    [i, j]
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

/* #endregion */

/* #region integral block copy */

#[inline(always)]
pub(crate) fn copy_3d_s1<T>(
    out: &mut [T],
    out_offsets: &[usize; 3],
    out_shape: &[usize; 3],
    buf: &[T],
    buf_shape: &[usize; 3],
) where
    T: Copy,
{
    for c in 0..buf_shape[2] {
        for j in 0..buf_shape[1] {
            let out_indices = [out_offsets[0], out_offsets[1] + j, out_offsets[2] + c];
            let buf_indices = [0, j, c];
            let out_start = get_f_index_3d(&out_indices, out_shape);
            let buf_start = get_f_index_3d(&buf_indices, buf_shape);
            let out_stop = out_start + buf_shape[0];
            let buf_stop = buf_start + buf_shape[0];
            out[out_start..out_stop].copy_from_slice(&buf[buf_start..buf_stop]);
        }
    }
}

#[inline(always)]
pub(crate) fn copy_4d_s1<T>(
    out: &mut [T],
    out_offsets: &[usize; 4],
    out_shape: &[usize; 4],
    buf: &[T],
    buf_shape: &[usize; 4],
) where
    T: Copy,
{
    for c in 0..buf_shape[3] {
        for k in 0..buf_shape[2] {
            for j in 0..buf_shape[1] {
                let out_indices =
                    [out_offsets[0], out_offsets[1] + j, out_offsets[2] + k, out_offsets[3] + c];
                let buf_indices = [0, j, k, c];
                let out_start = get_f_index_4d(&out_indices, out_shape);
                let buf_start = get_f_index_4d(&buf_indices, buf_shape);
                let out_stop = out_start + buf_shape[0];
                let buf_stop = buf_start + buf_shape[0];
                out[out_start..out_stop].copy_from_slice(&buf[buf_start..buf_stop]);
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_5d_s1<T>(
    out: &mut [T],
    out_offsets: &[usize; 5],
    out_shape: &[usize; 5],
    buf: &[T],
    buf_shape: &[usize; 5],
) where
    T: Copy,
{
    for c in 0..buf_shape[4] {
        for l in 0..buf_shape[3] {
            for k in 0..buf_shape[2] {
                for j in 0..buf_shape[1] {
                    let out_indices = [
                        out_offsets[0],
                        out_offsets[1] + j,
                        out_offsets[2] + k,
                        out_offsets[3] + l,
                        out_offsets[4] + c,
                    ];
                    let buf_indices = [0, j, k, l, c];
                    let out_start = get_f_index_5d(&out_indices, out_shape);
                    let buf_start = get_f_index_5d(&buf_indices, buf_shape);
                    let out_stop = out_start + buf_shape[0];
                    let buf_stop = buf_start + buf_shape[0];
                    out[out_start..out_stop].copy_from_slice(&buf[buf_start..buf_stop]);
                }
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_3d_s2ij<T>(
    out: &mut [T],
    out_offsets: &[usize; 3],
    out_shape: &[usize; 2],
    buf: &[T],
    buf_shape: &[usize; 3],
) where
    T: Copy,
{
    if out_offsets[0] != out_offsets[1] {
        // offdiag parts
        for c in 0..buf_shape[2] {
            for j in 0..buf_shape[1] {
                let out_indices = [out_offsets[0], out_offsets[1] + j, out_offsets[2] + c];
                let buf_indices = [0, j, c];
                let out_start = get_f_index_3d_s2ij(&out_indices, out_shape);
                let buf_start = get_f_index_3d(&buf_indices, buf_shape);
                let out_stop = out_start + buf_shape[0];
                let buf_stop = buf_start + buf_shape[0];
                out[out_start..out_stop].copy_from_slice(&buf[buf_start..buf_stop]);
            }
        }
    } else {
        // diag parts
        for c in 0..buf_shape[2] {
            for j in 0..buf_shape[1] {
                let out_indices = [out_offsets[0], out_offsets[1] + j, out_offsets[2] + c];
                let buf_indices = [0, j, c];
                let out_start = get_f_index_3d_s2ij(&out_indices, out_shape);
                let buf_start = get_f_index_3d(&buf_indices, buf_shape);
                let out_stop = out_start + j + 1;
                let buf_stop = buf_start + j + 1;
                out[out_start..out_stop].copy_from_slice(&buf[buf_start..buf_stop]);
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_4d_s2ij<T>(
    out: &mut [T],
    out_offsets: &[usize; 4],
    out_shape: &[usize; 3],
    buf: &[T],
    buf_shape: &[usize; 4],
) where
    T: Copy,
{
    if out_offsets[0] != out_offsets[1] {
        // offdiag parts
        for c in 0..buf_shape[3] {
            for k in 0..buf_shape[2] {
                for j in 0..buf_shape[1] {
                    let out_indices = [
                        out_offsets[0],
                        out_offsets[1] + j,
                        out_offsets[2] + k,
                        out_offsets[3] + c,
                    ];
                    let buf_indices = [0, j, k, c];
                    let out_start = get_f_index_4d_s2ij(&out_indices, out_shape);
                    let buf_start = get_f_index_4d(&buf_indices, buf_shape);
                    let out_stop = out_start + buf_shape[0];
                    let buf_stop = buf_start + buf_shape[0];
                    out[out_start..out_stop].copy_from_slice(&buf[buf_start..buf_stop]);
                }
            }
        }
    } else {
        // diag parts
        for c in 0..buf_shape[3] {
            for k in 0..buf_shape[2] {
                for j in 0..buf_shape[1] {
                    let out_indices = [
                        out_offsets[0],
                        out_offsets[1] + j,
                        out_offsets[2] + k,
                        out_offsets[3] + c,
                    ];
                    let buf_indices = [0, j, k, c];
                    let out_start = get_f_index_4d_s2ij(&out_indices, out_shape);
                    let buf_start = get_f_index_4d(&buf_indices, buf_shape);
                    let out_stop = out_start + j + 1;
                    let buf_stop = buf_start + j + 1;
                    out[out_start..out_stop].copy_from_slice(&buf[buf_start..buf_stop]);
                }
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_5d_s2ij<T>(
    out: &mut [T],
    out_offsets: &[usize; 5],
    out_shape: &[usize; 4],
    buf: &[T],
    buf_shape: &[usize; 5],
) where
    T: Copy,
{
    if out_offsets[0] != out_offsets[1] {
        // offdiag parts
        for c in 0..buf_shape[4] {
            for l in 0..buf_shape[3] {
                for k in 0..buf_shape[2] {
                    for j in 0..buf_shape[1] {
                        let out_indices = [
                            out_offsets[0],
                            out_offsets[1] + j,
                            out_offsets[2] + k,
                            out_offsets[3] + l,
                            out_offsets[4] + c,
                        ];
                        let buf_indices = [0, j, k, l, c];
                        let out_start = get_f_index_5d_s2ij(&out_indices, out_shape);
                        let buf_start = get_f_index_5d(&buf_indices, buf_shape);
                        let out_stop = out_start + buf_shape[0];
                        let buf_stop = buf_start + buf_shape[0];
                        out[out_start..out_stop].copy_from_slice(&buf[buf_start..buf_stop]);
                    }
                }
            }
        }
    } else {
        // diag parts
        for c in 0..buf_shape[4] {
            for l in 0..buf_shape[3] {
                for k in 0..buf_shape[2] {
                    for j in 0..buf_shape[1] {
                        let out_indices = [
                            out_offsets[0],
                            out_offsets[1] + j,
                            out_offsets[2] + k,
                            out_offsets[3] + l,
                            out_offsets[4] + c,
                        ];
                        let buf_indices = [0, j, k, l, c];
                        let out_start = get_f_index_5d_s2ij(&out_indices, out_shape);
                        let buf_start = get_f_index_5d(&buf_indices, buf_shape);
                        let out_stop = out_start + j + 1;
                        let buf_stop = buf_start + j + 1;
                        out[out_start..out_stop].copy_from_slice(&buf[buf_start..buf_stop]);
                    }
                }
            }
        }
    }
}

/* #endregion */

/* #region aligned alloc */

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
/// <https://users.rust-lang.org/t/how-can-i-allocate-aligned-memory-in-rust/33293>
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
/// # Safety
///
/// Caller must ensure that the vector is properly initialized before using it.
///
/// This is not a very good function, since `set_len` on uninitialized memory is
/// undefined-behavior (UB).
/// Nevertheless, if `T` is some type of `MaybeUninit`, then this will not UB.
#[inline]
pub unsafe fn aligned_uninitialized_vec<T>(size: usize) -> Vec<T> {
    const MIN_ALIGN: usize = 64; // minimal number of elements in vector to be aligned
    const ALIGNMENT: usize = 64; // 64 bytes alignment (minimal requirement for AVX-512)

    if size == 0 {
        vec![]
    } else if size < MIN_ALIGN {
        unsafe { unaligned_uninitialized_vec(size) }
    } else {
        let sizeof = core::mem::size_of::<T>();
        let pointer = aligned_alloc(size * sizeof, ALIGNMENT);
        if let Some(pointer) = pointer {
            let mut v = unsafe { Vec::from_raw_parts(pointer.as_ptr() as *mut T, size, size) };
            unsafe { v.set_len(size) };
            return v;
        } else {
            panic!("Allocation failed (probably due to out-of-memory, of size {size})")
        }
    }
}

/* #endregion */
