//! Utility functions and structs (not for users).

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

/* #region indices computation (col-major) */

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
    indices[0] + shape[0] * (indices[1] + shape[1] * (indices[2] + shape[2] * (indices[3] + shape[3] * (indices[4]))))
}

#[inline(always)]
pub(crate) fn get_f_index_3d_s2ij(indices: &[usize; 3], shape: &[usize; 2]) -> usize {
    indices[0] + indices[1] * (indices[1] + 1) / 2 + shape[0] * (indices[2])
}

#[inline(always)]
pub(crate) fn get_f_index_4d_s2ij(indices: &[usize; 4], shape: &[usize; 3]) -> usize {
    indices[0] + indices[1] * (indices[1] + 1) / 2 + shape[0] * (indices[2] + shape[1] * (indices[3]))
}

#[inline(always)]
pub(crate) fn get_f_index_5d_s2ij(indices: &[usize; 5], shape: &[usize; 4]) -> usize {
    indices[0] + indices[1] * (indices[1] + 1) / 2 + shape[0] * (indices[2] + shape[1] * (indices[3] + shape[2] * (indices[4])))
}

#[inline(always)]
pub(crate) fn get_f_index_5d_s2kl(indices: &[usize; 5], shape: &[usize; 4]) -> usize {
    let [i0, i1, i2, i3, i4] = indices;
    let [s0, s1, s2, _] = shape;
    let indices_kl = i2 + i3 * (i3 + 1) / 2;
    i0 + s0 * (i1 + s1 * (indices_kl + s2 * i4))
}

#[inline(always)]
pub(crate) fn get_f_index_5d_s4(indices: &[usize; 5], shape: &[usize; 3]) -> usize {
    let [i0, i1, i2, i3, i4] = indices;
    let [s0, s1, _] = shape;
    let indices_ij = i0 + i1 * (i1 + 1) / 2;
    let indices_kl = i2 + i3 * (i3 + 1) / 2;
    indices_ij + s0 * (indices_kl + s1 * i4)
}

#[inline(always)]
pub(crate) fn get_f_index_5d_s8(indices: &[usize; 5], shape: &[usize; 2]) -> usize {
    let [i0, i1, i2, i3, i4] = indices;
    let [s0, _] = shape;
    let indices_ij = i0 + i1 * (i1 + 1) / 2;
    let indices_kl = i2 + i3 * (i3 + 1) / 2;
    indices_ij + indices_kl * (indices_kl + 1) / 2 + s0 * i4
}

/* #endregion */

/* #region integral block copy (col-major) */

#[inline]
pub(crate) fn copy_f_3d_s1<T>(out: &mut [T], out_offsets: &[usize; 3], out_shape: &[usize; 3], buf: &[T], buf_shape: &[usize; 3])
where
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

#[inline]
pub(crate) fn copy_f_4d_s1<T>(out: &mut [T], out_offsets: &[usize; 4], out_shape: &[usize; 4], buf: &[T], buf_shape: &[usize; 4])
where
    T: Copy,
{
    for c in 0..buf_shape[3] {
        for k in 0..buf_shape[2] {
            for j in 0..buf_shape[1] {
                let out_indices = [out_offsets[0], out_offsets[1] + j, out_offsets[2] + k, out_offsets[3] + c];
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

#[inline]
pub(crate) fn copy_f_5d_s1<T>(out: &mut [T], out_offsets: &[usize; 5], out_shape: &[usize; 5], buf: &[T], buf_shape: &[usize; 5])
where
    T: Copy,
{
    for c in 0..buf_shape[4] {
        for l in 0..buf_shape[3] {
            for k in 0..buf_shape[2] {
                for j in 0..buf_shape[1] {
                    let out_indices = [out_offsets[0], out_offsets[1] + j, out_offsets[2] + k, out_offsets[3] + l, out_offsets[4] + c];
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

#[inline]
pub(crate) fn copy_f_3d_s2ij<T>(out: &mut [T], out_offsets: &[usize; 3], out_shape: &[usize; 2], buf: &[T], buf_shape: &[usize; 3])
where
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

#[inline]
pub(crate) fn copy_f_4d_s2ij<T>(out: &mut [T], out_offsets: &[usize; 4], out_shape: &[usize; 3], buf: &[T], buf_shape: &[usize; 4])
where
    T: Copy,
{
    if out_offsets[0] != out_offsets[1] {
        // offdiag parts
        for c in 0..buf_shape[3] {
            for k in 0..buf_shape[2] {
                for j in 0..buf_shape[1] {
                    let out_indices = [out_offsets[0], out_offsets[1] + j, out_offsets[2] + k, out_offsets[3] + c];
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
                    let out_indices = [out_offsets[0], out_offsets[1] + j, out_offsets[2] + k, out_offsets[3] + c];
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

#[inline]
pub(crate) fn copy_f_5d_s2ij<T>(out: &mut [T], out_offsets: &[usize; 5], out_shape: &[usize; 4], buf: &[T], buf_shape: &[usize; 5])
where
    T: Copy,
{
    if out_offsets[0] != out_offsets[1] {
        // offdiag parts
        for c in 0..buf_shape[4] {
            for l in 0..buf_shape[3] {
                for k in 0..buf_shape[2] {
                    for j in 0..buf_shape[1] {
                        let out_indices = [out_offsets[0], out_offsets[1] + j, out_offsets[2] + k, out_offsets[3] + l, out_offsets[4] + c];
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
                        let out_indices = [out_offsets[0], out_offsets[1] + j, out_offsets[2] + k, out_offsets[3] + l, out_offsets[4] + c];
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

#[inline]
pub(crate) fn copy_f_5d_s2kl<T>(out: &mut [T], out_offsets: &[usize; 5], out_shape: &[usize; 4], buf: &[T], buf_shape: &[usize; 5])
where
    T: Copy,
{
    let is_diag_kl = out_offsets[2] == out_offsets[3];
    for c in 0..buf_shape[4] {
        for l in 0..buf_shape[3] {
            let k_max = if is_diag_kl { l + 1 } else { buf_shape[2] };
            for k in 0..k_max {
                for j in 0..buf_shape[1] {
                    let out_indices = [out_offsets[0], out_offsets[1] + j, out_offsets[2] + k, out_offsets[3] + l, out_offsets[4] + c];
                    let buf_indices = [0, j, k, l, c];
                    let out_start = get_f_index_5d_s2kl(&out_indices, out_shape);
                    let buf_start = get_f_index_5d(&buf_indices, buf_shape);
                    let out_stop = out_start + buf_shape[0];
                    let buf_stop = buf_start + buf_shape[0];
                    out[out_start..out_stop].copy_from_slice(&buf[buf_start..buf_stop]);
                }
            }
        }
    }
}

#[inline]
pub(crate) fn copy_f_5d_s4<T>(out: &mut [T], out_offsets: &[usize; 5], out_shape: &[usize; 3], buf: &[T], buf_shape: &[usize; 5])
where
    T: Copy,
{
    let is_diag_ij = out_offsets[0] == out_offsets[1];
    let is_diag_kl = out_offsets[2] == out_offsets[3];
    for c in 0..buf_shape[4] {
        for l in 0..buf_shape[3] {
            let k_max = if is_diag_kl { l + 1 } else { buf_shape[2] };
            for k in 0..k_max {
                for j in 0..buf_shape[1] {
                    let i_max = if is_diag_ij { j + 1 } else { buf_shape[0] };
                    let out_indices = [out_offsets[0], out_offsets[1] + j, out_offsets[2] + k, out_offsets[3] + l, out_offsets[4] + c];
                    let buf_indices = [0, j, k, l, c];
                    let out_start = get_f_index_5d_s4(&out_indices, out_shape);
                    let buf_start = get_f_index_5d(&buf_indices, buf_shape);
                    let out_stop = out_start + i_max;
                    let buf_stop = buf_start + i_max;
                    out[out_start..out_stop].copy_from_slice(&buf[buf_start..buf_stop]);
                }
            }
        }
    }
}

#[inline]
pub(crate) fn copy_f_5d_s8<T>(out: &mut [T], out_offsets: &[usize; 5], out_shape: &[usize; 2], buf: &[T], buf_shape: &[usize; 5])
where
    T: Copy,
{
    let is_diag_ij = out_offsets[0] == out_offsets[1];
    let is_diag_kl = out_offsets[2] == out_offsets[3];
    let is_diag_jl = out_offsets[1] == out_offsets[3];
    let is_i_gt_k = out_offsets[0] > out_offsets[2];
    let is_diag_ik = out_offsets[0] == out_offsets[2];
    for c in 0..buf_shape[4] {
        for l in 0..buf_shape[3] {
            let k_max = if is_diag_kl { l + 1 } else { buf_shape[2] };
            for k in 0..k_max {
                let j_max = if is_diag_jl { l + 1 } else { buf_shape[1] };
                for j in 0..j_max {
                    if is_diag_jl && j == l && is_i_gt_k {
                        continue;
                    }
                    let i_max = if is_diag_ij { j + 1 } else { buf_shape[0] };
                    let i_max = if is_diag_jl && j == l && is_diag_ik { i_max.min(k + 1) } else { i_max };

                    let out_indices = [out_offsets[0], out_offsets[1] + j, out_offsets[2] + k, out_offsets[3] + l, out_offsets[4] + c];
                    let buf_indices = [0, j, k, l, c];
                    let out_start = get_f_index_5d_s8(&out_indices, out_shape);
                    let buf_start = get_f_index_5d(&buf_indices, buf_shape);
                    let out_stop = out_start + i_max;
                    let buf_stop = buf_start + i_max;
                    out[out_start..out_stop].copy_from_slice(&buf[buf_start..buf_stop]);
                }
            }
        }
    }
}

/* #endregion */

/* #region indices computation (row-major) */

#[inline(always)]
pub(crate) fn get_c_index_3d(indices: &[usize; 3], shape: &[usize; 3]) -> usize {
    let [i0, i1, i2] = indices;
    let [_, s1, s2] = shape;
    i2 + s2 * (i1 + s1 * i0)
}

#[inline(always)]
pub(crate) fn get_c_index_4d(indices: &[usize; 4], shape: &[usize; 4]) -> usize {
    let [i0, i1, i2, i3] = indices;
    let [_, s1, s2, s3] = shape;
    i3 + s3 * (i2 + s2 * (i1 + s1 * i0))
}

#[inline(always)]
pub(crate) fn get_c_index_5d(indices: &[usize; 5], shape: &[usize; 5]) -> usize {
    let [i0, i1, i2, i3, i4] = indices;
    let [_, s1, s2, s3, s4] = shape;
    i4 + s4 * (i3 + s3 * (i2 + s2 * (i1 + s1 * i0)))
}

#[inline(always)]
pub(crate) fn get_c_index_3d_s2ij(indices: &[usize; 3], shape: &[usize; 2]) -> usize {
    let [i0, i1, i2] = indices;
    let [_, s1] = shape;
    i2 + i1 * (i1 + 1) / 2 + i0 * s1
}

#[inline(always)]
pub(crate) fn get_c_index_4d_s2ij(indices: &[usize; 4], shape: &[usize; 3]) -> usize {
    let [i0, i1, i2, i3] = indices;
    let [_, s1, s2] = shape;
    i3 + s2 * (i2 + i1 * (i1 + 1) / 2 + i0 * s1)
}

#[inline(always)]
pub(crate) fn get_c_index_5d_s2ij(indices: &[usize; 5], shape: &[usize; 4]) -> usize {
    let [i0, i1, i2, i3, i4] = indices;
    let [_, s1, s2, s3] = shape;
    i4 + s3 * (i3 + s2 * (i2 + i1 * (i1 + 1) / 2 + i0 * s1))
}

/* #endregion */

/* #region integral block copy (row out from col buffer) */

pub(crate) fn copy_c_3d_s1<T>(out: &mut [T], out_offsets: &[usize; 3], out_shape: &[usize; 3], buf: &[T], buf_shape: &[usize; 3])
where
    T: Copy,
{
    let buf_stride = buf_shape[0];
    for c in 0..buf_shape[2] {
        for i in 0..buf_shape[0] {
            let out_indices = [out_offsets[0] + c, out_offsets[1] + i, out_offsets[2]];
            let buf_indices = [i, 0, c];
            let out_idx = get_c_index_3d(&out_indices, out_shape);
            let buf_idx = get_f_index_3d(&buf_indices, buf_shape);
            for j in 0..buf_shape[1] {
                out[out_idx + j] = buf[buf_idx + j * buf_stride];
            }
        }
    }
}

pub(crate) fn copy_c_4d_s1<T>(out: &mut [T], out_offsets: &[usize; 4], out_shape: &[usize; 4], buf: &[T], buf_shape: &[usize; 4])
where
    T: Copy,
{
    let buf_stride = buf_shape[0] * buf_shape[1];
    for c in 0..buf_shape[3] {
        for i in 0..buf_shape[0] {
            for j in 0..buf_shape[1] {
                let out_indices = [out_offsets[0] + c, out_offsets[1] + i, out_offsets[2] + j, out_offsets[3]];
                let buf_indices = [i, j, 0, c];
                let out_idx = get_c_index_4d(&out_indices, out_shape);
                let buf_idx = get_f_index_4d(&buf_indices, buf_shape);
                for k in 0..buf_shape[2] {
                    out[out_idx + k] = buf[buf_idx + k * buf_stride];
                }
            }
        }
    }
}

pub(crate) fn copy_c_5d_s1<T>(out: &mut [T], out_offsets: &[usize; 5], out_shape: &[usize; 5], buf: &[T], buf_shape: &[usize; 5])
where
    T: Copy,
{
    let buf_stride = buf_shape[0] * buf_shape[1] * buf_shape[2];
    for c in 0..buf_shape[4] {
        for i in 0..buf_shape[0] {
            for j in 0..buf_shape[1] {
                for k in 0..buf_shape[2] {
                    let out_indices = [out_offsets[0] + c, out_offsets[1] + i, out_offsets[2] + j, out_offsets[3] + k, out_offsets[4]];
                    let buf_indices = [i, j, k, 0, c];
                    let out_idx = get_c_index_5d(&out_indices, out_shape);
                    let buf_idx = get_f_index_5d(&buf_indices, buf_shape);
                    for l in 0..buf_shape[3] {
                        out[out_idx + l] = buf[buf_idx + l * buf_stride];
                    }
                }
            }
        }
    }
}

#[inline]
pub(crate) fn copy_c_3d_s2ij<T>(out: &mut [T], out_offsets: &[usize; 3], out_shape: &[usize; 2], buf: &[T], buf_shape: &[usize; 3])
where
    T: Copy,
{
    let buf_stride = buf_shape[0];
    let is_diag_ij = out_offsets[1] == out_offsets[2];
    for c in 0..buf_shape[2] {
        for i in 0..buf_shape[0] {
            let out_indices = [out_offsets[0] + c, out_offsets[1] + i, out_offsets[2]];
            let buf_indices = [i, 0, c];
            let out_idx = get_c_index_3d_s2ij(&out_indices, out_shape);
            let buf_idx = get_f_index_3d(&buf_indices, buf_shape);
            let j_max = if is_diag_ij { i + 1 } else { buf_shape[1] };
            for j in 0..j_max {
                out[out_idx + j] = buf[buf_idx + j * buf_stride];
            }
        }
    }
}

#[inline]
pub(crate) fn copy_c_4d_s2ij<T>(out: &mut [T], out_offsets: &[usize; 4], out_shape: &[usize; 3], buf: &[T], buf_shape: &[usize; 4])
where
    T: Copy,
{
    let buf_stride = buf_shape[0] * buf_shape[1];
    let is_diag_ij = out_offsets[1] == out_offsets[2];
    for c in 0..buf_shape[3] {
        for i in 0..buf_shape[0] {
            let j_max = if is_diag_ij { i + 1 } else { buf_shape[1] };
            for j in 0..j_max {
                let out_indices = [out_offsets[0] + c, out_offsets[1] + i, out_offsets[2] + j, out_offsets[3]];
                let buf_indices = [i, j, 0, c];
                let out_idx = get_c_index_4d_s2ij(&out_indices, out_shape);
                let buf_idx = get_f_index_4d(&buf_indices, buf_shape);
                for k in 0..buf_shape[2] {
                    out[out_idx + k] = buf[buf_idx + k * buf_stride];
                }
            }
        }
    }
}

#[inline]
pub(crate) fn copy_c_5d_s2ij<T>(out: &mut [T], out_offsets: &[usize; 5], out_shape: &[usize; 4], buf: &[T], buf_shape: &[usize; 5])
where
    T: Copy,
{
    let buf_stride = buf_shape[0] * buf_shape[1] * buf_shape[2];
    let is_diag_ij = out_offsets[1] == out_offsets[2];
    for c in 0..buf_shape[4] {
        for i in 0..buf_shape[0] {
            let j_max = if is_diag_ij { i + 1 } else { buf_shape[1] };
            for j in 0..j_max {
                for k in 0..buf_shape[2] {
                    let out_indices = [out_offsets[0] + c, out_offsets[1] + i, out_offsets[2] + j, out_offsets[3] + k, out_offsets[4]];
                    let buf_indices = [i, j, k, 0, c];
                    let out_idx = get_c_index_5d_s2ij(&out_indices, out_shape);
                    let buf_idx = get_f_index_5d(&buf_indices, buf_shape);
                    for l in 0..buf_shape[3] {
                        out[out_idx + l] = buf[buf_idx + l * buf_stride];
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
            v
        } else {
            panic!("Allocation failed (probably due to out-of-memory, of size {size})")
        }
    }
}

/* #endregion */
