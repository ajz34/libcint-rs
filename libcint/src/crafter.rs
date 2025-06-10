#![allow(dead_code)]

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

/* #endregion */
