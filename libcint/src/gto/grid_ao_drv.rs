use crate::gto::prelude_dev::*;

pub fn gto_nabla1(
    fx1: &mut [f64blk],
    fy1: &mut [f64blk],
    fz1: &mut [f64blk],
    fx0: &[f64blk],
    fy0: &[f64blk],
    fz0: &[f64blk],
    l: usize,
    a: f64,
) {
    let a2 = -2.0 * a;
    // first derivative
    for n in 0..BLKSIZE {
        fx1[0][n] = a2 * fx0[1][n];
        fy1[0][n] = a2 * fy0[1][n];
        fz1[0][n] = a2 * fz0[1][n];
    }
    // recursive derivatives
    for i in 1..=l {
        for n in 0..BLKSIZE {
            fx1[i][n] = i as f64 * fx0[i - 1][n] + a2 * fx0[i + 1][n];
            fy1[i][n] = i as f64 * fy0[i - 1][n] + a2 * fy0[i + 1][n];
            fz1[i][n] = i as f64 * fz0[i - 1][n] + a2 * fz0[i + 1][n];
        }
    }
}

pub fn gto_x1(
    fx1: &mut [f64blk],
    fy1: &mut [f64blk],
    fz1: &mut [f64blk],
    fx0: &[f64blk],
    fy0: &[f64blk],
    fz0: &[f64blk],
    l: usize,
    ri: [f64; 3],
) {
    for i in 0..=l {
        for n in 0..BLKSIZE {
            fx1[i][n] = ri[0] * fx0[i][n] + fx0[i + 1][n];
            fy1[i][n] = ri[1] * fy0[i][n] + fy0[i + 1][n];
            fz1[i][n] = ri[2] * fz0[i][n] + fz0[i + 1][n];
        }
    }
}

pub fn gto_prim_exp(eprim: &mut [f64blk], coord: &[f64blk; 3], alpha: &[f64], _coeff: &[f64], _nprim: usize, nctr: usize, fac: f64) {
    let mut rr = unsafe { f64blk::uninit() };
    for i in 0..BLKSIZE {
        rr[i] = coord[0][i] * coord[0][i] + coord[1][i] * coord[1][i] + coord[2][i] * coord[2][i];
    }

    for j in 0..nctr {
        for i in 0..BLKSIZE {
            let arr = alpha[j] * rr[i];
            eprim[j][i] = fac * (-arr).exp();
        }
    }
}
