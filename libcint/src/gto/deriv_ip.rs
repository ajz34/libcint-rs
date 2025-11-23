use crate::gto::prelude_dev::*;

pub fn gto_contract_exp1(
    // arguments
    ectr: &mut [f64blk],
    coord: &[f64blk; 3],
    alpha: &[f64],
    coeff: &[f64],
    fac: f64,
    // dimensions
    nprim: usize,
    nctr: usize,
) {
    const X: usize = 0;
    const Y: usize = 1;
    const Z: usize = 2;

    let (ectr, ectr_2a) = ectr.split_at_mut(nctr);

    // r² = x² + y² + z²
    let mut rr = unsafe { f64blk::uninit() };
    for g in 0..BLKSIZE {
        rr[g] = coord[X][g].mul_add(coord[X][g], coord[Y][g].mul_add(coord[Y][g], coord[Z][g] * coord[Z][g]));
    }
    // zero ectr
    for k in 0..nctr {
        for g in 0..BLKSIZE {
            ectr[k][g] = 0.0;
            ectr_2a[k][g] = 0.0;
        }
    }
    // -2 alpha_i C_ij exp(-alpha_i r_k^2)
    for p in 0..nprim {
        let mut eprim = unsafe { f64blk::uninit() };
        for g in 0..BLKSIZE {
            let arr = alpha[p] * rr[g];
            eprim[g] = (-arr).exp() * fac;
        }
        for k in 0..nctr {
            let coeff_2a = -2.0 * alpha[p] * coeff[k * nprim + p];
            for g in 0..BLKSIZE {
                ectr[k][g] = eprim[g].mul_add(coeff[k * nprim + p], ectr[k][g]);
                ectr_2a[k][g] = eprim[g].mul_add(coeff_2a, ectr_2a[k][g]);
            }
        }
    }
}

pub fn gto_shell_eval_grid_cart_deriv1(
    // arguments
    gto: &mut [f64blk],
    exps: &[f64blk],
    coord: &[f64blk; 3],
    l: usize,
    // dimensions
    nctr: usize,
) {
    // const ANG_MAX: usize = crate::ffi::cint_ffi::ANG_MAX as usize;

    const D1_0: usize = 0;
    const D1_X: usize = 1;
    const D1_Y: usize = 2;
    const D1_Z: usize = 3;

    let (exps, exp_2a) = exps.split_at(nctr);
    let mut gto = gto.chunks_exact_mut(4 * nctr).collect_vec();

    match l {
        0 => {
            for k in 0..nctr {
                for g in 0..BLKSIMD {
                    let e = exps[k].get_simd(g);
                    let e_2a = exp_2a[k].get_simd(g);
                    let x = coord[X].get_simd(g);
                    let y = coord[Y].get_simd(g);
                    let z = coord[Z].get_simd(g);

                    *gto[D1_0][k].get_simd_mut(g) = e;
                    *gto[D1_X][k].get_simd_mut(g) = e_2a * x;
                    *gto[D1_Y][k].get_simd_mut(g) = e_2a * y;
                    *gto[D1_Z][k].get_simd_mut(g) = e_2a * z;
                }
            }
        },
        1 => {
            for k in 0..nctr {
                for g in 0..BLKSIMD {
                    let e = exps[k].get_simd(g);
                    let e_2a = exp_2a[k].get_simd(g);
                    let x = coord[X].get_simd(g);
                    let y = coord[Y].get_simd(g);
                    let z = coord[Z].get_simd(g);

                    *gto[D1_0][3 * k + X].get_simd_mut(g) = e * x;
                    *gto[D1_0][3 * k + Y].get_simd_mut(g) = e * y;
                    *gto[D1_0][3 * k + Z].get_simd_mut(g) = e * z;

                    *gto[D1_X][3 * k + X].get_simd_mut(g) = e_2a * x * x + e_2a;
                    *gto[D1_X][3 * k + Y].get_simd_mut(g) = e_2a * x * y;
                    *gto[D1_X][3 * k + Z].get_simd_mut(g) = e_2a * x * z;

                    *gto[D1_Y][3 * k + X].get_simd_mut(g) = e_2a * y * x;
                    *gto[D1_Y][3 * k + Y].get_simd_mut(g) = e_2a * y * y + e_2a;
                    *gto[D1_Y][3 * k + Z].get_simd_mut(g) = e_2a * y * z;

                    *gto[D1_Z][3 * k + X].get_simd_mut(g) = e_2a * z * x;
                    *gto[D1_Z][3 * k + Y].get_simd_mut(g) = e_2a * z * y;
                    *gto[D1_Z][3 * k + Z].get_simd_mut(g) = e_2a * z * z + e_2a;
                }
            }
        },
        2 => {
            for k in 0..nctr {
                for g in 0..BLKSIMD {
                    let e = exps[k].get_simd(g);
                    let e_2a = exp_2a[k].get_simd(g);
                    let x = coord[X].get_simd(g);
                    let y = coord[Y].get_simd(g);
                    let z = coord[Z].get_simd(g);

                    let ex = e * x;
                    let ey = e * y;
                    let ez = e * z;
                    let xx = x * x;
                    let xy = x * y;
                    let xz = x * z;
                    let yy = y * y;
                    let yz = y * z;
                    let zz = z * z;
                    let ax = e_2a * x;
                    let ay = e_2a * y;
                    let az = e_2a * z;

                    *gto[D1_0][6 * k + XX].get_simd_mut(g) = e * xx;
                    *gto[D1_0][6 * k + XY].get_simd_mut(g) = e * xy;
                    *gto[D1_0][6 * k + XZ].get_simd_mut(g) = e * xz;
                    *gto[D1_0][6 * k + YY].get_simd_mut(g) = e * yy;
                    *gto[D1_0][6 * k + YZ].get_simd_mut(g) = e * yz;
                    *gto[D1_0][6 * k + ZZ].get_simd_mut(g) = e * zz;

                    *gto[D1_X][6 * k + XX].get_simd_mut(g) = ax * xx + ex * 2.0;
                    *gto[D1_X][6 * k + XY].get_simd_mut(g) = ax * xy + ey;
                    *gto[D1_X][6 * k + XZ].get_simd_mut(g) = ax * xz + ez;
                    *gto[D1_X][6 * k + YY].get_simd_mut(g) = ax * yy;
                    *gto[D1_X][6 * k + YZ].get_simd_mut(g) = ax * yz;
                    *gto[D1_X][6 * k + ZZ].get_simd_mut(g) = ax * zz;

                    *gto[D1_Y][6 * k + XX].get_simd_mut(g) = ay * xx;
                    *gto[D1_Y][6 * k + XY].get_simd_mut(g) = ay * xy + ex;
                    *gto[D1_Y][6 * k + XZ].get_simd_mut(g) = ay * xz;
                    *gto[D1_Y][6 * k + YY].get_simd_mut(g) = ay * yy + ey * 2.0;
                    *gto[D1_Y][6 * k + YZ].get_simd_mut(g) = ay * yz + ez;
                    *gto[D1_Y][6 * k + ZZ].get_simd_mut(g) = ay * zz;

                    *gto[D1_Z][6 * k + XX].get_simd_mut(g) = az * xx;
                    *gto[D1_Z][6 * k + XY].get_simd_mut(g) = az * xy;
                    *gto[D1_Z][6 * k + XZ].get_simd_mut(g) = az * xz + ex;
                    *gto[D1_Z][6 * k + YY].get_simd_mut(g) = az * yy;
                    *gto[D1_Z][6 * k + YZ].get_simd_mut(g) = az * yz + ey;
                    *gto[D1_Z][6 * k + ZZ].get_simd_mut(g) = az * zz + ez * 2.0;
                }
            }
        },
        3 => {
            for k in 0..nctr {
                for g in 0..BLKSIMD {
                    let e = exps[k].get_simd(g);
                    let e_2a = exp_2a[k].get_simd(g);
                    let x = coord[X].get_simd(g);
                    let y = coord[Y].get_simd(g);
                    let z = coord[Z].get_simd(g);

                    let exx = e * x * x;
                    let exy = e * x * y;
                    let exz = e * x * z;
                    let eyy = e * y * y;
                    let eyz = e * y * z;
                    let ezz = e * z * z;

                    let xxx = x * x * x;
                    let xxy = x * x * y;
                    let xxz = x * x * z;
                    let xyy = x * y * y;
                    let xyz = x * y * z;
                    let xzz = x * z * z;
                    let yyy = y * y * y;
                    let yyz = y * y * z;
                    let yzz = y * z * z;
                    let zzz = z * z * z;

                    let ax = e_2a * x;
                    let ay = e_2a * y;
                    let az = e_2a * z;

                    *gto[D1_0][10 * k + XXX].get_simd_mut(g) = e * xxx;
                    *gto[D1_0][10 * k + XXY].get_simd_mut(g) = e * xxy;
                    *gto[D1_0][10 * k + XXZ].get_simd_mut(g) = e * xxz;
                    *gto[D1_0][10 * k + XYY].get_simd_mut(g) = e * xyy;
                    *gto[D1_0][10 * k + XYZ].get_simd_mut(g) = e * xyz;
                    *gto[D1_0][10 * k + XZZ].get_simd_mut(g) = e * xzz;
                    *gto[D1_0][10 * k + YYY].get_simd_mut(g) = e * yyy;
                    *gto[D1_0][10 * k + YYZ].get_simd_mut(g) = e * yyz;
                    *gto[D1_0][10 * k + YZZ].get_simd_mut(g) = e * yzz;
                    *gto[D1_0][10 * k + ZZZ].get_simd_mut(g) = e * zzz;

                    *gto[D1_X][10 * k + XXX].get_simd_mut(g) = ax * xxx + exx * 3.0;
                    *gto[D1_X][10 * k + XXY].get_simd_mut(g) = ax * xxy + exy * 2.0;
                    *gto[D1_X][10 * k + XXZ].get_simd_mut(g) = ax * xxz + exz * 2.0;
                    *gto[D1_X][10 * k + XYY].get_simd_mut(g) = ax * xyy + eyy;
                    *gto[D1_X][10 * k + XYZ].get_simd_mut(g) = ax * xyz + eyz;
                    *gto[D1_X][10 * k + XZZ].get_simd_mut(g) = ax * xzz + ezz;
                    *gto[D1_X][10 * k + YYY].get_simd_mut(g) = ax * yyy;
                    *gto[D1_X][10 * k + YYZ].get_simd_mut(g) = ax * yyz;
                    *gto[D1_X][10 * k + YZZ].get_simd_mut(g) = ax * yzz;
                    *gto[D1_X][10 * k + ZZZ].get_simd_mut(g) = ax * zzz;

                    *gto[D1_Y][10 * k + XXX].get_simd_mut(g) = ay * xxx;
                    *gto[D1_Y][10 * k + XXY].get_simd_mut(g) = ay * xxy + exx;
                    *gto[D1_Y][10 * k + XXZ].get_simd_mut(g) = ay * xxz;
                    *gto[D1_Y][10 * k + XYY].get_simd_mut(g) = ay * xyy + exy * 2.0;
                    *gto[D1_Y][10 * k + XYZ].get_simd_mut(g) = ay * xyz + exz;
                    *gto[D1_Y][10 * k + XZZ].get_simd_mut(g) = ay * xzz;
                    *gto[D1_Y][10 * k + YYY].get_simd_mut(g) = ay * yyy + eyy * 3.0;
                    *gto[D1_Y][10 * k + YYZ].get_simd_mut(g) = ay * yyz + eyz * 2.0;
                    *gto[D1_Y][10 * k + YZZ].get_simd_mut(g) = ay * yzz + ezz;
                    *gto[D1_Y][10 * k + ZZZ].get_simd_mut(g) = ay * zzz;

                    *gto[D1_Z][10 * k + XXX].get_simd_mut(g) = az * xxx;
                    *gto[D1_Z][10 * k + XXY].get_simd_mut(g) = az * xxy;
                    *gto[D1_Z][10 * k + XXZ].get_simd_mut(g) = az * xxz + exx;
                    *gto[D1_Z][10 * k + XYY].get_simd_mut(g) = az * xyy;
                    *gto[D1_Z][10 * k + XYZ].get_simd_mut(g) = az * xyz + exy;
                    *gto[D1_Z][10 * k + XZZ].get_simd_mut(g) = az * xzz + exz * 2.0;
                    *gto[D1_Z][10 * k + YYY].get_simd_mut(g) = az * yyy;
                    *gto[D1_Z][10 * k + YYZ].get_simd_mut(g) = az * yyz + eyy;
                    *gto[D1_Z][10 * k + YZZ].get_simd_mut(g) = az * yzz + eyz * 2.0;
                    *gto[D1_Z][10 * k + ZZZ].get_simd_mut(g) = az * zzz + ezz * 3.0;
                }
            }
        },
        _ => panic!("gto_shell_eval_grid_cart_deriv1: l > 3 not implemented"),
    }
}
