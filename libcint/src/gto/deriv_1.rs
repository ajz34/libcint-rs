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
    const ANG_MAX: usize = crate::ffi::cint_ffi::ANG_MAX as usize;

    const D1_0: usize = 0;
    const D1_X: usize = 1;
    const D1_Y: usize = 2;
    const D1_Z: usize = 3;

    let ncart = (l + 1) * (l + 2) / 2;
    let (exps, exps_2a) = exps.split_at(nctr);
    let mut gto = gto.chunks_exact_mut(nctr * ncart).collect_vec();

    match l {
        0 => {
            for k in 0..nctr {
                for g in 0..BLKSIMD {
                    let e = exps[k].get_simd(g);
                    let e_2a = exps_2a[k].get_simd(g);
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
                    let e_2a = exps_2a[k].get_simd(g);
                    let x = coord[X].get_simd(g);
                    let y = coord[Y].get_simd(g);
                    let z = coord[Z].get_simd(g);

                    *gto[D1_0][3 * k + X].get_simd_mut(g) = e * x;
                    *gto[D1_0][3 * k + Y].get_simd_mut(g) = e * y;
                    *gto[D1_0][3 * k + Z].get_simd_mut(g) = e * z;

                    *gto[D1_X][3 * k + X].get_simd_mut(g) = e_2a * x * x + e;
                    *gto[D1_X][3 * k + Y].get_simd_mut(g) = e_2a * x * y;
                    *gto[D1_X][3 * k + Z].get_simd_mut(g) = e_2a * x * z;

                    *gto[D1_Y][3 * k + X].get_simd_mut(g) = e_2a * y * x;
                    *gto[D1_Y][3 * k + Y].get_simd_mut(g) = e_2a * y * y + e;
                    *gto[D1_Y][3 * k + Z].get_simd_mut(g) = e_2a * y * z;

                    *gto[D1_Z][3 * k + X].get_simd_mut(g) = e_2a * z * x;
                    *gto[D1_Z][3 * k + Y].get_simd_mut(g) = e_2a * z * y;
                    *gto[D1_Z][3 * k + Z].get_simd_mut(g) = e_2a * z * z + e;
                }
            }
        },
        2 => {
            for k in 0..nctr {
                for g in 0..BLKSIMD {
                    let e = exps[k].get_simd(g);
                    let e_2a = exps_2a[k].get_simd(g);
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
                    let e_2a = exps_2a[k].get_simd(g);
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
        _ => {
            // initialize pows buffer
            let mut pows_buffer = unsafe { [[f64blk::uninit(); 3]; ANG_MAX + 2] };
            for g in 0..BLKSIMD {
                pows_buffer[0][X].get_simd_mut(g).fill(0.0);
                pows_buffer[0][Y].get_simd_mut(g).fill(0.0);
                pows_buffer[0][Z].get_simd_mut(g).fill(0.0);
                pows_buffer[1][X].get_simd_mut(g).fill(1.0);
                pows_buffer[1][Y].get_simd_mut(g).fill(1.0);
                pows_buffer[1][Z].get_simd_mut(g).fill(1.0);
            }
            let pows = &mut pows_buffer[1..];
            for lx in 0..l {
                for g in 0..BLKSIMD {
                    *pows[lx + 1][X].get_simd_mut(g) = pows[lx][X].get_simd(g) * coord[X].get_simd(g);
                    *pows[lx + 1][Y].get_simd_mut(g) = pows[lx][Y].get_simd(g) * coord[Y].get_simd(g);
                    *pows[lx + 1][Z].get_simd_mut(g) = pows[lx][Z].get_simd(g) * coord[Z].get_simd(g);
                }
            }

            let pows = &pows_buffer[1..];
            let pows_1less = &pows_buffer;

            for k in 0..nctr {
                let mut icart = 0;
                for lx in (0..=l).rev() {
                    for ly in (0..=(l - lx)).rev() {
                        let lz = l - lx - ly;
                        for g in 0..BLKSIMD {
                            let e = exps[k].get_simd(g);
                            let e_2a = exps_2a[k].get_simd(g);
                            let x = coord[X].get_simd(g);
                            let y = coord[Y].get_simd(g);
                            let z = coord[Z].get_simd(g);
                            let pows_x = pows[lx][X].get_simd(g);
                            let pows_y = pows[ly][Y].get_simd(g);
                            let pows_z = pows[lz][Z].get_simd(g);
                            let pows_1less_x = pows_1less[lx][X].get_simd(g);
                            let pows_1less_y = pows_1less[ly][Y].get_simd(g);
                            let pows_1less_z = pows_1less[lz][Z].get_simd(g);
                            let pows_xyz = pows_x * pows_y * pows_z;

                            let offset = ncart * k + icart;
                            *gto[D1_0][offset].get_simd_mut(g) = e * pows_xyz;
                            *gto[D1_X][offset].get_simd_mut(g) = e_2a * pows_xyz * x + e * (lx as f64) * pows_1less_x * pows_y * pows_z;
                            *gto[D1_Y][offset].get_simd_mut(g) = e_2a * pows_xyz * y + e * (ly as f64) * pows_x * pows_1less_y * pows_z;
                            *gto[D1_Z][offset].get_simd_mut(g) = e_2a * pows_xyz * z + e * (lz as f64) * pows_x * pows_y * pows_1less_z;
                        }
                        icart += 1;
                    }
                }
            }
        },
    }
}

pub struct GTOEvalDeriv1 {}
impl GTOEvalAPI for GTOEvalDeriv1 {
    fn ne1(&self) -> usize {
        1
    }
    fn ntensor(&self) -> usize {
        4
    }
    fn gto_exp(
        &self,
        // arguments
        ebuf: &mut [f64blk],
        coord: &[f64blk; 3],
        alpha: &[f64],
        coeff: &[f64],
        fac: f64,
        // dimensions
        nprim: usize,
        nctr: usize,
    ) {
        let ectr = ebuf;
        gto_contract_exp1(ectr, coord, alpha, coeff, fac, nprim, nctr);
    }
    fn gto_shell_eval(
        &self,
        // arguments
        gto: &mut [f64blk],
        exps: &[f64blk],
        coord: &[f64blk; 3],
        l: usize,
        // dimensions
        nctr: usize,
    ) {
        gto_shell_eval_grid_cart_deriv1(gto, exps, coord, l, nctr);
    }
}

#[test]
fn playground_cart() {
    let mut cint_data = init_h2o_def2_tzvp();
    cint_data.set_cint_type("cart");
    let ngrid = 2048;
    let nao = cint_data.nao();
    let mut ao = vec![0.0f64; 4 * nao * ngrid];
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();

    gto_eval_loop::<true>(
        &GTOEvalDeriv1 {},
        &mut ao,
        &coord,
        1.0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        None,
        true,
        &cint_data.atm,
        &cint_data.bas,
        &cint_data.env,
    );

    println!("fp(ao): {:}", cint_fp(&ao));
}

#[test]
fn test_c10h22_sph() {
    let cint_data = init_c10h22_def2_qzvp();
    let ngrid = 2048;
    let nao = cint_data.nao();
    println!("nao: {:}", nao);
    let mut ao = vec![0.0f64; 4 * nao * ngrid];
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    gto_eval_loop::<false>(
        &GTOEvalDeriv1 {},
        &mut ao,
        &coord,
        1.0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        None,
        true,
        &cint_data.atm,
        &cint_data.bas,
        &cint_data.env,
    );
    let fp_ao = cint_fp(&ao);
    println!("fp_ao = {:}", fp_ao);
    assert!((fp_ao - 1158.9931602933664).abs() < 1e-10);

    let ngrid = 262144;
    let mut ao = vec![0.0f64; 4 * nao * ngrid];
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    let time = std::time::Instant::now();
    gto_eval_loop::<false>(
        &GTOEvalDeriv1 {},
        &mut ao,
        &coord,
        1.0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        None,
        true,
        &cint_data.atm,
        &cint_data.bas,
        &cint_data.env,
    );
    let elapsed = time.elapsed();
    println!("time: {:.3} s", elapsed.as_secs_f64());
    let fp_ao = cint_fp(&ao);
    println!("fp_ao (1M grids) = {:}", fp_ao);
}

#[test]
fn test_c10h22_sph_non0tab() {
    let cint_data = init_c10h22_def2_qzvp();
    let ngrid = 2048;
    let nao = cint_data.nao();
    println!("nao: {:}", nao);
    let mut ao = vec![f64::NAN; 4 * nao * ngrid];
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    let (non0tab, _) = gto_screen_index(&coord, [0, cint_data.nbas()], None, Some(0.01), &cint_data.atm, &cint_data.bas, &cint_data.env);
    gto_eval_loop::<false>(
        &GTOEvalDeriv1 {},
        &mut ao,
        &coord,
        1.0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        Some(&non0tab),
        true,
        &cint_data.atm,
        &cint_data.bas,
        &cint_data.env,
    );
    let fp_ao = cint_fp(&ao);
    println!("fp_ao = {:}", fp_ao);
    assert!((fp_ao - 1150.7578110811928).abs() < 1e-10);
}

#[test]
fn test_c10h22_sph_non0tab_timing() {
    let time = std::time::Instant::now();
    let cint_data = init_c10h22_def2_qzvp();
    let nao = cint_data.nao();
    let ngrid = 262144;
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    let mut ao = vec![0.0; 4 * nao * ngrid];
    let (non0tab, _) = gto_screen_index(&coord, [0, cint_data.nbas()], None, Some(0.01), &cint_data.atm, &cint_data.bas, &cint_data.env);
    gto_eval_loop::<false>(
        &GTOEvalDeriv1 {},
        &mut ao,
        &coord,
        1.0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        Some(&non0tab),
        false,
        &cint_data.atm,
        &cint_data.bas,
        &cint_data.env,
    );
    let elapsed = time.elapsed();
    println!("time: {:.3} s", elapsed.as_secs_f64());
}
