use crate::gto::prelude_dev::*;

pub fn gto_shell_eval_grid_cart_ipsp(
    // arguments
    gto: &mut [f64blk],
    eprim: &[f64blk],
    coord: &[f64blk; 3],
    alpha: &[f64],
    coeff: &[f64],
    l: usize,
    // dimensions
    shl_shape: [usize; 2],
) {
    const ANG_MAX: usize = crate::ffi::cint_ffi::ANG_MAX as usize;

    const COMP_NUM: usize = 12;

    let [nctr, nprim] = shl_shape;
    let ncart = (l + 1) * (l + 2) / 2;
    let nao_to_set = nctr * ncart;
    let mut f0 = [[f64simd::zero(); 3]; ANG_MAX + 3];
    let mut f1 = [[f64simd::zero(); 3]; ANG_MAX + 3];
    let mut f2 = [[f64simd::zero(); 3]; ANG_MAX + 3];
    let mut f3 = [[f64simd::zero(); 3]; ANG_MAX + 3];
    let mut s = [[f64simd::zero(); 3]; 3];
    let mut buf = [f64simd::zero(); COMP_NUM];
    let mut gto = gto.chunks_exact_mut(nao_to_set).collect_vec();

    // zero out the output buffer
    for icomp in 0..COMP_NUM {
        for mu in 0..nao_to_set {
            for g in 0..BLKSIMDD {
                gto[icomp][mu].get_simdd_mut(g).fill(0.0);
            }
        }
    }

    for g in 0..BLKSIMDD {
        let x = coord[0].get_simdd(g);
        let y = coord[1].get_simdd(g);
        let z = coord[2].get_simdd(g);

        for p in 0..nprim {
            let e = eprim[p].get_simdd(g);
            if e.is_gto_zero() {
                continue;
            }
            f0[0][X] = f64simd::splat(1.0);
            f0[0][Y] = f64simd::splat(1.0);
            f0[0][Z] = f64simd::splat(1.0);
            for ll in 1..=l + 2 {
                f0[ll][X] = f0[ll - 1][X] * x;
                f0[ll][Y] = f0[ll - 1][Y] * y;
                f0[ll][Z] = f0[ll - 1][Z] * z;
            }
            gto_nabla1_simdd(&mut f1, &f0, l, alpha[p]);
            gto_nabla1_simdd(&mut f2, &f0, l + 1, alpha[p]);
            gto_nabla1_simdd(&mut f3, &f2, l, alpha[p]);
            for (icart, (lx, ly, lz)) in gto_l_iter(l).enumerate() {
                s[X][X] = e * f3[lx][X] * f0[ly][Y] * f0[lz][Z];
                s[X][Y] = e * f2[lx][X] * f1[ly][Y] * f0[lz][Z];
                s[X][Z] = e * f2[lx][X] * f0[ly][Y] * f1[lz][Z];
                s[Y][X] = e * f1[lx][X] * f2[ly][Y] * f0[lz][Z];
                s[Y][Y] = e * f0[lx][X] * f3[ly][Y] * f0[lz][Z];
                s[Y][Z] = e * f0[lx][X] * f2[ly][Y] * f1[lz][Z];
                s[Z][X] = e * f1[lx][X] * f0[ly][Y] * f2[lz][Z];
                s[Z][Y] = e * f0[lx][X] * f1[ly][Y] * f2[lz][Z];
                s[Z][Z] = e * f0[lx][X] * f0[ly][Y] * f3[lz][Z];

                buf[0] = -s[X][X];
                buf[1] = -s[X][Y];
                buf[2] = -s[X][Z];
                buf[3] = f64simd::zero();
                buf[4] = -s[Y][X];
                buf[5] = -s[Y][Y];
                buf[6] = -s[Y][Z];
                buf[7] = f64simd::zero();
                buf[8] = -s[Z][X];
                buf[9] = -s[Z][Y];
                buf[10] = -s[Z][Z];
                buf[11] = f64simd::zero();

                for k in 0..nctr {
                    let c = f64simd::splat(coeff[k * nprim + p]);
                    for icomp in 0..COMP_NUM {
                        gto[icomp][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[icomp]);
                    }
                }
            }
        }
    }
}

pub struct GtoEvalDerivIpSp;
impl GtoEvalAPI for GtoEvalDerivIpSp {
    fn ne1(&self) -> usize {
        4
    }
    fn ntensor(&self) -> usize {
        3
    }
    fn gto_exp(&self, ebuf: &mut [f64blk], coord: &[f64blk; 3], alpha: &[f64], _coeff: &[f64], fac: f64, shl_shape: [usize; 2]) {
        let eprim = ebuf;
        let [_nctr, nprim] = shl_shape;
        gto_prim_exp(eprim, coord, alpha, fac, nprim);
    }
    fn gto_shell_eval(
        &self,
        // arguments
        gto: &mut [f64blk],
        ebuf: &[f64blk],
        coord: &[f64blk; 3],
        alpha: &[f64],
        coeff: &[f64],
        l: usize,
        _center: [f64; 3],
        // dimensions
        shl_shape: [usize; 2],
    ) {
        gto_shell_eval_grid_cart_ipsp(gto, ebuf, coord, alpha, coeff, l, shl_shape);
    }
    unsafe fn cint_c2s_ket_spinor(
        &self,
        gspa: *mut Complex<f64>,
        gspb: *mut Complex<f64>,
        gcart: *const f64,
        lds: c_int,
        ldc: c_int,
        nctr: c_int,
        kappa: c_int,
        l: c_int,
    ) {
        crate::ffi::cint_ffi::CINTc2s_ket_spinor_si1(gspa as _, gspb as _, gcart as _, lds, ldc, nctr, kappa, l);
    }
}

#[test]
fn test_usual_case() {
    // this will test usual case (naive `eval_gto` call)
    let mol = init_c10h22_def2_qzvp();
    let ngrid = 2048;
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    let (out, shape) = mol.eval_gto("GTOval_ipsp_sph", &coord).into();
    let out_fp = cint_fp(&out);
    println!("out_fp = {}", out_fp);
    let nao = shape[1];
    let out_fp0 = cint_fp(&out[0..3 * nao * ngrid]);
    println!("out_fp0 = {}", out_fp0);
    assert!((out_fp0 - 3846.4314775880575).abs() < 1e-8);
}
