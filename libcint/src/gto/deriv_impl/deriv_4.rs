use crate::gto::prelude_dev::*;

pub fn gto_shell_eval_grid_cart_deriv4(
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

    const COMP_NUM: usize = 35;
    const COMP: [[usize; 3]; COMP_NUM] = [
        [0, 0, 0],
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [2, 0, 0],
        [1, 1, 0],
        [1, 0, 1],
        [0, 2, 0],
        [0, 1, 1],
        [0, 0, 2],
        [3, 0, 0],
        [2, 1, 0],
        [2, 0, 1],
        [1, 2, 0],
        [1, 1, 1],
        [1, 0, 2],
        [0, 3, 0],
        [0, 2, 1],
        [0, 1, 2],
        [0, 0, 3],
        [4, 0, 0],
        [3, 1, 0],
        [3, 0, 1],
        [2, 2, 0],
        [2, 1, 1],
        [2, 0, 2],
        [1, 3, 0],
        [1, 2, 1],
        [1, 1, 2],
        [1, 0, 3],
        [0, 4, 0],
        [0, 3, 1],
        [0, 2, 2],
        [0, 1, 3],
        [0, 0, 4],
    ];

    let [nctr, nprim] = shl_shape;
    let ncart = (l + 1) * (l + 2) / 2;
    let nao_to_set = nctr * ncart;
    let mut f0 = [[f64simd::zero(); 3]; ANG_MAX + 4];
    let mut f1 = [[f64simd::zero(); 3]; ANG_MAX + 4];
    let mut f2 = [[f64simd::zero(); 3]; ANG_MAX + 4];
    let mut f3 = [[f64simd::zero(); 3]; ANG_MAX + 4];
    let mut f4 = [[f64simd::zero(); 3]; ANG_MAX + 4];
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
            for ll in 1..=l + 4 {
                f0[ll][X] = f0[ll - 1][X] * x;
                f0[ll][Y] = f0[ll - 1][Y] * y;
                f0[ll][Z] = f0[ll - 1][Z] * z;
            }
            gto_nabla1_simdd(&mut f1, &f0, l + 3, alpha[p]);
            gto_nabla1_simdd(&mut f2, &f1, l + 2, alpha[p]);
            gto_nabla1_simdd(&mut f3, &f2, l + 1, alpha[p]);
            gto_nabla1_simdd(&mut f4, &f3, l, alpha[p]);

            let f = [&f0, &f1, &f2, &f3, &f4];
            for (icart, (lx, ly, lz)) in gto_l_iter(l).enumerate() {
                for (icomp, &[ix, iy, iz]) in COMP.iter().enumerate() {
                    buf[icomp] = e * f[ix][lx][X] * f[iy][ly][Y] * f[iz][lz][Z];
                }
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

pub struct GtoEvalDeriv4;
impl GtoEvalAPI for GtoEvalDeriv4 {
    fn ne1(&self) -> usize {
        1
    }
    fn ntensor(&self) -> usize {
        35
    }
    fn gto_exp(
        &self,
        // arguments
        ebuf: &mut [f64blk],
        coord: &[f64blk; 3],
        alpha: &[f64],
        _coeff: &[f64],
        fac: f64,
        // dimensions
        shl_shape: [usize; 2],
    ) {
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
        gto_shell_eval_grid_cart_deriv4(gto, ebuf, coord, alpha, coeff, l, shl_shape);
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
        crate::ffi::cint_ffi::CINTc2s_ket_spinor_sf1(gspa as _, gspb as _, gcart as _, lds, ldc, nctr, kappa, l);
    }
}
