use crate::gto::prelude_dev::*;

#[allow(non_upper_case_globals)]
pub fn gto_shell_eval_grid_cart_deriv2(
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

    const COMP_NUM: usize = 10;
    const COMP: [[usize; 3]; COMP_NUM] =
        [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1], [2, 0, 0], [1, 1, 0], [1, 0, 1], [0, 2, 0], [0, 1, 1], [0, 0, 2]];

    let [nctr, nprim] = shl_shape;
    let ncart = (l + 1) * (l + 2) / 2;
    let nao_to_set = nctr * ncart;
    let mut f0 = [[f64simd::zero(); 3]; ANG_MAX + 2];
    let mut f1 = [[f64simd::zero(); 3]; ANG_MAX + 2];
    let mut f2 = [[f64simd::zero(); 3]; ANG_MAX + 2];
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
            gto_nabla1_simd(&mut f1, &f0, l + 1, alpha[p]);
            gto_nabla1_simd(&mut f2, &f1, l, alpha[p]);

            let f = [&f0, &f1, &f2];
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

pub struct GtoEvalDeriv2;
impl GtoEvalAPI for GtoEvalDeriv2 {
    fn ne1(&self) -> usize {
        1
    }
    fn ntensor(&self) -> usize {
        10
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
        // dimensions
        shl_shape: [usize; 2],
    ) {
        gto_shell_eval_grid_cart_deriv2(gto, ebuf, coord, alpha, coeff, l, shl_shape);
    }
}
