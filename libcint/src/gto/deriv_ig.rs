use crate::gto::prelude_dev::*;

pub fn gto_shell_eval_grid_cart_ig(
    // arguments
    gto: &mut [f64blk],
    eprim: &[f64blk],
    coord: &[f64blk; 3],
    coeff: &[f64],
    l: usize,
    center: [f64; 3],
    // dimensions
    shl_shape: [usize; 2],
) {
    const ANG_MAX: usize = crate::ffi::cint_ffi::ANG_MAX as usize;

    const COMP_NUM: usize = 3;

    let [nctr, nprim] = shl_shape;
    let ncart = (l + 1) * (l + 2) / 2;
    let nao_to_set = nctr * ncart;
    let c = [f64simd::splat(-center[X]), f64simd::splat(-center[Y]), f64simd::splat(-center[Z])];
    let mut f0 = [[f64simd::zero(); 3]; ANG_MAX + 2];
    let mut f1 = [[f64simd::zero(); 3]; ANG_MAX + 2];
    let mut s = [f64simd::zero(); 3];
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
            for ll in 1..=l + 1 {
                f0[ll][X] = f0[ll - 1][X] * x;
                f0[ll][Y] = f0[ll - 1][Y] * y;
                f0[ll][Z] = f0[ll - 1][Z] * z;
            }
            gto_x1_simdd(&mut f1, &f0, l, center);
            for (icart, (lx, ly, lz)) in gto_l_iter(l).enumerate() {
                s[X] = e * f1[lx][X] * f0[ly][Y] * f0[lz][Z];
                s[Y] = e * f0[lx][X] * f1[ly][Y] * f0[lz][Z];
                s[Z] = e * f0[lx][X] * f0[ly][Y] * f1[lz][Z];
                buf[X] = s[Y] * c[Z] - s[Z] * c[Y];
                buf[Y] = s[Z] * c[X] - s[X] * c[Z];
                buf[Z] = s[X] * c[Y] - s[Y] * c[X];
                for k in 0..nctr {
                    let c = f64simd::splat(coeff[k * nprim + p]);
                    gto[0][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[X]);
                    gto[1][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[Y]);
                    gto[2][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[Z]);
                }
            }
        }
    }
}

pub struct GtoEvalDerivIg;
impl GtoEvalAPI for GtoEvalDerivIg {
    fn ne1(&self) -> usize {
        1
    }
    fn ntensor(&self) -> usize {
        3
    }
    fn gto_exp(&self, ebuf: &mut [f64blk], coord: &[f64blk; 3], alpha: &[f64], _coeff: &[f64], fac: f64, shl_shape: [usize; 2]) {
        let eprim = ebuf;
        let [_nctr, nprim] = shl_shape;
        const IG_FAC: f64 = 0.5;
        gto_prim_exp(eprim, coord, alpha, IG_FAC * fac, nprim);
    }
    fn gto_shell_eval(
        &self,
        // arguments
        gto: &mut [f64blk],
        ebuf: &[f64blk],
        coord: &[f64blk; 3],
        _alpha: &[f64],
        coeff: &[f64],
        l: usize,
        center: [f64; 3],
        // dimensions
        shl_shape: [usize; 2],
    ) {
        gto_shell_eval_grid_cart_ig(gto, ebuf, coord, coeff, l, center, shl_shape);
    }
}
