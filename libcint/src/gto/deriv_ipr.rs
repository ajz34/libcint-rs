use crate::gto::prelude_dev::*;

pub fn gto_shell_eval_grid_cart_ipr(
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

    const COMP_NUM: usize = 9;

    let [nctr, nprim] = shl_shape;
    let ncart = (l + 1) * (l + 2) / 2;
    let nao_to_set = nctr * ncart;
    let mut f0 = [[f64simd::zero(); 3]; ANG_MAX + 3];
    let mut f1 = [[f64simd::zero(); 3]; ANG_MAX + 3];
    let mut f2 = [[f64simd::zero(); 3]; ANG_MAX + 3];
    let mut f3 = [[f64simd::zero(); 3]; ANG_MAX + 3];
    let mut buf = [[f64simd::zero(); 3]; 3];
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

            gto_r_simdd(&mut f1, &f0, l);
            gto_nabla1_simdd(&mut f2, &f0, l + 1, alpha[p]);
            gto_r_simdd(&mut f3, &f2, l);

            for (icart, (lx, ly, lz)) in gto_l_iter(l).enumerate() {
                buf[X][X] = e * f3[lx][X] * f0[ly][Y] * f0[lz][Z];
                buf[X][Y] = e * f2[lx][X] * f1[ly][Y] * f0[lz][Z];
                buf[X][Z] = e * f2[lx][X] * f0[ly][Y] * f1[lz][Z];
                buf[Y][X] = e * f1[lx][X] * f2[ly][Y] * f0[lz][Z];
                buf[Y][Y] = e * f0[lx][X] * f3[ly][Y] * f0[lz][Z];
                buf[Y][Z] = e * f0[lx][X] * f2[ly][Y] * f1[lz][Z];
                buf[Z][X] = e * f1[lx][X] * f0[ly][Y] * f2[lz][Z];
                buf[Z][Y] = e * f0[lx][X] * f1[ly][Y] * f2[lz][Z];
                buf[Z][Z] = e * f0[lx][X] * f0[ly][Y] * f3[lz][Z];

                for k in 0..nctr {
                    let c = f64simd::splat(coeff[k * nprim + p]);
                    gto[0][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[X][X]);
                    gto[1][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[X][Y]);
                    gto[2][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[X][Z]);
                    gto[3][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[Y][X]);
                    gto[4][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[Y][Y]);
                    gto[5][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[Y][Z]);
                    gto[6][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[Z][X]);
                    gto[7][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[Z][Y]);
                    gto[8][k * ncart + icart].get_simdd_mut(g).fma_from(c, buf[Z][Z]);
                }
            }
        }
    }
}

pub struct GtoEvalDerivIpR;
impl GtoEvalAPI for GtoEvalDerivIpR {
    fn ne1(&self) -> usize {
        1
    }
    fn ntensor(&self) -> usize {
        9
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
        gto_shell_eval_grid_cart_ipr(gto, ebuf, coord, alpha, coeff, l, shl_shape);
    }
}
