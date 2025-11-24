use crate::gto::prelude_dev::*;
use rayon::prelude::*;

pub trait GTOEvalAPI: Send + Sync {
    fn ne1(&self) -> usize;
    fn ntensor(&self) -> usize;
    fn ncomp(&self) -> usize {
        self.ne1() * self.ntensor()
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
    );
    fn gto_shell_eval(
        &self,
        // arguments
        gto: &mut [f64blk],
        exps: &[f64blk],
        coord: &[f64blk; 3],
        l: usize,
        // dimensions
        nctr: usize,
    );
}

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

pub fn gto_copy_grids<const CART: bool>(
    // arguments
    gto: &[f64blk],
    ao: &mut [f64],
    l: usize,
    // dimensions
    nctr: usize,
    ncomp: usize,
    nao: usize,
    ngrid: usize,
    // offsets
    iao: usize,
    igrid: usize,
) {
    let nang = if CART { (l + 1) * (l + 2) / 2 } else { 2 * l + 1 };
    let nshlbas = nctr * nang;
    if igrid + BLKSIZE < ngrid {
        for a in 0..ncomp {
            for mu in 0..nshlbas {
                let ptr = a * nao * ngrid + (iao + mu) * ngrid + igrid;
                let dst = &mut ao[ptr..ptr + BLKSIZE];
                let src = &gto[a * nshlbas + mu];
                unsafe { src.write_ensure(dst) }
            }
        }
    } else {
        let bgrid = ngrid - igrid;
        for a in 0..ncomp {
            for mu in 0..nshlbas {
                let ptr = a * nao * ngrid + (iao + mu) * ngrid + igrid;
                let dst = &mut ao[ptr..ptr + bgrid];
                let src = &gto[a * nshlbas + mu];
                src.write(dst);
            }
        }
    }
}

/// Get shell location indices where atom changes.
///
/// This will categorize shells by atom, useful for atom-wise operations.
///
/// Note that the returned `shloc` contains the end index of the last shell, so
/// length is $n_\mathrm{atom} + 1$.
pub fn gto_shloc_by_atom(shls_slice: [usize; 2], bas: &[[c_int; BAS_SLOTS as usize]]) -> Vec<usize> {
    const ATOM_OF: usize = crate::ffi::cint_ffi::ATOM_OF as usize;
    let [sh0, sh1] = shls_slice;
    let mut shloc = vec![sh0];
    let mut lastatm = bas[sh0][ATOM_OF] as usize;
    for ish in sh0..sh1 {
        if lastatm != bas[ish][ATOM_OF] as usize {
            lastatm = bas[ish][ATOM_OF] as usize;
            shloc.push(ish);
        }
    }
    shloc.push(sh1);
    shloc
}

/// Compute grid coordinates relative to atoms.
///
/// # Formula
///
/// $$
/// \bm r_g \leftarrow \bm r_g - \bm R_A
/// $$
///
/// # Indices Table
///
/// | index | size | notes |
/// |--|--|--|
/// | $A$ | `natm` | atom index |
/// | $g$ | `bgrid` | grids, usually equal but may be smaller than [`BLKSIZE`] |
///
/// # Argument Table
///
/// | variable | formula | dimension | shape | notes |
/// |--|--|--|--|--|
/// | `grid2atm` | $\bm r_g - \bm R_A$ | $(A, 3, g)$ | `(natm, 3, bgrid)` | grid coordinates relative to atom |
/// | `coord` | $\bm r_g$ | $(g, 3)$ | `(bgrid, 3)` | coordinates of grids |
/// | `atm` | $\bm R_A$ | $(A,)$ | `(natm, 3)` | coordinates of atoms (contained by `atm` field of [`CInt`]) |
///
/// Note of `atm`: Only the necessary slice of `atm` should be passed in, do not
/// pass the entire slice in most use cases.
///
/// # Offset Table
///
/// | variable | notes |
/// |--|--|
/// | `igrid` | starting index $g$ of the grids in `coord` |
pub fn gto_fill_grid2atm(
    // arguments
    grid2atm: &mut [[f64blk; 3]],
    coord: &[[f64; 3]],
    // dimensions
    bgrid: usize,
    // cint data
    atm: &[[c_int; ATM_SLOTS as usize]],
    env: &[f64],
) {
    const PTR_COORD: usize = crate::ffi::cint_ffi::PTR_COORD as usize;
    unsafe { core::hint::assert_unchecked(bgrid <= BLKSIZE) }

    for iatm in 0..atm.len() {
        let coord_ptr = atm[iatm][PTR_COORD] as usize;
        let ratm: [f64; 3] = env[coord_ptr..coord_ptr + 3].try_into().unwrap();
        for g in 0..bgrid {
            grid2atm[iatm][X][g] = coord[g][X] - ratm[X];
            grid2atm[iatm][Y][g] = coord[g][Y] - ratm[Y];
            grid2atm[iatm][Z][g] = coord[g][Z] - ratm[Z];
        }
    }
}

pub fn gto_nctr_cart_max(shls_slice: [usize; 2], bas: &[[c_int; BAS_SLOTS as usize]]) -> usize {
    const ANG_OF: usize = crate::ffi::cint_ffi::ANG_OF as usize;
    const NCTR_OF: usize = crate::ffi::cint_ffi::NCTR_OF as usize;
    let [sh0, sh1] = shls_slice;
    (sh0..sh1)
        .map(|ish| {
            let l = bas[ish][ANG_OF] as usize;
            let nctr = bas[ish][NCTR_OF] as usize;
            let ncart = (l + 1) * (l + 2) / 2;
            nctr * ncart
        })
        .max()
        .unwrap_or(0)
}

pub fn gto_nctr_max(shls_slice: [usize; 2], bas: &[[c_int; BAS_SLOTS as usize]]) -> usize {
    const NCTR_OF: usize = crate::ffi::cint_ffi::NCTR_OF as usize;
    let [sh0, sh1] = shls_slice;
    (sh0..sh1).map(|ish| bas[ish][NCTR_OF] as usize).max().unwrap_or(0)
}

pub fn gto_nprim_max(shls_slice: [usize; 2], bas: &[[c_int; BAS_SLOTS as usize]]) -> usize {
    const NPRIM_OF: usize = crate::ffi::cint_ffi::NPRIM_OF as usize;
    let [sh0, sh1] = shls_slice;
    (sh0..sh1).map(|ish| bas[ish][NPRIM_OF] as usize).max().unwrap_or(0)
}

#[allow(clippy::too_many_arguments)]
pub fn gto_eval_iter<const CART: bool>(
    evaluator: &dyn GTOEvalAPI,
    // arguments
    ao: &mut [f64], // (ncomp, nao, ngrid)
    coord: &[[f64; 3]],
    fac: f64,
    shls_slice: [usize; 2],
    ao_loc: &[usize],
    // buffer
    buf: &mut [f64blk],
    // dimensions
    nao: usize,
    ngrid: usize,
    // offsets
    iao: usize,
    igrid: usize,
    // cint data
    atm: &[[c_int; ATM_SLOTS as usize]],
    bas: &[[c_int; BAS_SLOTS as usize]],
    env: &[f64],
) {
    const ATOM_OF: usize = crate::ffi::cint_ffi::ATOM_OF as usize;
    const NPRIM_OF: usize = crate::ffi::cint_ffi::NPRIM_OF as usize;
    const NCTR_OF: usize = crate::ffi::cint_ffi::NCTR_OF as usize;
    const ANG_OF: usize = crate::ffi::cint_ffi::ANG_OF as usize;
    const PTR_EXP: usize = crate::ffi::cint_ffi::PTR_EXP as usize;
    const PTR_COEFF: usize = crate::ffi::cint_ffi::PTR_COEFF as usize;

    let ncomp = evaluator.ncomp();
    let nctr_cart_max = gto_nctr_cart_max(shls_slice, bas);
    let nctr_max = gto_nctr_max(shls_slice, bas);
    let nprim_max = gto_nprim_max(shls_slice, bas);

    let [sh0, sh1] = shls_slice;
    let [atm0, atm1] = [bas[sh0][ATOM_OF] as usize, bas[sh1 - 1][ATOM_OF] as usize + 1];
    let atm_count = atm1 - atm0;
    let nbuf_grid2atm = atm_count * 3;
    let nbuf_eprim = 2 * nctr_max.max(nprim_max);
    let nbuf_gto = ncomp * nctr_cart_max;
    let (grid2atm, buf) = buf.split_at_mut(nbuf_grid2atm);
    let (eprim, buf) = buf.split_at_mut(nbuf_eprim);
    let (gto_cart, buf) = buf.split_at_mut(nbuf_gto);
    let (gto_sph, _) = buf.split_at_mut(nbuf_gto);
    let bgrid = (ngrid - igrid).min(BLKSIZE);

    let grid2atm = unsafe { grid2atm.as_chunks_unchecked_mut::<3>() };
    gto_fill_grid2atm(grid2atm, coord, bgrid, &atm[atm0..atm1], env);
    for bas_id in sh0..sh1 {
        let nprim = bas[bas_id][NPRIM_OF] as usize;
        let nctr = bas[bas_id][NCTR_OF] as usize;
        let l = bas[bas_id][ANG_OF] as usize;
        let fac1 = fac * cint_common_fac_sp(l as c_int);

        let p_exp = bas[bas_id][PTR_EXP] as usize;
        let p_coeff = bas[bas_id][PTR_COEFF] as usize;
        let atm_id = bas[bas_id][ATOM_OF] as usize;
        let alpha = &env[p_exp..p_exp + nprim];
        let coeff = &env[p_coeff..p_coeff + nprim * nctr];
        let coord = &grid2atm[atm_id - atm0];
        let iao = iao + ao_loc[bas_id] - ao_loc[sh0];

        if CART || l <= 1 {
            evaluator.gto_exp(eprim, coord, alpha, coeff, fac1, nprim, nctr);
            evaluator.gto_shell_eval(gto_cart, eprim, coord, l, nctr);
            gto_copy_grids::<true>(gto_cart, ao, l, nctr, ncomp, nao, ngrid, iao, igrid);
        } else {
            evaluator.gto_exp(eprim, coord, alpha, coeff, fac1, nprim, nctr);
            evaluator.gto_shell_eval(gto_cart, eprim, coord, l, nctr);
            let ncart = (l + 1) * (l + 2) / 2;
            let nsph = 2 * l + 1;
            for a in 0..ncomp {
                for k in 0..nctr {
                    let ptr_sph = &gto_sph[a * nctr * nsph + k * nsph] as *const _ as *mut f64;
                    let ptr_cart = &gto_cart[a * nctr * ncart + k * ncart..] as *const _ as *mut f64;
                    unsafe { CINTc2s_ket_sph1(ptr_sph, ptr_cart, BLKSIZE as c_int, BLKSIZE as c_int, l as c_int) };
                }
            }
            gto_copy_grids::<false>(gto_sph, ao, l, nctr, ncomp, nao, ngrid, iao, igrid);
        }
    }
}

pub fn gto_eval_loop<const CART: bool>(
    evaluator: &dyn GTOEvalAPI,
    ao: &mut [f64], // (ncomp, nao, ngrid)
    coord: &[[f64; 3]],
    fac: f64,
    shls_slice: [usize; 2],
    ao_loc: &[usize],
    atm: &[[c_int; ATM_SLOTS as usize]],
    bas: &[[c_int; BAS_SLOTS as usize]],
    env: &[f64],
) {
    let ncomp = evaluator.ncomp();
    let [sh0, sh1] = shls_slice;
    let nao = ao_loc[sh1] - ao_loc[sh0];
    let ngrid = coord.len();
    assert!(ao.len() >= ncomp * nao * ngrid);

    // split shells by atoms
    let shloc = gto_shloc_by_atom(shls_slice, bas);
    let nshblk = shloc.len() - 1; // usually number of atoms

    // split grids to blocks
    let nblk = ngrid.div_ceil(BLKSIZE);

    // grid2atm: usually 3, since we always split shells by atoms
    // eprim or ectr: 2 * nctr
    // gto_cart: ncomp_max * nctr_cart
    // gto_sph: ncomp_max * nctr_cart
    let nctr_cart_max = gto_nctr_cart_max(shls_slice, bas);
    let nctr_max = gto_nctr_max(shls_slice, bas);
    let nprim_max = gto_nprim_max(shls_slice, bas);
    let nbuf_grid2atm = 3;
    let nbuf_eprim = 2 * nctr_max.max(nprim_max);
    let nbuf_gto = 2 * ncomp * nctr_cart_max;
    let ncache = nbuf_eprim + nbuf_gto + nbuf_grid2atm;
    let nthreads = rayon::current_num_threads();
    let thread_cache = (0..nthreads).map(|_| unsafe { vec![f64blk::uninit(); ncache] }).collect_vec();

    (0..nblk * nshblk).into_par_iter().for_each(|niter| {
        let ishblk = niter / nblk;
        let iblk = niter % nblk;
        let thread_id = rayon::current_thread_index().unwrap_or(0);

        let ish0 = shloc[ishblk];
        let ish1 = shloc[ishblk + 1];
        let igrid = iblk * BLKSIZE;
        let iao = ao_loc[ish0] - ao_loc[sh0];
        let bgrid = (ngrid - igrid).min(BLKSIZE);

        let ao = unsafe { cast_mut_slice(ao) };
        let coord = &coord[igrid..igrid + bgrid];
        let cache = unsafe { cast_mut_slice(&thread_cache[thread_id]) };
        gto_eval_iter::<CART>(evaluator, ao, coord, fac, [ish0, ish1], ao_loc, cache, nao, ngrid, iao, igrid, atm, bas, env);
    });
}
