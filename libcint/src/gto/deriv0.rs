use std::cell::UnsafeCell;

use crate::gto::prelude_dev::*;
use rayon::prelude::*;

/// Compute contracted GTO value, zeroth order.
///
/// This function only handles batch of grids with fixed size [`BLKSIZE`], to be
/// better SIMD optimized.
///
/// This only handles $l = 0$ orbitals.
///
/// # Formula
///
/// $$ $$
///
/// $$
/// \begin{align*}
/// \mathscr G_k (\bm r_g) &= \mathrm{fac} \times \sum_p C_{kp} \mathscr G_p
/// (\bm r_g)
/// \\\\
/// \mathscr G_p (\bm r_g) &= \exp(- \alpha_p r_g^2)
/// \\\\
/// r_g^2 &= x_g^2 + y_g^2 + z_g^2
/// \end{align*}
/// $$
///
/// # Indices Table
///
/// | index | size | notes |
/// |--|--|--|
/// | $k$ | `nctr` | contracted AO |
/// | $p$ | `nprim` | primitive AO |
/// | $g$ | [`BLKSIZE`] | grids |
///
/// # Argument Table
///
/// | variable | formula | dimension | shape | notes |
/// |--|--|--|--|--|
/// | `ectr` | $\phi_k (r_g)$ | $(k, g)$ | `(nctr, BLKSIZE)` | GTO value of contracted AO |
/// | `coord` | $(x_g, y_g, z_g)$ | $(3, g)$ | `(3, BLKSIZE)` | coordinates of grids (taking AO's center as origin) |
/// | `alpha` | $\alpha_p$ | $(p,)$ | `nprim` | exponents of primitive AO |
/// | `coeff` | $C_{kp}$ | $(k, p)$ | `(nctr, BLKSIZE)` | cGTO contraction coefficients |
/// | `fac` | $\mathrm{fac}$ | scalar | | Additional factor (also see [`cint_common_fac_sp`]) |
pub fn gto_contract_exp0(
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

    // r² = x² + y² + z²
    let mut rr = unsafe { f64blk::uninit() };
    for g in 0..BLKSIZE {
        rr[g] = coord[X][g].mul_add(coord[X][g], coord[Y][g].mul_add(coord[Y][g], coord[Z][g] * coord[Z][g]));
    }
    // zero ectr
    for k in 0..nctr {
        for g in 0..BLKSIZE {
            ectr[k][g] = 0.0;
        }
    }
    for p in 0..nprim {
        let mut eprim = unsafe { f64blk::uninit() };
        for g in 0..BLKSIZE {
            let arr = alpha[p] * rr[g];
            eprim[g] = (-arr).exp() * fac;
        }
        // we make the grid index to be the least iterated, different to PySCF
        for k in 0..nctr {
            for g in 0..BLKSIZE {
                ectr[k][g] = eprim[g].mul_add(coeff[k * nprim + p], ectr[k][g]);
            }
        }
    }
}

/// Evaluate GTO shell in Cartesian basis at grids, at arbitary angular momentum
/// $l$.
///
/// # Formula
///
/// $$
/// \phi_\mu (\bm r_g) = x_g^{\mu_x} y_g^{\mu_y} z_g^{\mu_z} \times
/// \mathscr G_{\mu_k} (\bm r_g)
/// $$
///
/// Notes on index $\mu$:
/// - This index represents the cartesian atomic orbital basis. The size of
///   $\mu$ is the multiplication of $n_\mathrm{cart} = (l + 1) (l + 2) / 2$,
///   and $n_\mathrm{ctr}$ (the number of contracted functions in the shell).
/// - $\mu_x$, $\mu_y$, $\mu_z$ are the cartesian exponents of the basis
///   function, satisfying $\mu_x + \mu_y + \mu_z = l$. Note that for a shell of
///   basis, the angular momentum $l$ is fixed for all contracted functions.
/// - The ordering of $(\mu_x, \mu_y, \mu_z)$ follows the convention in PySCF,
///   i.e., the most rapidly changing index is $\mu_z$, then $\mu_y$, and
///   finally $\mu_x$.
/// - $\mu_k$ is the index of contracted function.
///
/// For example of angular momentum $l=2$ (d shell) with $n_\mathrm{ctr} = 2$,
/// the ordering of $\mu$ is:
///
/// | $\mu$ | $(\mu_x, \mu_y, \mu_z)$ | $\mu_k$ || $\mu$ | $(\mu_x, \mu_y, \mu_z)$ | $\mu_k$ |
/// |--|--|--|--|--|--|--|
/// |  0 | (2, 0, 0) | 0 ||  6 | (2, 0, 0) | 1 |
/// |  1 | (1, 1, 0) | 0 ||  7 | (1, 1, 0) | 1 |
/// |  2 | (1, 0, 1) | 0 ||  8 | (1, 0, 1) | 1 |
/// |  3 | (0, 2, 0) | 0 ||  9 | (0, 2, 0) | 1 |
/// |  4 | (0, 1, 1) | 0 || 10 | (0, 1, 1) | 1 |
/// |  5 | (0, 0, 2) | 0 || 11 | (0, 0, 2) | 1 |
///
/// # Indices Table
///
/// | index | size | notes |
/// |--|--|--|
/// | $\mu$ | `ncart * nctr` | cartesian AO basis function |
/// | $\mu_k$ | `nctr` | contracted AO function |
/// | $g$ | `bgrid` | grids, usually equal but may be smaller than [`BLKSIZE`] |
///
/// # Argument Table
///
/// | variable | formula | dimension | shape | notes |
/// |--|--|--|--|--|
/// | `gto` | $\phi_\mu (\bm r_g)$ | $(\mathbb{A}, \mu, g)$ | `(1, nao, ngrid)` | GTO value of cartesian AO basis |
/// | `exps` | $\mathscr G_{\mu_k} (\bm r_g)$ | $(\mu_k, g)$ | `(nctr, BLKSIZE)` | cGTO value (given by function [`gto_contract_exp0`]) |
/// | `coord` | $(x_g, y_g, z_g)$ | $(3, g)$ | `(3, BLKSIZE)` | coordinates of grids (taking AO's center as origin) |
/// | `l` | $l$ | scalar | | angular momentum of the shell |
///
/// Note of `gto`:
/// - The first dimension `\mathbb{A}` is for properties, currently always 1
///   since we are evaluating only the basis value, not its derivatives.
///
/// # Offset Table
///
/// | variable | notes |
/// |--|--|
/// | `iao` | starting index $\mu$ of the shell in `gto` |
/// | `igrid` | starting index $g$ of the grids in `gto` |
pub fn gto_shell_eval_grid_cart(
    // arguments
    gto: &mut [&mut [&mut [f64]]],
    exps: &[f64blk],
    coord: &[f64blk; 3],
    l: usize,
    // dimensions
    nctr: usize,
    bgrid: usize,
    // offsets
    iao: usize,
    igrid: usize,
) {
    const ANG_MAX: usize = crate::ffi::cint_ffi::ANG_MAX as usize;
    unsafe { core::hint::assert_unchecked(bgrid <= BLKSIZE) }

    match l {
        0 => {
            for k in 0..nctr {
                for i in 0..bgrid {
                    gto[0][iao + k][igrid + i] = exps[k][i];
                }
            }
        },
        1 => {
            for k in 0..nctr {
                for g in 0..bgrid {
                    let e = exps[k][g];
                    let x = coord[X][g];
                    let y = coord[Y][g];
                    let z = coord[Z][g];
                    gto[0][iao + 3 * k + X][igrid + g] = x * e;
                    gto[0][iao + 3 * k + Y][igrid + g] = y * e;
                    gto[0][iao + 3 * k + Z][igrid + g] = z * e;
                }
            }
        },
        2 => {
            for k in 0..nctr {
                for g in 0..bgrid {
                    let e = exps[k][g];
                    let x = coord[X][g];
                    let y = coord[Y][g];
                    let z = coord[Z][g];
                    gto[0][iao + 6 * k + XX][igrid + g] = x * x * e;
                    gto[0][iao + 6 * k + XY][igrid + g] = x * y * e;
                    gto[0][iao + 6 * k + XZ][igrid + g] = x * z * e;
                    gto[0][iao + 6 * k + YY][igrid + g] = y * y * e;
                    gto[0][iao + 6 * k + YZ][igrid + g] = y * z * e;
                    gto[0][iao + 6 * k + ZZ][igrid + g] = z * z * e;
                }
            }
        },
        3 => {
            for k in 0..nctr {
                for g in 0..bgrid {
                    let e = exps[k][g];
                    let x = coord[X][g];
                    let y = coord[Y][g];
                    let z = coord[Z][g];
                    gto[0][iao + 10 * k + XXX][igrid + g] = x * x * x * e;
                    gto[0][iao + 10 * k + XXY][igrid + g] = x * x * y * e;
                    gto[0][iao + 10 * k + XXZ][igrid + g] = x * x * z * e;
                    gto[0][iao + 10 * k + XYY][igrid + g] = x * y * y * e;
                    gto[0][iao + 10 * k + XYZ][igrid + g] = x * y * z * e;
                    gto[0][iao + 10 * k + XZZ][igrid + g] = x * z * z * e;
                    gto[0][iao + 10 * k + YYY][igrid + g] = y * y * y * e;
                    gto[0][iao + 10 * k + YYZ][igrid + g] = y * y * z * e;
                    gto[0][iao + 10 * k + YZZ][igrid + g] = y * z * z * e;
                    gto[0][iao + 10 * k + ZZZ][igrid + g] = z * z * z * e;
                }
            }
        },
        _ => {
            let mut pows = unsafe { [[f64blk::uninit(); 3]; ANG_MAX + 1] };
            let ncart = (l + 1) * (l + 2) / 2;
            for i in 0..bgrid {
                pows[0][X][i] = 1.0;
                pows[0][Y][i] = 1.0;
                pows[0][Z][i] = 1.0;
            }
            for lx in 1..=l {
                for i in 0..bgrid {
                    pows[lx][X][i] = pows[lx - 1][X][i] * coord[X][i];
                    pows[lx][Y][i] = pows[lx - 1][Y][i] * coord[Y][i];
                    pows[lx][Z][i] = pows[lx - 1][Z][i] * coord[Z][i];
                }
            }
            for k in 0..nctr {
                let mut icart = 0;
                for lx in (0..=l).rev() {
                    for ly in (0..=(l - lx)).rev() {
                        let lz = l - lx - ly;
                        for i in 0..bgrid {
                            gto[0][iao + ncart * k + icart][igrid + i] = exps[k][i] * pows[lx][X][i] * pows[ly][Y][i] * pows[lz][Z][i];
                        }
                        icart += 1;
                    }
                }
            }
        },
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
/// | `coord` | $\bm r_g$ | $(g, 3)$ | `(ngrid, 3)` | coordinates of grids |
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

#[allow(clippy::too_many_arguments)]
pub fn gto_eval_cart_iter(
    // arguments
    ao: &mut [&mut [&mut [f64]]], // (nprop, nao, ngrid)
    coord: &[[f64; 3]],
    fac: f64,
    shls_slice: [usize; 2],
    ao_loc: &[usize],
    // buffer
    buf: &mut [f64blk],
    // dimensions
    bgrid: usize,
    // offsets
    iao: usize,
    igrid: usize,
    // cint data
    atm: &[[c_int; ATM_SLOTS as usize]],
    bas: &[[c_int; BAS_SLOTS as usize]],
    env: &[f64],
) {
    const ATOM_OF: usize = crate::ffi::cint_ffi::ATOM_OF as usize;
    const NPRIM_MAX: usize = crate::ffi::cint_ffi::NPRIM_MAX as usize;
    const NCTR_MAX: usize = crate::ffi::cint_ffi::NCTR_MAX as usize;
    const NPRIM_OF: usize = crate::ffi::cint_ffi::NPRIM_OF as usize;
    const NCTR_OF: usize = crate::ffi::cint_ffi::NCTR_OF as usize;
    const ANG_OF: usize = crate::ffi::cint_ffi::ANG_OF as usize;
    const PTR_EXP: usize = crate::ffi::cint_ffi::PTR_EXP as usize;
    const PTR_COEFF: usize = crate::ffi::cint_ffi::PTR_COEFF as usize;

    let [sh0, sh1] = shls_slice;
    let [atm0, atm1] = [bas[sh0][ATOM_OF] as usize, bas[sh1 - 1][ATOM_OF] as usize + 1];
    let atm_count = atm1 - atm0;
    let nbuf_grid2atm = atm_count * 3;
    let nbuf_eprim = NPRIM_MAX.max(NCTR_MAX);
    let (grid2atm, buf) = buf.split_at_mut(nbuf_grid2atm);
    let (eprim, _) = buf.split_at_mut(nbuf_eprim);

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

        gto_contract_exp0(eprim, coord, alpha, coeff, fac1, nprim, nctr);
        gto_shell_eval_grid_cart(ao, &eprim[0..nctr], coord, l, nctr, bgrid, iao, igrid);
    }
}

pub fn gto_eval_loop<const NE1: usize, const NTENSOR: usize>(
    ao: &mut [f64], // (nprop, nao, ngrids)
    fac: f64,
    shls_slice: [usize; 2],
    ao_loc: &[usize],
    ngrid: usize,
    coord: &[[f64; 3]],
    atm: &[[c_int; ATM_SLOTS as usize]],
    bas: &[[c_int; BAS_SLOTS as usize]],
    env: &[f64],
) {
    let [sh0, sh1] = shls_slice;
    let nao = ao_loc[sh1] - ao_loc[sh0];
    let nprop = NE1 * NTENSOR;

    // convert `ao` from 1-D `&mut [f64]` to 3-D `&mut [&mut [&mut [f64]]]`
    assert!(ao.len() == nprop * nao * ngrid);
    let mut ao: Vec<Vec<&mut [f64]>> = ao.chunks_exact_mut(nao * ngrid).map(|chunk| chunk.chunks_exact_mut(ngrid).collect()).collect();
    let ao: Vec<&mut [&mut [f64]]> = ao.iter_mut().map(|v| v.as_mut_slice()).collect::<Vec<_>>();
    let ao: &[&mut [&mut [f64]]] = ao.as_slice();

    // split shells by atoms
    let shloc = gto_shloc_by_atom(shls_slice, bas);
    let nshblk = shloc.len() - 1; // usually number of atoms

    // split grids to blocks
    let nblk = ngrid.div_ceil(BLKSIZE);

    const NBUF: usize = if NPRIM_MAX < NCTR_MAX { NCTR_MAX } else { NPRIM_MAX } as usize + NCTR_CART + 3;

    (0..nblk * nshblk).into_par_iter().for_each(|niter| {
        let ishblk = niter / nblk;
        let iblk = niter % nblk;
        let mut buf = [f64blk::zero(); NBUF];

        let [ish0, ish1] = [shloc[ishblk], shloc[ishblk + 1]];
        let bgrid = if (iblk + 1) * BLKSIZE <= ngrid { BLKSIZE } else { ngrid - iblk * BLKSIZE };
        let igrid = iblk * BLKSIZE;
        let iao = ao_loc[ish0] - ao_loc[sh0];

        let ao = unsafe { cast_mut_slice(ao) };
        gto_eval_cart_iter(ao, &coord[igrid..igrid + bgrid], fac, [ish0, ish1], ao_loc, &mut buf, bgrid, iao, igrid, atm, bas, env);
    });
}

#[test]
fn playground() {
    let mut cint_data = init_h2o_def2_tzvp();
    cint_data.set_cint_type("cart");
    let ngrid = 1024;
    let nao = cint_data.nao();
    let mut ao = vec![0.0f64; nao * ngrid];
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();

    gto_eval_loop::<1, 1>(
        &mut ao,
        1.0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        ngrid,
        &coord,
        &cint_data.atm,
        &cint_data.bas,
        &cint_data.env,
    );

    println!("ao[  0..10]: {:?}", &ao[0..10]);
    println!("ao[  0..10]: {:?}", &ao[500..510]);
    println!("ao[nao..10]: {:?}", &ao[ngrid..ngrid + 10]);
    println!("ao[nao..10]: {:?}", &ao[45 * ngrid..45 * ngrid + 10]);
}
