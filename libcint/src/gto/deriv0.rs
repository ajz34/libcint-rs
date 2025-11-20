use crate::gto::prelude_dev::*;

const ANG_MAX_USIZE: usize = ANG_MAX as usize;

pub fn gto_contract_exp0(ectr: &mut [BlkF64], coord: &[BlkF64; 3], alpha: &[f64], coeff: &[f64], nprim: usize, nctr: usize, fac: f64) {
    const X: usize = 0;
    const Y: usize = 1;
    const Z: usize = 2;

    let mut rr = BlkF64::zero();
    for i in 0..BLKSIZE {
        rr[i] = coord[X][i] * coord[X][i] + coord[Y][i] * coord[Y][i] + coord[Z][i] * coord[Z][i];
    }
    for k in 0..nctr {
        for i in 0..BLKSIZE {
            ectr[k][i] = 0.0;
        }
    }
    for j in 0..nprim {
        for i in 0..BLKSIZE {
            let arr = alpha[j] * rr[i];
            let eprim = (-arr).exp() * fac;
            for k in 0..nctr {
                ectr[k][i] += eprim * coeff[k * nprim + j];
            }
        }
    }
}

#[allow(clippy::too_many_arguments)]
pub fn gto_shell_eval_grid_cart(
    gto: &mut [&mut [&mut [f64]]],
    exps: &[BlkF64],
    coord: &[BlkF64; 3],
    l: usize,
    nctr: usize,
    iao: usize,
    igrid: usize,
    bgrid: usize,
) {
    unsafe { core::hint::assert_unchecked(bgrid <= BLKSIZE) }

    let mut pows = unsafe { [[BlkF64::uninit(); ANG_MAX_USIZE + 1]; 3] };

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
                for i in 0..bgrid {
                    gto[0][iao + 3 * k + X][igrid + i] = coord[X][i] * exps[k][i];
                    gto[0][iao + 3 * k + Y][igrid + i] = coord[Y][i] * exps[k][i];
                    gto[0][iao + 3 * k + Z][igrid + i] = coord[Z][i] * exps[k][i];
                }
            }
        },
        2 => {
            for k in 0..nctr {
                for i in 0..bgrid {
                    let e = exps[k][i];
                    let x = coord[X][i];
                    let y = coord[Y][i];
                    let z = coord[Z][i];
                    gto[0][iao + 6 * k + XX][igrid + i] = e * x * x;
                    gto[0][iao + 6 * k + XY][igrid + i] = e * x * y;
                    gto[0][iao + 6 * k + XZ][igrid + i] = e * x * z;
                    gto[0][iao + 6 * k + YY][igrid + i] = e * y * y;
                    gto[0][iao + 6 * k + YZ][igrid + i] = e * y * z;
                    gto[0][iao + 6 * k + ZZ][igrid + i] = e * z * z;
                }
            }
        },
        _ => {
            let ncart = (l + 1) * (l + 2) / 2;
            for k in 0..nctr {
                for i in 0..bgrid {
                    pows[X][0][i] = 1.0;
                    pows[Y][0][i] = 1.0;
                    pows[Z][0][i] = 1.0;
                }
                for lx in 1..=l {
                    for i in 0..bgrid {
                        pows[X][lx][i] = pows[X][lx - 1][i] * coord[X][i];
                        pows[Y][lx][i] = pows[Y][lx - 1][i] * coord[Y][i];
                        pows[Z][lx][i] = pows[Z][lx - 1][i] * coord[Z][i];
                    }
                }
                let mut icart = 0;
                for lx in (0..=l).rev() {
                    for ly in (0..=(l - lx)).rev() {
                        let lz = l - lx - ly;
                        for i in 0..bgrid {
                            gto[0][iao + ncart * k + icart][igrid + i] = exps[k][i] * pows[X][lx][i] * pows[Y][ly][i] * pows[Z][lz][i];
                        }
                        icart += 1;
                    }
                }
            }
        },
    }
}

pub fn gto_shloc_by_atom(shls_slice: [usize; 2], bas: &[[c_int; BAS_SLOTS as usize]]) -> Vec<usize> {
    const ATOM_OF: usize = crate::ffi::cint_ffi::ATOM_OF as usize;
    let [sh0, sh1] = shls_slice;
    let mut shloc = vec![];
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

pub fn gto_fill_grid2atm(
    grid2atm: &mut [[BlkF64; 3]],
    coord: &[[f64; 3]],
    igrid: usize,
    bgrid: usize,
    atm: &[[c_int; ATM_SLOTS as usize]],
    env: &[f64],
) {
    const PTR_COORD: usize = crate::ffi::cint_ffi::PTR_COORD as usize;
    unsafe { core::hint::assert_unchecked(bgrid <= BLKSIZE) }

    for atm_id in 0..atm.len() {
        let coord_ptr = atm[atm_id][PTR_COORD] as usize;
        let ratm: [f64; 3] = env[coord_ptr..coord_ptr + 3].try_into().unwrap();
        for i in 0..bgrid {
            grid2atm[atm_id][0][i] = coord[igrid + i][0] - ratm[0];
            grid2atm[atm_id][1][i] = coord[igrid + i][1] - ratm[1];
            grid2atm[atm_id][2][i] = coord[igrid + i][2] - ratm[2];
        }
    }
}

pub fn gto_eval_cart_iter(
    ao: &mut [&mut [&mut [f64]]], // (nprop, nao, ngrids)
    fac: f64,
    ngrids: usize,
    igrid: usize,
    shls_slice: [usize; 2],
    ao_loc: &[usize],
    coord: &[[f64; 3]],
    atm: &[[c_int; ATM_SLOTS as usize]],
    bas: &[[c_int; BAS_SLOTS as usize]],
    env: &[f64],
    buf: &mut [BlkF64],
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
    // let mut grid2atm = unsafe { vec![[BlkF64::uninit(); 3]; atm_count] };
    let nbuf_grid2atm = atm_count * 3;
    let nbuf_eprim = NPRIM_MAX.max(NCTR_MAX);
    let (grid2atm, buf) = buf.split_at_mut(nbuf_grid2atm);
    let (eprim, _) = buf.split_at_mut(nbuf_eprim);

    let grid2atm = unsafe { grid2atm.as_chunks_unchecked_mut::<3>() };
    gto_fill_grid2atm(grid2atm, coord, igrid, ngrids, &atm[atm0..atm1], env);
    for bas_id in sh0..sh1 {
        let nprim = bas[bas_id][NPRIM_OF] as usize;
        let nctr = bas[bas_id][NCTR_OF] as usize;
        let l = bas[bas_id][ANG_OF] as usize;
        let fac1 = fac * unsafe { CINTcommon_fac_sp(l as c_int) };

        let p_exp = bas[bas_id][PTR_EXP] as usize;
        let p_coeff = bas[bas_id][PTR_COEFF] as usize;
        let atm_id = bas[bas_id][ATOM_OF] as usize;
        let alpha = &env[p_exp..p_exp + nprim];
        let coeff = &env[p_coeff..p_coeff + nprim * nctr];
        let coord = &grid2atm[atm_id - atm0];
        let iao = ao_loc[bas_id] - ao_loc[sh0];

        gto_contract_exp0(eprim, coord, alpha, coeff, nprim, nctr, fac1);
        gto_shell_eval_grid_cart(ao, &eprim[0..nctr], coord, l, nctr, iao, igrid, ngrids);
    }
}

#[test]
fn playground() {
    let mut cint_data = init_h2o_def2_tzvp();
    cint_data.set_cint_type("cart");
    let ngrids = 56;
    let nao = cint_data.nao();
    let mut ao = vec![0.0f64; nao * ngrids];
    let mut ao_slices: Vec<&mut [f64]> = ao.chunks_exact_mut(ngrids).collect();
    let ao_refs: &mut [&mut [&mut [f64]]] = &mut [ao_slices.as_mut_slice()];
    let coord: Vec<[f64; 3]> = (0..ngrids).map(|i| [i as f64 * 0.1, i as f64 * 0.2, i as f64 * 0.3]).collect();
    let mut buf = vec![BlkF64::zero(); (cint_data.natm() * 3) + NPRIM_MAX.max(NCTR_MAX) as usize];
    gto_eval_cart_iter(
        ao_refs,
        1.0,
        ngrids,
        0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        &coord,
        &cint_data.atm,
        &cint_data.bas,
        &cint_data.env,
        &mut buf,
    );
    println!("ao[0..10]: {:?}", &ao[0..10]);
}

pub fn gto_eval_loop(
    ao: &mut [f64], // (nprop, nao, ngrids)
    fac: f64,
    ngrids: usize,
    param: [usize; 2],
    shls_slice: [usize; 2],
    ao_loc: &[usize],
    coord: &[[f64; 3]],
    atm: &[[c_int; ATM_SLOTS as usize]],
    bas: &[[c_int; BAS_SLOTS as usize]],
    env: &[f64],
) {
    let [sh0, sh1] = shls_slice;
    let nao = ao_loc[sh1] - ao_loc[sh0];
    let nprop = param[POS_E1] * param[TENSOR];

    assert!(ao.len() == nprop * nao * ngrids);
    let mut ao_slices: Vec<Vec<&mut [f64]>> =
        ao.chunks_exact_mut(nao * ngrids).map(|chunk| chunk.chunks_exact_mut(ngrids).collect()).collect();
    let ao_mut_refs: &mut [&mut [&mut [f64]]] = ao_slices.iter_mut().map(|v| v.as_mut_slice()).collect::<Vec<_>>().as_mut_slice();

    let shloc = gto_shloc_by_atom(shls_slice, bas);
    let nshblk = shloc.len() - 1;
    let nblk = ngrids.div_ceil(BLKSIZE);

    let nbuf = NPRIM_MAX.max(NCTR_MAX) as usize + NCTR_CART * nprop;

    todo!()
}
