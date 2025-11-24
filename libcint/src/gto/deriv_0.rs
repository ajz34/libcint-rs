use crate::gto::prelude_dev::*;

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
    gto: &mut [f64blk],
    exps: &[f64blk],
    coord: &[f64blk; 3],
    l: usize,
    // dimensions
    nctr: usize,
) {
    const ANG_MAX: usize = crate::ffi::cint_ffi::ANG_MAX as usize;

    match l {
        0 => {
            for k in 0..nctr {
                for g in 0..BLKSIZE {
                    gto[k][g] = exps[k][g];
                }
            }
        },
        1 => {
            for k in 0..nctr {
                for g in 0..BLKSIMD {
                    let e = exps[k].get_simd(g);
                    let x = coord[X].get_simd(g);
                    let y = coord[Y].get_simd(g);
                    let z = coord[Z].get_simd(g);

                    *gto[3 * k + X].get_simd_mut(g) = x * e;
                    *gto[3 * k + Y].get_simd_mut(g) = y * e;
                    *gto[3 * k + Z].get_simd_mut(g) = z * e;
                }
            }
        },
        2 => {
            for k in 0..nctr {
                for g in 0..BLKSIMD {
                    let e = exps[k].get_simd(g);
                    let x = coord[X].get_simd(g);
                    let y = coord[Y].get_simd(g);
                    let z = coord[Z].get_simd(g);

                    *gto[6 * k + XX].get_simd_mut(g) = x * x * e;
                    *gto[6 * k + XY].get_simd_mut(g) = x * y * e;
                    *gto[6 * k + XZ].get_simd_mut(g) = x * z * e;
                    *gto[6 * k + YY].get_simd_mut(g) = y * y * e;
                    *gto[6 * k + YZ].get_simd_mut(g) = y * z * e;
                    *gto[6 * k + ZZ].get_simd_mut(g) = z * z * e;
                }
            }
        },
        3 => {
            for k in 0..nctr {
                for g in 0..BLKSIMD {
                    let e = exps[k].get_simd(g);
                    let x = coord[X].get_simd(g);
                    let y = coord[Y].get_simd(g);
                    let z = coord[Z].get_simd(g);

                    *gto[10 * k + XXX].get_simd_mut(g) = x * x * x * e;
                    *gto[10 * k + XXY].get_simd_mut(g) = x * x * y * e;
                    *gto[10 * k + XXZ].get_simd_mut(g) = x * x * z * e;
                    *gto[10 * k + XYY].get_simd_mut(g) = x * y * y * e;
                    *gto[10 * k + XYZ].get_simd_mut(g) = x * y * z * e;
                    *gto[10 * k + XZZ].get_simd_mut(g) = x * z * z * e;
                    *gto[10 * k + YYY].get_simd_mut(g) = y * y * y * e;
                    *gto[10 * k + YYZ].get_simd_mut(g) = y * y * z * e;
                    *gto[10 * k + YZZ].get_simd_mut(g) = y * z * z * e;
                    *gto[10 * k + ZZZ].get_simd_mut(g) = z * z * z * e;
                }
            }
        },
        _ => {
            let mut pows = unsafe { [[f64blk::uninit(); 3]; ANG_MAX + 1] };
            let ncart = (l + 1) * (l + 2) / 2;
            for g in 0..BLKSIMD {
                pows[0][X].get_simd_mut(g).fill(1.0);
                pows[0][Y].get_simd_mut(g).fill(1.0);
                pows[0][Z].get_simd_mut(g).fill(1.0);
            }
            for lx in 1..=l {
                for g in 0..BLKSIMD {
                    *pows[lx][X].get_simd_mut(g) = pows[lx - 1][X].get_simd(g) * coord[X].get_simd(g);
                    *pows[lx][Y].get_simd_mut(g) = pows[lx - 1][Y].get_simd(g) * coord[Y].get_simd(g);
                    *pows[lx][Z].get_simd_mut(g) = pows[lx - 1][Z].get_simd(g) * coord[Z].get_simd(g);
                }
            }
            for k in 0..nctr {
                let mut icart = 0;
                for lx in (0..=l).rev() {
                    for ly in (0..=(l - lx)).rev() {
                        let lz = l - lx - ly;
                        for g in 0..BLKSIMD {
                            *gto[ncart * k + icart].get_simd_mut(g) =
                                exps[k].get_simd(g) * pows[lx][X].get_simd(g) * pows[ly][Y].get_simd(g) * pows[lz][Z].get_simd(g);
                        }
                        icart += 1;
                    }
                }
            }
        },
    }
}

pub struct GTOEvalDeriv0 {}
impl GTOEvalAPI for GTOEvalDeriv0 {
    fn ne1(&self) -> usize {
        1
    }
    fn ntensor(&self) -> usize {
        1
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
        gto_contract_exp0(ectr, coord, alpha, coeff, fac, nprim, nctr);
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
        gto_shell_eval_grid_cart(gto, exps, coord, l, nctr);
    }
}

#[test]
fn playground_cart() {
    let mut cint_data = init_h2o_def2_tzvp();
    cint_data.set_cint_type("cart");
    let ngrid = 1024;
    let nao = cint_data.nao();
    let mut ao = vec![0.0f64; nao * ngrid];
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();

    gto_eval_loop::<true>(
        &GTOEvalDeriv0 {},
        &mut ao,
        &coord,
        1.0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        None,
        &cint_data.atm,
        &cint_data.bas,
        &cint_data.env,
    );

    println!("ao[  0..10]: {:?}", &ao[0..10]);
    println!("ao[  0..10]: {:?}", &ao[50..60]);
    println!("ao[nao..10]: {:?}", &ao[45 * ngrid..45 * ngrid + 10]);
    println!("fp(ao): {:}", cint_fp(&ao));
}

#[test]
fn playground_sph() {
    let cint_data = init_h2o_def2_tzvp();
    let ngrid = 1024;
    let nao = cint_data.nao();
    let mut ao = vec![0.0f64; nao * ngrid];
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    println!("ao_loc: {:?}", cint_data.ao_loc());

    gto_eval_loop::<false>(
        &GTOEvalDeriv0 {},
        &mut ao,
        &coord,
        1.0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        None,
        &cint_data.atm,
        &cint_data.bas,
        &cint_data.env,
    );

    for i in 0..40 {
        println!("ao[{}]: {:?}", i, &ao[i * ngrid..i * ngrid + 10]);
    }
    println!("fp(ao): {:}", cint_fp(&ao));
    println!("nao: {:}", nao);
}

#[test]
fn test_c10h22_cart() {
    let mut cint_data = init_c10h22_def2_qzvp();
    cint_data.set_cint_type("cart");
    let ngrid = 2048;
    let nao = cint_data.nao();
    let mut ao = vec![0.0f64; nao * ngrid];
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    gto_eval_loop::<true>(
        &GTOEvalDeriv0 {},
        &mut ao,
        &coord,
        1.0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        None,
        &cint_data.atm,
        &cint_data.bas,
        &cint_data.env,
    );
    let fp_ao = cint_fp(&ao);
    println!("fp_ao = {:}", fp_ao);
    println!("ao[  0..10]: {:?}", &ao[0..10]);
    println!("ao[  0..10]: {:?}", &ao[500..510]);
    println!("ao[nao..10]: {:?}", &ao[ngrid..ngrid + 10]);
    println!("ao[nao..10]: {:?}", &ao[45 * ngrid..45 * ngrid + 10]);
    assert!((fp_ao - 198.05403491229913).abs() < 1e-10);

    let ngrid = 1048576;
    let mut ao = vec![0.0f64; nao * ngrid];
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    let time = std::time::Instant::now();
    gto_eval_loop::<true>(
        &GTOEvalDeriv0 {},
        &mut ao,
        &coord,
        1.0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        None,
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
fn test_c10h22_sph() {
    let cint_data = init_c10h22_def2_qzvp();
    let ngrid = 2048;
    let nao = cint_data.nao();
    println!("nao: {:}", nao);
    let mut ao = vec![0.0f64; nao * ngrid];
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    gto_eval_loop::<false>(
        &GTOEvalDeriv0 {},
        &mut ao,
        &coord,
        1.0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        None,
        &cint_data.atm,
        &cint_data.bas,
        &cint_data.env,
    );
    let fp_ao = cint_fp(&ao);
    println!("fp_ao = {:}", fp_ao);
    println!("ao[  0..10]: {:?}", &ao[0..10]);
    println!("ao[  0..10]: {:?}", &ao[500..510]);
    println!("ao[nao..10]: {:?}", &ao[ngrid..ngrid + 10]);
    println!("ao[nao..10]: {:?}", &ao[45 * ngrid..45 * ngrid + 10]);
    assert!((fp_ao - 378.9344483361458).abs() < 1e-10);

    let ngrid = 1048576;
    let mut ao = vec![0.0f64; nao * ngrid];
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    let time = std::time::Instant::now();
    gto_eval_loop::<false>(
        &GTOEvalDeriv0 {},
        &mut ao,
        &coord,
        1.0,
        [0, cint_data.nbas()],
        &cint_data.ao_loc(),
        None,
        &cint_data.atm,
        &cint_data.bas,
        &cint_data.env,
    );
    let elapsed = time.elapsed();
    println!("time: {:.3} s", elapsed.as_secs_f64());
    let fp_ao = cint_fp(&ao);
    println!("fp_ao (1M grids) = {:}", fp_ao);
}
