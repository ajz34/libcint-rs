use crate::gto::prelude_dev::*;

/// Compute **early-contracted** exponental part of GTO without angular momentum
/// and normalization in SIMD blocks of a shell.
///
/// # Formula
///
/// This function evaluates only the exponent part (without polynomial part) of
/// GTO, which is also the same to the $l = 0$ case.
///
/// Given $l = 0$, the GTO value of contracted function $k$ at grid point $g$ is
///
/// $$
/// \phi_{k g}^{l = 0} = \phi_k^{l = 0} (\bm r_g) = \sum_p C_{k p} | \bm 0,
/// \alpha_p, \bm r_g \rangle = \sum_p C_{k p} \exp(- \alpha_p r_g^2)
/// $$
///
/// where $r_g^2 = x_g^2 + y_g^2 + z_g^2$.
///
/// This function only handles $l = 0$ case. For even higher angular momenta
/// $\bm l = (l_x, l_y, l_z)$, the contracted GTO can be evaluated by
/// multiplying the polynomial part (also see [`gto_shell_eval_grid_cart`]):
///
/// $$
/// \phi_{k g} = x_g^{l_x} y_g^{l_y} z_g^{l_z} \times \phi_{k g}^{l = 0}
/// $$
///
/// This function also allows to include an additional factor `fac`.
///
/// # Indices Table
///
/// | index | size | notes |
/// |--|--|--|
/// | $k$ | `nctr` | GTO contraction |
/// | $p$ | `nprim` | primitive GTO |
/// | $g$ | [`BLKSIZE`] | grid block in iterations |
///
/// # Argument Table
///
/// | variable | formula | dimension | shape | notes |
/// |--|--|--|--|--|
/// | `ectr` | $\phi_{kg}^{l = 0}$ | $(k, g)$ | `(nctr, BLKSIZE)` | contracted exponential part of GTO |
/// | `coord` | $(x_g, y_g, z_g)$ | $(3, g)$ | `(3, BLKSIZE)` | coordinates of grids (taking GTO's center as origin) |
/// | `alpha` | $\alpha_p$ | $(p,)$ | `nprim` | exponents of primitive GTO |
/// | `coeff` | $C_{kp}$ | $(k, p)$ | `(nctr, BLKSIZE)` | cGTO contraction coefficients |
/// | `fac` | $\mathrm{fac}$ | | scalar | Additional factor (also see [`cint_common_fac_sp`]) |
///
/// # See also
///
/// This function is a low-level function to implement the
/// [`GtoEvalAPI::gto_exp`] trait method.
///
/// Also see [`gto_contract_exp1`] for handling first derivative of
/// early-contracted GTO exponents.
///
/// This function should be used together with [`gto_shell_eval_grid_cart`].
///
/// # PySCF equivalent
///
/// `libcgto.so`: `int GTOcontract_exp0`
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
/// This functions starts from contracted GTO exponents (without polynomial
/// part):
///
/// $$
/// \phi_{kg}^{l = 0} = \sum_p C_{k p} | \bm 0, \alpha_p, \bm r_g \rangle
/// $$
///
/// The full GTO value is evaluated by multiplying the polynomial part:
///
/// $$
/// \begin{align*}
/// \phi_{\mu g} &= \phi_\mu (\bm r_g) = \sum_p C_{k p} | \bm l, \alpha_p, \bm
/// r_g \rangle = \sum_p C_{k p} x^{l_x} y^{l_y} z^{l_z} | \bm 0, \alpha_p, \bm
/// r_g \rangle
/// \\\\
/// &= x_g^{l_x} y_g^{l_y} z_g^{l_z} \times \phi_{kg}^{l = 0}
/// \end{align*}
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
/// | $\mu$ | `ncart * nctr` | cartesian GTO basis function |
/// | $k$ | `nctr` | GTO contraction |
/// | $g$ | [`BLKSIZE`] | grid block in iterations |
///
/// # Argument Table
///
/// | variable | formula | dimension | shape | notes |
/// |--|--|--|--|--|
/// | `gto` | $\phi_{\mu g}$ | $(\mu, g)$ | `(ncart * nctr, ngrid)` | GTO value by grid block |
/// | `exps` | $\phi_{kg}^{l = 0}$ | $(k, g)$ | `(nctr, BLKSIZE)` | early-contracted GTO exponents <br>(given by function [`gto_contract_exp0`]) |
/// | `coord` | $(x_g, y_g, z_g)$ | $(3, g)$ | `(3, BLKSIZE)` | coordinates of grids (taking GTO's center as origin) |
/// | `l` | $l$ | | scalar | angular momentum of the shell |
///
/// # See also
///
/// This function is a low-level function to implement the
/// [`GtoEvalAPI::gto_shell_eval`] trait method.
///
/// This function should be used together with [`gto_contract_exp0`].
///
/// # PySCF equivalent
///
/// `libcgto.so`: `int GTOshell_eval_grid_cart`
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
                for g in 0..BLKSIMDD {
                    let e = exps[k].get_simdd(g);
                    let x = coord[X].get_simdd(g);
                    let y = coord[Y].get_simdd(g);
                    let z = coord[Z].get_simdd(g);

                    *gto[3 * k + X].get_simdd_mut(g) = x * e;
                    *gto[3 * k + Y].get_simdd_mut(g) = y * e;
                    *gto[3 * k + Z].get_simdd_mut(g) = z * e;
                }
            }
        },
        2 => {
            for k in 0..nctr {
                for g in 0..BLKSIMDD {
                    let e = exps[k].get_simdd(g);
                    let x = coord[X].get_simdd(g);
                    let y = coord[Y].get_simdd(g);
                    let z = coord[Z].get_simdd(g);

                    *gto[6 * k + XX].get_simdd_mut(g) = x * x * e;
                    *gto[6 * k + XY].get_simdd_mut(g) = x * y * e;
                    *gto[6 * k + XZ].get_simdd_mut(g) = x * z * e;
                    *gto[6 * k + YY].get_simdd_mut(g) = y * y * e;
                    *gto[6 * k + YZ].get_simdd_mut(g) = y * z * e;
                    *gto[6 * k + ZZ].get_simdd_mut(g) = z * z * e;
                }
            }
        },
        3 => {
            for k in 0..nctr {
                for g in 0..BLKSIMDD {
                    let e = exps[k].get_simdd(g);
                    let x = coord[X].get_simdd(g);
                    let y = coord[Y].get_simdd(g);
                    let z = coord[Z].get_simdd(g);

                    *gto[10 * k + XXX].get_simdd_mut(g) = x * x * x * e;
                    *gto[10 * k + XXY].get_simdd_mut(g) = x * x * y * e;
                    *gto[10 * k + XXZ].get_simdd_mut(g) = x * x * z * e;
                    *gto[10 * k + XYY].get_simdd_mut(g) = x * y * y * e;
                    *gto[10 * k + XYZ].get_simdd_mut(g) = x * y * z * e;
                    *gto[10 * k + XZZ].get_simdd_mut(g) = x * z * z * e;
                    *gto[10 * k + YYY].get_simdd_mut(g) = y * y * y * e;
                    *gto[10 * k + YYZ].get_simdd_mut(g) = y * y * z * e;
                    *gto[10 * k + YZZ].get_simdd_mut(g) = y * z * z * e;
                    *gto[10 * k + ZZZ].get_simdd_mut(g) = z * z * z * e;
                }
            }
        },
        _ => {
            let mut pows = unsafe { [[f64blk::uninit(); 3]; ANG_MAX + 1] };
            let ncart = (l + 1) * (l + 2) / 2;
            for g in 0..BLKSIMDD {
                pows[0][X].get_simdd_mut(g).fill(1.0);
                pows[0][Y].get_simdd_mut(g).fill(1.0);
                pows[0][Z].get_simdd_mut(g).fill(1.0);
            }
            for lx in 1..=l {
                for g in 0..BLKSIMDD {
                    *pows[lx][X].get_simdd_mut(g) = pows[lx - 1][X].get_simdd(g) * coord[X].get_simdd(g);
                    *pows[lx][Y].get_simdd_mut(g) = pows[lx - 1][Y].get_simdd(g) * coord[Y].get_simdd(g);
                    *pows[lx][Z].get_simdd_mut(g) = pows[lx - 1][Z].get_simdd(g) * coord[Z].get_simdd(g);
                }
            }
            for k in 0..nctr {
                for (icart, (lx, ly, lz)) in gto_l_iter(l).enumerate() {
                    for g in 0..BLKSIMDD {
                        *gto[ncart * k + icart].get_simdd_mut(g) =
                            exps[k].get_simdd(g) * pows[lx][X].get_simdd(g) * pows[ly][Y].get_simdd(g) * pows[lz][Z].get_simdd(g);
                    }
                }
            }
        },
    }
}

pub struct GtoEvalDeriv0;
impl GtoEvalAPI for GtoEvalDeriv0 {
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
        shl_shape: [usize; 2],
    ) {
        let ectr = ebuf;
        let [nctr, nprim] = shl_shape;
        gto_contract_exp0(ectr, coord, alpha, coeff, fac, nprim, nctr);
    }
    fn gto_shell_eval(
        &self,
        // arguments
        gto: &mut [f64blk],
        exps: &[f64blk],
        coord: &[f64blk; 3],
        _alpha: &[f64],
        _coeff: &[f64],
        l: usize,
        _center: [f64; 3],
        // dimensions
        shl_shape: [usize; 2],
    ) {
        let [nctr, _nprim] = shl_shape;
        gto_shell_eval_grid_cart(gto, exps, coord, l, nctr);
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

#[test]
fn playground_cart() {
    let mut mol = init_h2o_def2_tzvp();
    mol.set_cint_type("cart");
    let ngrid = 1024;
    let nao = mol.nao();
    let mut ao = vec![0.0f64; nao * ngrid];
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();

    gto_eval_loop(&mol, &GtoEvalDeriv0, &mut ao, &coord, 1.0, [0, mol.nbas()], None, true).unwrap();

    println!("ao[  0..10]: {:?}", &ao[0..10]);
    println!("ao[  0..10]: {:?}", &ao[50..60]);
    println!("ao[nao..10]: {:?}", &ao[45 * ngrid..45 * ngrid + 10]);
    println!("fp(ao): {:}", cint_fp(&ao));
}
