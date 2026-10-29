//! Regression test for the adaptive radial quadrature in the ECP integrals.
//!
//! Adapted from PySCF `pyscf/gto/test/test_ecp.py` (commit 322e5fa02,
//! PR #3250). At large combined exponents the ECP radial integrand is sharply
//! peaked at small r, and two successive coarse Gauss-Chebyshev rules could
//! agree to 1e-12 even when both were under-resolved. The ECP code therefore
//! has to require two consecutive converged levels before dropping a
//! primitive pair from refinement.
//!
//! For a single-primitive same-center AO with angular momentum `l` and
//! exponent `alpha`, and a one-term local or semilocal ECP channel
//! `c * r^(n-2) * exp(-g r^2)`, the radial integral factorises so that
//!
//! ```text
//! I(g1) / I(g2) = ((2 alpha + g2) / (2 alpha + g1)) ** ((n + 2 l + 1) / 2)
//! ```
//!
//! independently of the AO normalisation. This provides a cheap closed-form
//! check on the radial quadrature.

use libcint::prelude::*;

/// Matrix element `<0|ECPscalar|0>` for a single primitive AO and one ECP
/// projector term on a single atom at the origin.
///
/// `l_ecp = -1` is the local (UL) channel; otherwise it is the semilocal
/// channel for angular momentum `l_ecp`.
fn ecp_scalar(l_ao: i32, alpha: f64, n: i32, g: f64, l_ecp: i32) -> f64 {
    // atm: [charge, ptr_coord, nuc_mod, ptr_zeta, 0, 0]
    // bas: [atom, l, nprim, nctr, kappa, ptr_exp, ptr_coeff, 0]
    // ecpbas: [atom, l, nprim, radipower, so_type, ptr_exp, ptr_coeff, 0]
    let atm = vec![[0, 20, 1, 23, 0, 0]];
    let bas = vec![[0, l_ao, 1, 1, 0, 24, 25, 0]];
    let ecpbas = vec![[0, l_ecp, 1, n, 0, 26, 27, 0]];
    let mut env = vec![0.0; 28];
    env[20] = 0.0; // x
    env[21] = 0.0; // y
    env[22] = 0.0; // z
    env[23] = 0.0; // zeta (unused for point nucleus)
    env[24] = alpha; // AO exponent
    env[25] = 1.0; // AO coefficient
    env[26] = g; // ECP exponent
    env[27] = 1.0; // ECP coefficient
    let cint = CInt { atm, bas, ecpbas, env, cint_type: CIntType::Spheric };
    let (out, _) = cint.integrate("ECPscalar", None, None).into();
    out[0]
}

#[test]
fn ecp_radial_large_exponent_closed_form() {
    for l in 0..=2 {
        for alpha in [1e0f64, 1e2, 1e4, 3e5] {
            for n in 1..=2 {
                for (g1, g2) in [(1e1f64, 1e7f64), (1e3, 1e7), (1e5, 1e7)] {
                    let p: f64 = (n + 2 * l + 1) as f64 / 2.0;
                    let expected = ((2.0 * alpha + g2) / (2.0 * alpha + g1)).powf(p);
                    for l_ecp in [-1, l] {
                        let m1 = ecp_scalar(l, alpha, n, g1, l_ecp);
                        let m2 = ecp_scalar(l, alpha, n, g2, l_ecp);
                        let ratio = m1 / m2;
                        let err = (ratio / expected - 1.0).abs();
                        assert!(
                            err < 1e-10,
                            "l={l} alpha={alpha} n={n} g1={g1} g2={g2} l_ecp={l_ecp}: \
                             ratio={ratio} expected={expected} err={err}"
                        );
                    }
                }
            }
        }
    }
}
