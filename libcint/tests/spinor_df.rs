//! Density-fitting integrals with a spinor bra pair and a scalar auxiliary
//! center.
//!
//! `int3c2e*_spinor` evaluates
//! `(mu_spinor nu_spinor | L_sph)`, so its three centers do not share one
//! representation. These tests pin the per-center layout metadata, the cross
//! integration path that uses it, and the deliberate rejection of the ordinary
//! (non-cross) path.
//!
//! Numerical `s2ij` references are compared through `integrate_row_major`:
//! PySCF returns C-order full tensors and packs the lower triangle, which is
//! exactly the documented row-major layout. The column-major layout packs the
//! upper triangle; for the Hermitian bra pair it is the elementwise conjugate
//! of the row-major buffer.

use approx::assert_relative_eq;
use libcint::prelude::*;
use num::complex::ComplexFloat;
use num::Complex;

use CIntType::{Spheric, Spinor};

const MIXED_INTEGRATORS: [&str; 13] = [
    "int3c2e",
    "int3c2e_ip1",
    "int3c2e_ip2",
    "int3c2e_pvp1",
    "int3c2e_pvxp1",
    "int3c2e_ig1",
    "int3c2e_spsp1",
    "int3c2e_ipspsp1",
    "int3c2e_spsp1ip2",
    "int3c2e_ipip1",
    "int3c2e_ipip2",
    "int3c2e_ipvip1",
    "int3c2e_ip1ip2",
];

fn assert_fp(out: &[Complex<f64>], expected: Complex<f64>) {
    let fp = cint_fp(out);
    assert_relative_eq!(fp.re, expected.re, epsilon = 1e-10);
    assert_relative_eq!(fp.im, expected.im, epsilon = 1e-10);
}

fn h2o_spinor_and_scalar_aux() -> (CInt, CInt) {
    let mut bra = init_h2o_def2_tzvp();
    bra.cint_type = Spinor;
    let mut aux = init_h2o_def2_jk();
    aux.cint_type = Spheric;
    (bra, aux)
}

#[test]
fn mixed_layout_metadata() {
    for intor in MIXED_INTEGRATORS {
        let integrator = CInt::get_integrator(intor);
        let layout = integrator.center_layout().unwrap_or_else(|| panic!("{intor} has no layout"));
        assert_eq!(layout.driver_type, Spinor, "{intor}");
        assert_eq!(layout.center_types, [Spinor, Spinor, Spheric], "{intor}");
        assert_eq!(layout.center_types.len(), integrator.n_center(), "{intor}");
    }

    assert!(CInt::get_integrator("int2e").center_layout().is_none());
    assert!(CInt::get_integrator("int2c2e").center_layout().is_none());
}

#[test]
fn cross_spinor_accepts_default_scalar_molecules() {
    // Ordinary test molecules default to `Spheric`; the cross spinor API and
    // the preset mixed layout must still select spinor for centers 0/1.
    let bra = init_h2o_def2_tzvp();
    let aux = init_h2o_def2_jk();
    assert_eq!(bra.cint_type, Spheric);
    assert_eq!(aux.cint_type, Spheric);

    let (out, shape) = CInt::integrate_cross_spinor("int3c2e", [&bra, &bra, &aux], "s1", None).into();
    assert_eq!(shape, [86, 86, 113]);
    assert_fp(&out, Complex::new(13.409726582360673, -8.692140806901875));
}

#[test]
fn cross_ignores_argument_cint_type() {
    let (bra, aux) = h2o_spinor_and_scalar_aux();
    let mols = [&bra, &bra, &aux];

    for intor in ["int3c2e", "int3c2e_spinor", "int3c2e_sph"] {
        let (out, shape) = CInt::integrate_cross_spinor(intor, mols, "s1", None).into();
        assert_eq!(shape, [86, 86, 113], "{intor}");
        assert_fp(&out, Complex::new(13.409726582360673, -8.692140806901875));
    }
}

#[test]
fn cross_with_prebuilt_optimizer_and_row_major() {
    let (bra, aux) = h2o_spinor_and_scalar_aux();
    let mols = [&bra, &bra, &aux];

    // `get_cintopt_cross` must also accept the mixed input types.
    let opt = CInt::get_cintopt_cross("int3c2e_spinor", mols);
    let args = bra.integrate_cross_args_builder_spinor().intor("int3c2e").mols(&mols).aosym("s1").cintopt(&opt).build().unwrap();
    let (out, shape) = CInt::integrate_cross_with_args_spinor(args).into();
    assert_eq!(shape, [86, 86, 113]);
    assert_fp(&out, Complex::new(13.409726582360673, -8.692140806901875));

    let (out, shape) = CInt::integrate_cross_row_major_spinor("int3c2e", mols, "s2ij", None).into();
    assert_eq!(shape, [3741, 113]);
    assert_fp(&out, Complex::new(205.63015904614224, -1.4588346679304407));

    let (out, shape) = CInt::integrate_cross_row_major_spinor("int3c2e", mols, "s1", None).into();
    assert_eq!(shape, [86, 86, 113]);
    assert_fp(&out, Complex::new(162.91776666915572, 48.89793366217381));
}

#[test]
fn cross_row_major_numerical_references() {
    // PySCF's `int3c2e*_spinor` references use C-order full tensors and pack
    // the lower triangle for `s2ij`, which is exactly the row-major layout.
    // Generated with
    // `df.incore.aux_e2(mol, auxmol, intor + "_spinor", aosym)` and
    // `auxmol = df.addons.make_auxmol(mol, "def2-universal-jkfit")`.
    let (bra, aux) = h2o_spinor_and_scalar_aux();
    let mols = [&bra, &bra, &aux];

    let cases: [(&str, &str, &[usize], Complex<f64>); 8] = [
        ("int3c2e_spinor", "s1", &[86, 86, 113], Complex::new(162.91776666915572, 48.89793366217381)),
        ("int3c2e_spinor", "s2ij", &[3741, 113], Complex::new(205.63015904614224, -1.4588346679304407)),
        ("int3c2e_spsp1_spinor", "s1", &[86, 86, 113], Complex::new(1169.3027797396803, -1.1819567573974155)),
        ("int3c2e_spsp1_spinor", "s2ij", &[3741, 113], Complex::new(49.04933321951309, 79.63260644068524)),
        ("int3c2e_ip1_spinor", "s1", &[3, 86, 86, 113], Complex::new(26.315042258933516, 8.124342029226401)),
        ("int3c2e_ip1_spinor", "s2ij", &[3, 3741, 113], Complex::new(397.101858927927, -36.34893216730843)),
        ("int3c2e_ip2_spinor", "s1", &[3, 86, 86, 113], Complex::new(37.30509426228551, 30.567798214841673)),
        ("int3c2e_ip2_spinor", "s2ij", &[3, 3741, 113], Complex::new(12.77406409829333, -0.40258380777703895)),
    ];

    for (intor, aosym, expected_shape, expected_fp) in cases {
        let (out, shape) = CInt::integrate_cross_row_major_spinor(intor, mols, aosym, None).into();
        assert_eq!(shape, expected_shape, "{intor}/{aosym}");
        assert_fp(&out, expected_fp);
    }
}

#[test]
fn cross_s2ij_shape_and_value_row_major() {
    let (bra, aux) = h2o_spinor_and_scalar_aux();
    let (out, shape) = CInt::integrate_cross_row_major_spinor("int3c2e", [&bra, &bra, &aux], "s2ij", None).into();
    assert_eq!(shape, [3741, 113]);
    assert_fp(&out, Complex::new(205.63015904614224, -1.4588346679304407));
}

#[test]
fn cross_s2ij_col_major_is_conjugate_of_row_major() {
    // Column-major `s2ij` packs the upper triangle; the mixed 3c-2e spinor
    // tensor is Hermitian in the bra pair, so its values are the elementwise
    // conjugates of the row-major (lower-triangle, PySCF-compatible) buffer.
    let (bra, aux) = h2o_spinor_and_scalar_aux();
    let mols = [&bra, &bra, &aux];
    let (ntri, naux) = (3741usize, 113usize);

    let (col, col_shape) = CInt::integrate_cross_spinor("int3c2e", mols, "s2ij", None).into();
    let (row, row_shape) = CInt::integrate_cross_row_major_spinor("int3c2e", mols, "s2ij", None).into();
    assert_eq!(col_shape, [ntri, naux]);
    assert_eq!(row_shape, [ntri, naux]);
    for k in 0..ntri {
        for p in 0..naux {
            let c = col[k + ntri * p];
            let r = row[k * naux + p];
            assert!((c - r.conj()).norm() < 1e-10, "int3c2e s2ij k={k} p={p}: {c:?} vs {r:?}");
        }
    }

    let (col, col_shape) = CInt::integrate_cross_spinor("int3c2e_ip2", mols, "s2ij", None).into();
    let (row, row_shape) = CInt::integrate_cross_row_major_spinor("int3c2e_ip2", mols, "s2ij", None).into();
    assert_eq!(col_shape, [ntri, naux, 3]);
    assert_eq!(row_shape, [3, ntri, naux]);
    for c in 0..3 {
        for k in 0..ntri {
            for p in 0..naux {
                let cv = col[k + ntri * (p + naux * c)];
                let rv = row[c * ntri * naux + k * naux + p];
                assert!((cv - rv.conj()).norm() < 1e-10, "int3c2e_ip2 s2ij c={c} k={k} p={p}: {cv:?} vs {rv:?}");
            }
        }
    }
}

#[test]
fn cross_overlapping_centers_use_distinct_center_types() {
    let mut bra = init_h2o_def2_tzvp();
    bra.cint_type = Spinor;
    let mols = [&bra, &bra, &bra];

    let cases = [
        ("int3c2e_spinor", Complex::new(-12.697280181635538, -12.540288220212021)),
        ("int3c2e_spsp1_spinor", Complex::new(-1173.0182942325184, -52.05443849137967)),
    ];

    for (intor, expected_fp) in cases {
        let (out, shape) = CInt::integrate_cross_spinor(intor, mols, "s1", None).into();
        assert_eq!(shape, [86, 86, 43], "{intor}");
        assert_fp(&out, expected_fp);
    }
}

#[test]
fn cross_shapes_for_all_variants() {
    let (bra, aux) = h2o_spinor_and_scalar_aux();
    let mols = [&bra, &bra, &aux];

    for intor in MIXED_INTEGRATORS {
        let n_comp = CInt::get_integrator(intor).n_spinor_comp();
        let expected_s1 = if n_comp > 1 { vec![86, 86, 113, n_comp] } else { vec![86, 86, 113] };
        let expected_s2ij = if n_comp > 1 { vec![3741, 113, n_comp] } else { vec![3741, 113] };

        let (_, shape) = CInt::integrate_cross_spinor(intor, mols, "s1", None).into();
        assert_eq!(shape, expected_s1, "{intor}/s1");

        let (_, shape) = CInt::integrate_cross_spinor(intor, mols, "s2ij", None).into();
        assert_eq!(shape, expected_s2ij, "{intor}/s2ij");
    }
}

#[test]
fn non_cross_mixed_layout_is_rejected() {
    let mut bra = init_h2o_def2_tzvp();
    bra.cint_type = Spinor;

    for result in [
        bra.integrate_spinor_f("int3c2e_spinor", "s1", None).map(|out| out.shape),
        bra.integrate_row_major_spinor_f("int3c2e_spinor", "s1", None).map(|out| out.shape),
        bra.integrate_f("int3c2e", "s1", None).map(|out| out.shape),
    ] {
        let err = result.expect_err("mixed-layout non-cross integration must fail");
        assert!(matches!(err.kind, CIntErrorBase::IntegratorNotAvailable(_)), "{err:?}");
    }

    // Scalar variants of the same family must keep working.
    let (_, shape) = bra.integrate("int3c2e_sph", "s1", None).into();
    assert_eq!(shape, [43, 43, 43]);
    let (_, shape) = bra.integrate_row_major("int3c2e_sph", "s1", None).into();
    assert_eq!(shape, [43, 43, 43]);
    let (_, shape) = bra.integrate("int3c2e_cart", "s1", None).into();
    assert_eq!(shape, [48, 48, 48]);
}

#[test]
fn nonzero_kappa_keeps_scalar_aux_size() {
    let mut bra = init_h2o_def2_tzvp();
    bra.cint_type = Spinor;
    // KAPPA_OF == 4; a p shell with kappa = -2 has four spinors.
    bra.bas[5][4] = -2;
    assert_eq!(bra.nao(), 84);

    let (out, shape) = CInt::integrate_cross_spinor("int3c2e", [&bra, &bra, &bra], "s1", None).into();
    assert_eq!(shape, [84, 84, 43]);
    assert_fp(&out, Complex::new(90.85924581862466, 8.893668202986625));
}
