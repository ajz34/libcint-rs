//! `int3c2e*_ssc` integrals: scalar/spherical bra centers against a Cartesian
//! auxiliary center.
//!
//! The `_ssc` entry points are exported by libcint but omitted from its
//! `cint_funcs.h`; qcint does not provide them publicly. These tests pin the
//! scalar cross path and the per-center layout metadata.

// qcint does not export the `_ssc` entry points, so the integrator entries are
// compiled out under that feature; skip the whole test crate there.
#![cfg(not(feature = "qcint"))]

use approx::assert_relative_eq;
use libcint::prelude::*;

use CIntType::{Cartesian, Spheric, Spinor};

#[test]
fn ssc_layout_metadata() {
    let scalar = CInt::get_integrator("int3c2e_sph_ssc");
    let layout = scalar.center_layout().unwrap();
    assert_eq!(layout.driver_type, Spheric);
    assert_eq!(layout.center_types, [Spheric, Spheric, Cartesian]);
    assert!(scalar.is_sph_available());
    assert!(!scalar.is_cart_available());
    assert!(!scalar.is_spinor_available());

    for intor in ["int3c2e_spinor_ssc", "int3c2e_spsp1_spinor_ssc"] {
        let integrator = CInt::get_integrator(intor);
        let layout = integrator.center_layout().unwrap();
        assert_eq!(layout.driver_type, Spinor, "{intor}");
        assert_eq!(layout.center_types, [Spinor, Spinor, Cartesian], "{intor}");
        assert!(!integrator.is_sph_available(), "{intor}");
        assert!(!integrator.is_cart_available(), "{intor}");
        assert!(integrator.is_spinor_available(), "{intor}");
    }
    assert_eq!(CInt::get_integrator("int3c2e_sph_ssc").n_comp(), 1);
    assert_eq!(CInt::get_integrator("int3c2e_spsp1_spinor_ssc").n_comp(), 4);
    assert_eq!(CInt::get_integrator("int3c2e_spsp1_spinor_ssc").n_spinor_comp(), 1);
}

#[test]
fn scalar_ssc_cross_numerical_references() {
    let bra = init_h2o_def2_tzvp();
    let aux = init_h2o_def2_jk();
    assert_eq!(bra.cint_type, Spheric);
    assert_eq!(aux.cint_type, Spheric);

    let mols = [&bra, &bra, &aux];
    // PySCF's generic `getints3c` scalar `_ssc` s1 fill is inconsistent, so the
    // s1 references below are reconstructed by symmetrizing its (correct) s2ij
    // tensor. Both use the row-major C-order layout, which matches PySCF's
    // packed s2ij order.
    let (out, shape) = CInt::integrate_cross_row_major("int3c2e_sph_ssc", mols, "s1", None).into();
    assert_eq!(shape, [43, 43, 133]);
    assert_relative_eq!(cint_fp(&out), 40.46812503253784, epsilon = 1e-10);

    let (out, shape) = CInt::integrate_cross_row_major("int3c2e_sph_ssc", mols, "s2ij", None).into();
    assert_eq!(shape, [946, 133]);
    assert_relative_eq!(cint_fp(&out), -243.99421697179832, epsilon = 1e-10);

    // The preset per-center layout wins over the auxiliary molecule's own
    // `cint_type`; only the layout's Cartesian third center matters.
    let mut aux_cart = init_h2o_def2_jk();
    aux_cart.cint_type = Cartesian;
    let (out, shape) = CInt::integrate_cross_row_major("int3c2e_sph_ssc", [&bra, &bra, &aux_cart], "s1", None).into();
    assert_eq!(shape, [43, 43, 133]);
    assert_relative_eq!(cint_fp(&out), 40.46812503253784, epsilon = 1e-10);

    // Same molecule as the auxiliary basis: the auxiliary dimension is its
    // Cartesian AO count (48), not its spherical count (43).
    let (out, shape) = CInt::integrate_cross_row_major("int3c2e_sph_ssc", [&bra, &bra, &bra], "s1", None).into();
    assert_eq!(shape, [43, 43, 48]);
    assert_relative_eq!(cint_fp(&out), -15.314640888962277, epsilon = 1e-10);

    let (out, shape) = CInt::integrate_cross_row_major("int3c2e_sph_ssc", [&bra, &bra, &bra], "s2ij", None).into();
    assert_eq!(shape, [946, 48]);
    assert_relative_eq!(cint_fp(&out), 142.4267522201795, epsilon = 1e-10);
}

#[test]
fn scalar_ssc_optimizer_reuse() {
    let bra = init_h2o_def2_tzvp();
    let aux = init_h2o_def2_jk();
    let mols = [&bra, &bra, &aux];

    let opt = CInt::get_cintopt_cross("int3c2e_sph_ssc", mols);
    assert_eq!(opt.intor(), "int3c2e_sph_ssc");
    let args = bra.integrate_cross_args_builder().intor("int3c2e_sph_ssc").mols(&mols).aosym("s1").cintopt(&opt).build().unwrap();
    let (out, shape) = CInt::integrate_cross_with_args(args).into();
    assert_eq!(shape, [43, 43, 133]);
    assert_relative_eq!(cint_fp(&out), -413.59493458295447, epsilon = 1e-10);
}

#[test]
fn scalar_ssc_non_cross_is_rejected() {
    let bra = init_h2o_def2_tzvp();
    for result in [
        bra.integrate_f("int3c2e_sph_ssc", "s1", None).map(|out| out.shape),
        bra.integrate_row_major_f("int3c2e_sph_ssc", "s1", None).map(|out| out.shape),
    ] {
        let err = result.expect_err("mixed-layout non-cross integration must fail");
        assert!(matches!(err.kind, CIntErrorBase::IntegratorNotAvailable(_)), "{err:?}");
    }
}

#[test]
fn spinor_ssc_cross() {
    let bra = init_h2o_def2_tzvp();
    let aux = init_h2o_def2_jk();
    let mols = [&bra, &bra, &aux];

    let (out, shape) = CInt::integrate_cross_spinor("int3c2e_spinor_ssc", mols, "s1", None).into();
    assert_eq!(shape, [86, 86, 133]);
    let fp = cint_fp(&out);
    assert_relative_eq!(fp.re, 4.262858345776973, epsilon = 1e-10);
    assert_relative_eq!(fp.im, -39.15055162751742, epsilon = 1e-10);

    // PySCF references for `_ssc` use C-order (row-major) and pack the lower
    // triangle of `s2ij`, so compare the row-major output directly. Generated
    // with `moleintor.getints3c("int3c2e_spinor_ssc", atm, bas, env,
    // shls_slice, aosym="s2ij")` on the concatenated spinor/Cartesian-aux
    // environment (`auxmol.cart = True`).
    let (out, shape) = CInt::integrate_cross_row_major_spinor("int3c2e_spinor_ssc", mols, "s2ij", None).into();
    assert_eq!(shape, [3741, 133]);
    let fp = cint_fp(&out);
    assert_relative_eq!(fp.re, -131.1988242696546, epsilon = 1e-10);
    assert_relative_eq!(fp.im, 2.252888537377013, epsilon = 1e-10);

    let (_, shape) = CInt::integrate_cross_spinor("int3c2e_spsp1_spinor_ssc", mols, "s1", None).into();
    assert_eq!(shape, [86, 86, 133]);
}
