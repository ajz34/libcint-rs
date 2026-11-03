use approx::assert_relative_eq;
use libcint::prelude::*;
use num::Complex;

use CIntType::Spinor;

/// A single spinor shell's overlap with itself must be the identity, regardless
/// of `kappa`.
fn assert_identity(out: &[Complex<f64>], n: usize) {
    assert_eq!(out.len(), n * n);
    for i in 0..n {
        for j in 0..n {
            let expected = if i == j { 1.0 } else { 0.0 };
            assert_relative_eq!(out[i * n + j].re, expected, epsilon = 1e-12);
            assert_relative_eq!(out[i * n + j].im, 0.0, epsilon = 1e-12);
        }
    }
}

/// Replace the `kappa` of one shell and check all spinor-size dependent
/// quantities, then integrate that shell against itself.
fn check_nonzero_kappa(shl: usize, kappa: i32, expected_shl_nao: usize) {
    let mut data = init_h2o_def2_tzvp();
    let base_loc = data.make_loc_with_type(Spinor);
    let base_nao = *base_loc.last().unwrap();
    let base_shl_nao = base_loc[shl + 1] - base_loc[shl];

    // KAPPA_OF == 4
    data.bas[shl][4] = kappa;

    let loc = data.make_loc_with_type(Spinor);
    assert_eq!(loc[shl + 1] - loc[shl], expected_shl_nao);
    assert_eq!(data.nao_with_type(Spinor), base_nao - base_shl_nao + expected_shl_nao);

    let (out, shape) = data.integrate_row_major_spinor("int1e_ovlp_spinor", None, [[shl, shl + 1], [shl, shl + 1]]).into();
    assert_eq!(shape, vec![expected_shl_nao, expected_shl_nao]);
    assert_identity(&out, expected_shl_nao);
}

#[test]
fn len_spinor_matches_libcint_convention() {
    // kappa == 0: 4l + 2
    assert_eq!(CInt::len_spinor(0, 0), 2);
    assert_eq!(CInt::len_spinor(1, 0), 6);
    assert_eq!(CInt::len_spinor(2, 0), 10);

    // kappa < 0: 2l + 2
    assert_eq!(CInt::len_spinor(0, -1), 2);
    assert_eq!(CInt::len_spinor(1, -2), 4);
    assert_eq!(CInt::len_spinor(2, -3), 6);

    // kappa > 0: 2l
    assert_eq!(CInt::len_spinor(1, 1), 2);
    assert_eq!(CInt::len_spinor(2, 2), 4);
}

#[test]
fn nonzero_kappa_shell_lengths_and_integrals() {
    // Reported case: kappa = -1 s shell (j = 1/2) has 2 spinors, not 0.
    check_nonzero_kappa(0, -1, 2);
    // kappa < 0 for a p shell: j = 3/2 has 4 spinors.
    check_nonzero_kappa(5, -2, 4);
    // kappa > 0 for a p shell: j = 1/2 has 2 spinors.
    check_nonzero_kappa(5, 1, 2);
}
