use libcint::prelude::*;

#[test]
fn test_workable() {
    let cint_data = init_h2o_def2_tzvp();
    let (out, shape) = cint_data.integrate("int1e_ovlp", None, None).into();
    assert!((cint_fingerprint(&out) - 33.02101079126514).abs() < 1e-10);
    assert_eq!(shape, [43, 43]);
}
