#[cfg(test)]
mod test {
    use libcint::prelude::*;
    use rstest::rstest;

    #[rstest]
    #[case("GTOval_sph"                   ,   378.9344483361458 , [2048, 1464, 1])]
    #[case("GTOval_cart"                  ,   198.05403491229913, [2048, 1774, 1])]
    #[case("GTOval_sph_deriv1"            ,  1158.9931602933664 , [2048, 1464, 4])]
    #[case("GTOval_cart_deriv1"           , -1725.5519928866993 , [2048, 1774, 4])]
    fn test_usual_case(#[case] eval_name: &str, #[case] ref_fp: f64, #[case] ref_shape: impl AsRef<[usize]>) {
        // this will test usual case (naive `eval_gto` call)
        let mol = init_c10h22_def2_qzvp();
        let ngrid = 2048;
        let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
        let (out, shape) = mol.eval_gto(eval_name, &coord).into();
        let out_fp = cint_fp(&out);
        assert_eq!(shape, ref_shape.as_ref());
        assert!((out_fp - ref_fp).abs() < 1e-10);
    }
}
