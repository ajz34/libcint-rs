#[cfg(test)]
mod test {
    use libcint::{gto::gto_crafter::get_gto_eval_name_f, prelude::*};
    use rstest::rstest;

    #[rstest]
    #[case("GTOval_sph"                        ,    378.9344483361458 , [2048, 1464,  1])]
    #[case("GTOval_cart"                       ,    198.05403491229913, [2048, 1774,  1])]
    #[case("GTOval_sph_deriv1"                 ,   1158.9931602933664 , [2048, 1464,  4])]
    #[case("GTOval_cart_deriv1"                ,  -1725.5519928866993 , [2048, 1774,  4])]
    #[case("GTOval_sph_deriv2"                 ,   2777.1107970453154 , [2048, 1464, 10])]
    #[case("GTOval_cart_deriv2"                ,    718.8381138263526 , [2048, 1774, 10])]
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

    #[rstest]
    #[case("GTOval_sph"                        , 1e-8 ,   -89.26220691313141, [2048, 1358,  1])]
    #[case("GTOval_sph"                        , 1e-2 ,   -86.33740170264781, [2048, 1358,  1])]
    #[case("GTOval_cart"                       , 1e-5 ,   143.0683173879014 , [2048, 1649,  1])]
    #[case("GTOval_sph_deriv1"                 , 1e-2 ,  -926.9693102995727 , [2048, 1358,  4])]
    #[case("GTOval_cart_deriv1"                , 1e-2 , -2165.767335096756  , [2048, 1649,  4])]
    #[case("GTOval_sph_deriv2"                 , 1e-2 ,  -427.9248897536187 , [2048, 1358, 10])]
    #[case("GTOval_cart_deriv2"                , 1e-2 , -2661.780110374988  , [2048, 1649, 10])]
    fn test_with_args(#[case] eval_name: &str, #[case] cutoff: f64, #[case] ref_fp: f64, #[case] ref_shape: impl AsRef<[usize]>) {
        let mol = init_c10h22_def2_qzvp();
        let ngrid = 2048;
        let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();

        let shls_slice = [4, 432];
        let (evaluator, cint_type) = get_gto_eval_name_f(eval_name).unwrap();
        let ao_loc = mol.make_loc_with_type(cint_type.unwrap_or(mol.cint_type));
        let nao = ao_loc[shls_slice[1]] - ao_loc[shls_slice[0]];
        let ncomp = evaluator.ncomp();
        let mut out = vec![f64::NAN; ngrid * nao * ncomp];
        let args =
            mol.gto_args_builder().eval_name(eval_name).coord(&coord).shls_slice(shls_slice).out(&mut out).cutoff(cutoff).build().unwrap();
        let result = mol.eval_gto_with_args(args);
        let out_fp = cint_fp(&out);
        assert_eq!(result.shape, ref_shape.as_ref());
        assert!((out_fp - ref_fp).abs() < 1e-10);
    }
}
