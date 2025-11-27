#[cfg(test)]
mod test {
    use libcint::{gto::gto_crafter::get_gto_eval_name_f, prelude::*};
    use num::Complex;
    use rstest::rstest;

    #[rstest]
    #[case("GTOval_sph"                        ,     378.9344483361458 , [2048, 1464,  1])]
    #[case("GTOval_cart"                       ,     198.05403491229913, [2048, 1774,  1])]
    #[case("GTOval_sph_deriv1"                 ,    1158.9931602933664 , [2048, 1464,  4])]
    #[case("GTOval_cart_deriv1"                ,   -1725.5519928866993 , [2048, 1774,  4])]
    #[case("GTOval_sph_deriv2"                 ,    2777.1107970453154 , [2048, 1464, 10])]
    #[case("GTOval_cart_deriv2"                ,     718.8381138263526 , [2048, 1774, 10])]
    #[case("GTOval_sph_deriv3"                 ,    5879.015564180246  , [2048, 1464, 20])]
    #[case("GTOval_cart_deriv3"                ,   -2664.5459224494057 , [2048, 1774, 20])]
    #[case("GTOval_sph_deriv4"                 ,   30771.778920777477  , [2048, 1464, 35])]
    #[case("GTOval_cart_deriv4"                ,   33427.10951556449   , [2048, 1774, 35])]
    #[case("GTOval_ip_sph"                     ,     331.80038695474656, [2048, 1464,  3])]
    #[case("GTOval_ig_sph"                     ,    -751.1976925410391 , [2048, 1464,  3])]
    #[case("GTOval_ipig_sph"                   ,   -3052.9189303084386 , [2048, 1464,  9])]
    #[case("GTOval_ipr_sph"                    ,      54.16700105790881, [2048, 1464,  9])]
    #[case("GTOval_iprc_sph"                   ,    -901.1885030079377 , [2048, 1464,  9])]
    fn test_usual_case(#[case] eval_name: &str, #[case] ref_fp: f64, #[case] ref_shape: impl AsRef<[usize]>) {
        // this will test usual case (naive `eval_gto` call)
        let mol = init_c10h22_def2_qzvp();
        let ngrid = 2048;
        let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
        let (out, shape) = mol.eval_gto(eval_name, &coord).into();
        let out_fp = cint_fp(&out);
        assert_eq!(shape, ref_shape.as_ref());
        println!("out_fp = {}", out_fp);
        assert!((out_fp - ref_fp).abs() < 1e-10 || (out_fp / ref_fp - 1.0).abs() < 1e-7);
    }

    #[rstest]
    #[case("GTOval_sph"                        , 1e-8 ,    -89.26220691313141, [2048, 1358,  1])]
    #[case("GTOval_sph"                        , 1e-2 ,    -86.33740170264781, [2048, 1358,  1])]
    #[case("GTOval_cart"                       , 1e-5 ,    143.0683173879014 , [2048, 1649,  1])]
    #[case("GTOval_sph_deriv1"                 , 1e-2 ,   -926.9693102995727 , [2048, 1358,  4])]
    #[case("GTOval_cart_deriv1"                , 1e-2 ,  -2165.767335096756  , [2048, 1649,  4])]
    #[case("GTOval_sph_deriv2"                 , 1e-2 ,   -427.9248897536187 , [2048, 1358, 10])]
    #[case("GTOval_cart_deriv2"                , 1e-2 ,  -2661.780110374988  , [2048, 1649, 10])]
    #[case("GTOval_sph_deriv3"                 , 1e-2 ,  -4097.951202338538  , [2048, 1358, 20])]
    #[case("GTOval_cart_deriv3"                , 1e-2 ,  -2321.4395429846345 , [2048, 1649, 20])]
    #[case("GTOval_sph_deriv4"                 , 1e-2 ,   -974.7870302755828 , [2048, 1358, 35])]
    #[case("GTOval_cart_deriv4"                , 1e-2 , -31806.509015190735  , [2048, 1649, 35])]
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
        assert!((out_fp - ref_fp).abs() < 1e-10 || (out_fp / ref_fp - 1.0).abs() < 1e-7);
    }

    #[rstest]
    #[case("GTOval_spinor"                        ,    -125.69056292477536,    -993.648317458363  , [2048, 2928,  1,  2])]
    #[case("GTOval_spinor_deriv1"                 ,    -550.601092411041  ,     372.2702812731907 , [2048, 2928,  4,  2])]
    #[case("GTOval_spinor_deriv2"                 ,   -1018.4141009661886 ,   -1747.377520581361  , [2048, 2928, 10,  2])]
    #[case("GTOval_spinor_deriv3"                 ,   -6868.453294385325  ,      39.34436040140736, [2048, 2928, 20,  2])]
    #[case("GTOval_spinor_deriv4"                 ,   -17119.668571065467 ,   -4260.1550424077595 , [2048, 2928, 35,  2])]
    #[case("GTOval_spinor_ip"                     ,     1164.5992226651474,   -2165.5425060625184 , [2048, 2928,  3,  2])]
    #[case("GTOval_spinor_ig"                     ,      537.221323651945 ,   -1539.282115487602  , [2048, 2928,  3,  2])]
    #[case("GTOval_spinor_ipig"                   ,      335.2086997959833,    3726.8878032087882 , [2048, 2928,  9,  2])]
    #[case("GTOval_spinor_ipr"                    ,    -3162.527253858847 ,    1177.7605640921986 , [2048, 2928,  9,  2])]
    #[case("GTOval_spinor_iprc"                   ,      753.689035698545 ,   -1395.6027721906285 , [2048, 2928,  9,  2])]
    #[case("GTOval_spinor_sp"                     ,     -276.2710005825974,    -689.2353028137179 , [2048, 2928,  1,  2])]
    #[case("GTOval_spinor_ipsp"                   ,     2829.8957257389334,   -5169.514008491763  , [2048, 2928,  3,  2])]
    #[case("GTOval_spinor_ipipsp"                 ,     3318.597980963531 ,  -14499.595354496942  , [2048, 2928,  9,  2])]
    fn test_usual_case_spinor(
        #[case] eval_name: &str,
        #[case] ref_fp_re: f64,
        #[case] ref_fp_im: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        // this will test usual case (naive `eval_gto` call)
        let mol = init_c10h22_def2_qzvp();
        let ngrid = 2048;
        let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
        let (out, shape) = mol.eval_gto_spinor(eval_name, &coord).into();
        let out_fp = cint_fp(&out);
        assert_eq!(shape, ref_shape.as_ref());
        println!("out_fp = {}", out_fp);
        let ref_fp = Complex::<f64>::new(ref_fp_re, ref_fp_im);
        assert!((out_fp - ref_fp).norm() < 1e-10 || (out_fp / ref_fp - 1.0).norm() < 1e-7);
    }
}
