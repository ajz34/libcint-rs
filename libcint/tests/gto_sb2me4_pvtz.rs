#[cfg(test)]
mod test {
    use libcint::{gto::gto_crafter::get_gto_eval_name_f, prelude::*};
    use num::Complex;
    use rstest::rstest;

    #[rstest]
    #[case("GTOval_sph"                        ,    -291.50611912900354, [2048,  366,  1])]
    #[case("GTOval_cart"                       ,    -109.46020421466841, [2048,  410,  1])]
    #[case("GTOval_ig_sph"                     ,     180.75613409086753, [2048,  366,  3])]
    fn test_usual_case(#[case] eval_name: &str, #[case] ref_fp: f64, #[case] ref_shape: impl AsRef<[usize]>) {
        // this will test usual case (naive `eval_gto` call)
        let mol = init_sb2me4_cc_pvtz();
        let ngrid = 2048;
        let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
        let (out, shape) = mol.eval_gto(eval_name, &coord).into();
        let out_fp = cint_fp(&out);
        assert_eq!(shape, ref_shape.as_ref());
        println!("out_fp = {}", out_fp);
        assert!((out_fp - ref_fp).abs() < 1e-10 || (out_fp / ref_fp - 1.0).abs() < 1e-7);
    }

    #[rstest]
    #[case("GTOval_sph"                        , 1e-4 ,     20.06329070736091, [2048,  349,  1])]
    #[case("GTOval_sph_deriv1"                 , 1e-2 ,    228.2644447105396 , [2048,  349,  4])]
    #[case("GTOval_cart_deriv1"                , 1e-2 ,   -427.8204853612391 , [2048,  392,  4])]
    fn test_with_args(#[case] eval_name: &str, #[case] cutoff: f64, #[case] ref_fp: f64, #[case] ref_shape: impl AsRef<[usize]>) {
        let mol = init_sb2me4_cc_pvtz();
        let ngrid = 2048;
        let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();

        let shls_slice = [4, 126];
        let (evaluator, cint_type) = get_gto_eval_name_f(eval_name).unwrap();
        let ao_loc = mol.make_loc_with_type(cint_type.unwrap_or(mol.cint_type));
        let nao = ao_loc[shls_slice[1]] - ao_loc[shls_slice[0]];
        let ncomp = evaluator.ncomp();
        let mut out = vec![f64::NAN; ngrid * nao * ncomp];
        let args =
            mol.gto_args_builder().eval_name(eval_name).coord(&coord).shls_slice(shls_slice).out(&mut out).cutoff(cutoff).build().unwrap();
        let result = mol.eval_gto_with_args(args);
        let out_fp = cint_fp(&out);
        println!("out_fp = {}", out_fp);
        assert_eq!(result.shape, ref_shape.as_ref());
        assert!((out_fp - ref_fp).abs() < 1e-10 || (out_fp / ref_fp - 1.0).abs() < 1e-7);
    }

    #[rstest]
    #[case("GTOval_spinor"                        ,    -198.61510549238164,     257.7618955783744 , [2048,  732,  1,  2])]
    #[case("GTOval_spinor_ig"                     ,      26.57165093832205,     128.3701183295896 , [2048,  732,  3,  2])]
    #[case("GTOval_spinor_sp"                     ,     -41.42811657970221,     -13.38432159407   , [2048,  732,  1,  2])]
    fn test_usual_case_spinor(
        #[case] eval_name: &str,
        #[case] ref_fp_re: f64,
        #[case] ref_fp_im: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        // this will test usual case (naive `eval_gto` call)
        let mol = init_sb2me4_cc_pvtz();
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
