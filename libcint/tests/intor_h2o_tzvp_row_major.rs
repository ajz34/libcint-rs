use approx::assert_relative_eq;
use libcint::prelude::*;
use rstest::rstest;

#[cfg(test)]
mod test {
    use super::*;

    #[rstest]
    #[case("int2c2e"        , "s1"  , None                , 353.92327889575165    , [43, 43]      )] // no-deriv
    #[case("int1e_igovlp"   , "s1"  , None                , 3.357208986844009     , [3, 43, 43]   )] // deriv
    #[case("int1e_kin"      , "s1"  , [[0, 12], [8, 18]]  , 12.899828616445948    , [32, 26]      )] // no-deriv | shl
    #[case("int1e_ipnuc"    , "s1"  , [[0, 12], [8, 18]]  , -27.63666200750539    , [3, 32, 26]   )] // deriv | shl
    fn test_2c_sph(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl Into<ShlsSlice>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let cint_data = init_h2o_def2_tzvp();

        let (out, shape) = cint_data.integrate_row_major(intor, aosym, shls_slice).into();
        assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-12);
        assert_eq!(shape, ref_shape.as_ref());
    }

    #[rstest]
    #[case("int3c2e"        , "s1"  , None                         , 48.161159148027394    , [43, 43, 43]      )] // no-deriv
    #[case("int3c2e_ig1"    , "s1"  , None                         , 20.278756457464894    , [3, 43, 43, 43]   )] // deriv
    #[case("int3c2e"        , "s1"  , [[0, 12], [8, 18], [6, 15]]  , -25.296528540125045   , [32, 26, 29]      )] // no-deriv | shl
    #[case("int3c2e_ipip1"  , "s1"  , [[0, 12], [8, 18], [6, 15]]  , 968.5261168299187     , [9, 32, 26, 29]   )] // deriv | shl
    fn test_3c_sph(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl Into<ShlsSlice>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let cint_data = init_h2o_def2_tzvp();

        let (out, shape) = cint_data.integrate_row_major(intor, aosym, shls_slice).into();
        assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-12);
        assert_eq!(shape, ref_shape.as_ref());
    }

    #[rstest]
    #[case("int2e"              , "s1"  , None                                  , 70.00106603114827     , [43, 43, 43, 43]      )] // no-deriv
    #[case("int2e_ip1"          , "s1"  , None                                  , 26.41216164838567     , [3, 43, 43, 43, 43]   )] // deriv
    #[case("int2e"              , "s1"  , [[0, 12], [8, 18], [6, 15], [3, 10]]  , 7.039814725057187     , [32, 26, 29, 21]      )] // no-deriv | shl
    #[case("int2e_cg_ssa10ssp2" , "s1"  , [[0, 12], [8, 18], [6, 15], [3, 10]]  , 36.40584837312814     , [48, 32, 26, 29, 21]  )] // deriv | shl
    fn test_4c_sph(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl Into<ShlsSlice>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let cint_data = init_h2o_def2_tzvp();

        let (out, shape) = cint_data.integrate_row_major(intor, aosym, shls_slice).into();
        assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-11);
        assert_eq!(shape, ref_shape.as_ref());
    }
}
