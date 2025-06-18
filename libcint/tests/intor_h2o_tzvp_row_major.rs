// Reference source: libcint/scripts/pyscf_h2o_tzvp_row_major.ipynb

use approx::assert_relative_eq;
use libcint::prelude::*;
use rstest::rstest;

#[cfg(test)]
mod test {
    use super::*;

    #[rstest]
    // s1
    #[case("int2c2e"        , "s1"  , None                , 353.92327889575165    , [43, 43]      )] // no-deriv
    #[case("int1e_igovlp"   , "s1"  , None                , 3.357208986844009     , [3, 43, 43]   )] // deriv
    #[case("int1e_kin"      , "s1"  , [[0, 12], [8, 18]]  , 12.899828616445948    , [32, 26]      )] // no-deriv | shl
    #[case("int1e_ipnuc"    , "s1"  , [[0, 12], [8, 18]]  , -27.63666200750539    , [3, 32, 26]   )] // deriv | shl
    // s2ij
    #[case("int2c2e"        , "s2ij", None                , -176.5328705658656    , [946]         )] // no-deriv
    #[case("int1e_igovlp"   , "s2ij", None                , 1.3748304523870443    , [3, 946]      )] // deriv
    #[case("int1e_kin"      , "s2ij", [[8, 18], [8, 18]]  , -14.04470523493559    , [351]         )] // no-deriv | shl
    #[case("int1e_ipnuc"    , "s2ij", [[8, 18], [8, 18]]  , 23.32391114542229     , [3, 351]      )] // deriv | shl
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
    // s1
    #[case("int3c2e"        , "s1"  , None                         , 48.161159148027394    , [43, 43, 43]      )] // no-deriv
    #[case("int3c2e_ig1"    , "s1"  , None                         , 20.278756457464894    , [3, 43, 43, 43]   )] // deriv
    #[case("int3c2e"        , "s1"  , [[0, 12], [8, 18], [6, 15]]  , -25.296528540125045   , [32, 26, 29]      )] // no-deriv | shl
    #[case("int3c2e_ipip1"  , "s1"  , [[0, 12], [8, 18], [6, 15]]  , 968.5261168299187     , [9, 32, 26, 29]   )] // deriv | shl
    // s2ij
    #[case("int3c2e"        , "s2ij", None                         , -13.182225570003517   , [946, 43]         )] // no-deriv
    #[case("int3c2e_ip2"    , "s2ij", None                         , -1.656437965487193    , [3, 946, 43]      )] // deriv
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
    // s1
    #[case("int2e"              , "s1"  , None                                  , 70.00106603114827     , [43, 43, 43, 43]      )] // no-deriv
    #[case("int2e_ip1"          , "s1"  , None                                  , 26.41216164838567     , [3, 43, 43, 43, 43]   )] // deriv
    #[case("int2e"              , "s1"  , [[0, 12], [8, 18], [6, 15], [3, 10]]  , 7.039814725057187     , [32, 26, 29, 21]      )] // no-deriv | shl
    #[case("int2e_cg_ssa10ssp2" , "s1"  , [[0, 12], [8, 18], [6, 15], [3, 10]]  , 36.40584837312814     , [48, 32, 26, 29, 21]  )] // deriv | shl
    // s2ij
    #[case("int2e"              , "s2ij", None                                  , -135.4209188844083    , [946, 43, 43]         )] // no-deriv
    #[case("int2e_gg1"          , "s2ij", None                                  , 148.48845559162677    , [9, 946, 43, 43]      )] // deriv
    #[case("int2e"              , "s2ij", [[8, 18], [8, 18], [6, 15], [3, 10]]  , -6.438967844331271    , [351, 29, 21]         )] // no-deriv | shl
    #[case("int2e_gg1"          , "s2ij", [[8, 18], [8, 18], [6, 15], [3, 10]]  , -5.312151890309538    , [9, 351, 29, 21]      )] // deriv | shl
    // s2kl
    #[case("int2e"              , "s2kl", None                                  , 1.8306651475960938    , [43, 43, 946]         )] // no-deriv
    #[case("int2e_ip1"          , "s2kl", None                                  , -37.34064823381276    , [3, 43, 43, 946]      )] // deriv
    #[case("int2e"              , "s2kl", [[0, 12], [8, 18], [6, 15], [6, 15]]  , 5.793369823495306     , [32, 26, 435]         )] // no-deriv | shl
    #[case("int2e_ip1"          , "s2kl", [[0, 12], [8, 18], [6, 15], [6, 15]]  , -2.3043335457616916   , [3, 32, 26, 435]      )] // deriv | shl
    // s4
    #[case("int2e"              , "s4"  , None                                  , -2.1937308808534586   , [946, 946]            )] // no-deriv
    #[case("int2e"              , "s4"  , [[0, 12], [0, 12], [6, 15], [6, 15]]  , -10.689968606671224   , [528, 435]            )] // no-deriv | shl
    // s8
    #[case("int2e"              , "s8"  , None                                  , 21.90951424661062     , [447931]              )] // no-deriv
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
