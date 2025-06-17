use approx::assert_relative_eq;
use libcint::prelude::*;
use num::complex::ComplexFloat;
use rstest::rstest;

use CIntType::{Cartesian, Spheric, Spinor};

#[cfg(test)]
mod test {
    use super::*;

    /* #region 2-center */

    #[rstest]
    #[case("int2c2e"        , "s1"  , None                , 353.92327889575165    , [43, 43]      )] // no-deriv
    #[case("int1e_igovlp"   , "s1"  , None                , -3.3572089868440105   , [43, 43, 3]   )] // deriv
    #[case("int1e_kin"      , "s1"  , [[0, 12], [8, 18]]  , 0.8338209486794069    , [32, 26]      )] // no-deriv | shl
    #[case("int1e_ipnuc"    , "s1"  , [[0, 12], [8, 18]]  , 17.2214370364167      , [32, 26, 3]   )] // deriv | shl
    #[case("int2c2e"        , "s2ij", None                , -176.53287056586566   , [946]         )] // no-deriv | s2ij
    #[case("int1e_igovlp"   , "s2ij", None                , -1.374830452387045    , [946, 3]      )] // deriv | s2ij
    #[case("int1e_kin"      , "s2ij", [[8, 18], [8, 18]]  , -14.044705234935593   , [351]         )] // no-deriv | shl | s2ij
    #[case("int1e_ipnuc"    , "s2ij", [[8, 18], [8, 18]]  , -22.156443966726883   , [351, 3]      )] // deriv | shl | s2ij
    fn test_2c_sph(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl Into<ShlsSlice>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_h2o_def2_tzvp();

        cint_data.with_cint_type(Spheric, |cint_data| {
            let (out, shape) = cint_data.integrate(intor, aosym, shls_slice).into();
            assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-12);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    #[rstest]
    #[case("int2c2e"        , "s1"  , None                , 149.25926633235045    , [48, 48]      )] // no-deriv
    #[case("int1e_igovlp"   , "s1"  , None                , -0.7210509131585712   , [48, 48, 3]   )] // deriv
    #[case("int1e_kin"      , "s1"  , [[0, 12], [8, 18]]  , -33.34383660289398    , [37, 31]      )] // no-deriv | shl
    #[case("int1e_ipnuc"    , "s1"  , [[0, 12], [8, 18]]  , 19.297701826145243    , [37, 31, 3]   )] // deriv | shl
    #[case("int2c2e"        , "s2ij", None                , -329.69400031638395   , [1176]        )] // no-deriv | s2ij
    #[case("int1e_igovlp"   , "s2ij", None                , -0.775053635406958    , [1176, 3]     )] // deriv | s2ij
    #[case("int1e_kin"      , "s2ij", [[8, 18], [8, 18]]  , 2.0474763079697453    , [496]         )] // no-deriv | shl | s2ij
    #[case("int1e_ipnuc"    , "s2ij", [[8, 18], [8, 18]]  , 43.07817054892782     , [496, 3]      )] // deriv | shl | s2ij
    fn test_2c_cart(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl Into<ShlsSlice>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_h2o_def2_tzvp();

        cint_data.with_cint_type(Cartesian, |cint_data| {
            let (out, shape) = cint_data.integrate(intor, aosym, shls_slice).into();
            assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-12);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    #[rstest]
    #[case("int2c2e"        , "s1"  , None                , -82.66855640793396    , 0.0                   , [86, 86]      )] // no-deriv
    #[case("int1e_igovlp"   , "s1"  , None                , -1.989836876952396    , -2.5513536327277677   , [86, 86, 3]   )] // deriv
    #[case("int1e_kin"      , "s1"  , [[0, 12], [8, 18]]  , -14.325037173756517   , 0.0                   , [64, 52]      )] // no-deriv | shl
    #[case("int1e_ipnuc"    , "s1"  , [[0, 12], [8, 18]]  , -66.46315988601708    , 7.7000121009088165    , [64, 52, 3]   )] // deriv | shl
    #[case("int2c2e"        , "s2ij", None                , -108.42520152395204   , 0.0                   , [3741]        )] // no-deriv | s2ij
    #[case("int1e_igovlp"   , "s2ij", None                , -0.43544614746500226  , -4.7372500914607505   , [3741, 3]     )] // deriv | s2ij
    #[case("int1e_kin"      , "s2ij", [[8, 18], [8, 18]]  , -6.866472862080281    , 0.0                   , [1378]         )] // no-deriv | shl | s2ij
    #[case("int1e_ipnuc"    , "s2ij", [[8, 18], [8, 18]]  , -16.294988603262368   , 21.630652439596524    , [1378, 3]      )] // deriv | shl | s2ij
    fn test_2c_spinor(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl Into<ShlsSlice>,
        #[case] ref_fp_re: f64,
        #[case] ref_fp_im: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_h2o_def2_tzvp();

        cint_data.with_cint_type(Spinor, |cint_data| {
            let (out, shape) = cint_data.integrate_spinor(intor, aosym, shls_slice).into();
            let fp = cint_fp(&out);
            assert_relative_eq!(ref_fp_re, fp.re(), epsilon = 1e-12);
            assert_relative_eq!(ref_fp_im, fp.im(), epsilon = 1e-12);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    /* #endregion */

    /* #region 3-center */

    #[rstest]
    #[case("int3c2e"        , "s1"  , None                         , 385.7604793373525     , [43, 43, 43]      )] // no-deriv
    #[case("int3c2e_ig1"    , "s1"  , None                         , 17.403936599238815    , [43, 43, 43, 3]   )] // deriv
    #[case("int3c2e"        , "s1"  , [[0, 12], [8, 18], [6, 15]]  , 3.365987495786978     , [32, 26, 29]      )] // no-deriv | shl
    #[case("int3c2e_ipip1"  , "s1"  , [[0, 12], [8, 18], [6, 15]]  , 29.05796504400456     , [32, 26, 29, 9]   )] // deriv | shl
    #[case("int3c2e"        , "s2ij", None                         , 8.973911146875718     , [946, 43]         )] // no-deriv | s2ij
    #[case("int3c2e_ip2"    , "s2ij", None                         , -4.835348657362738    , [946, 43, 3]      )] // deriv | s2ij
    #[case("int3c2e"        , "s2ij", [[8, 18], [8, 18], [6, 15]]  , 30.29313979284533     , [351, 29]         )] // no-deriv | shl | s2ij
    #[case("int3c2e_ip1ip2" , "s2ij", [[8, 18], [8, 18], [6, 15]]  , -6.144703885051337    , [351, 29, 9]      )] // deriv | shl | s2ij
    fn test_3c_sph(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl Into<ShlsSlice>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_h2o_def2_tzvp();

        cint_data.with_cint_type(Spheric, |cint_data| {
            let (out, shape) = cint_data.integrate(intor, aosym, shls_slice).into();
            assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-12);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    #[rstest]
    #[case("int3c2e"        , "s1"  , None                         , -139.11319865592063   , [48, 48, 48]      )] // no-deriv
    #[case("int3c2e_ig1"    , "s1"  , None                         , 5.9023441440830915    , [48, 48, 48, 3]   )] // deriv
    #[case("int3c2e"        , "s1"  , [[0, 12], [8, 18], [6, 15]]  , 53.55901107470838     , [37, 31, 34]      )] // no-deriv | shl
    #[case("int3c2e_ipip1"  , "s1"  , [[0, 12], [8, 18], [6, 15]]  , 904.5370718832919     , [37, 31, 34, 9]   )] // deriv | shl
    #[case("int3c2e"        , "s2ij", None                         , 80.25669268335058     , [1176, 48]        )] // no-deriv | s2ij
    #[case("int3c2e_ip2"    , "s2ij", None                         , 1.3547693151160463    , [1176, 48, 3]     )] // deriv | s2ij
    #[case("int3c2e"        , "s2ij", [[8, 18], [8, 18], [6, 15]]  , -236.8269201943644    , [496, 34]         )] // no-deriv | shl | s2ij
    #[case("int3c2e_ip1ip2" , "s2ij", [[8, 18], [8, 18], [6, 15]]  , -34.75948218463372    , [496, 34, 9]      )] // deriv | shl | s2ij
    fn test_3c_cart(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl Into<ShlsSlice>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_h2o_def2_tzvp();

        cint_data.with_cint_type(Cartesian, |cint_data| {
            let (out, shape) = cint_data.integrate(intor, aosym, shls_slice).into();
            assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-12);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    /* #endregion */

    /* #region 4-center */

    #[rstest]
    #[case("int2e"              , "s1"  , None                                  , 70.00106603114827     , [43, 43, 43, 43]      )] // no-deriv
    #[case("int2e_ip1"          , "s1"  , None                                  , 118.50663866543611    , [43, 43, 43, 43, 3]   )] // deriv
    #[case("int2e"              , "s1"  , [[0, 12], [8, 18], [6, 15], [3, 10]]  , 1.5094112713177106    , [32, 26, 29, 21]      )] // no-deriv | shl
    #[case("int2e_cg_ssa10ssp2" , "s1"  , [[0, 12], [8, 18], [6, 15], [3, 10]]  , 36.85088072685591     , [32, 26, 29, 21, 48]  )] // deriv | shl
    #[case("int2e"              , "s2ij", None                                  , 1.8306651475961022    , [946, 43, 43]         )] // no-deriv | s2ij
    #[case("int2e_gg1"          , "s2ij", None                                  , -12.331846851227873   , [946, 43, 43, 9]      )] // deriv | s2ij
    #[case("int2e"              , "s2ij", [[8, 18], [8, 18], [6, 15], [3, 10]]  , 1.7525131209661868    , [351, 29, 21]         )] // no-deriv | shl | s2ij
    #[case("int2e_gg1"          , "s2ij", [[8, 18], [8, 18], [6, 15], [3, 10]]  , 6.391595640883066     , [351, 29, 21, 9]      )] // deriv | shl | s2ij
    fn test_4c_sph(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl Into<ShlsSlice>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_h2o_def2_tzvp();

        cint_data.with_cint_type(Spheric, |cint_data| {
            let (out, shape) = cint_data.integrate(intor, aosym, shls_slice).into();
            assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-11);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    #[rstest]
    #[case("int2e"              , "s1"  , None                                  , -27.40361459334603    , [48, 48, 48, 48]      )] // no-deriv
    #[case("int2e_ip1"          , "s1"  , None                                  , -102.69339786816539   , [48, 48, 48, 48, 3]   )] // deriv
    #[case("int2e"              , "s1"  , [[0, 12], [8, 18], [6, 15], [3, 10]]  , -3.3950467582181965   , [37, 31, 34, 23]      )] // no-deriv | shl
    #[case("int2e_cg_ssa10ssp2" , "s1"  , [[0, 12], [8, 18], [6, 15], [3, 10]]  , -23.165032397305538   , [37, 31, 34, 23, 48]  )] // deriv | shl
    #[case("int2e"              , "s2ij", None                                  , 40.18960951457926     , [1176, 48, 48]        )] // no-deriv | s2ij
    #[case("int2e_gg1"          , "s2ij", None                                  , 53.03167458535967     , [1176, 48, 48, 9]     )] // deriv | s2ij
    #[case("int2e"              , "s2ij", [[8, 18], [8, 18], [6, 15], [3, 10]]  , -84.37906159867836    , [496, 34, 23]         )] // no-deriv | shl | s2ij
    #[case("int2e_gg1"          , "s2ij", [[8, 18], [8, 18], [6, 15], [3, 10]]  , -18.774888488360734   , [496, 34, 23, 9]      )] // deriv | shl | s2ij
    fn test_4c_cart(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl Into<ShlsSlice>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_h2o_def2_tzvp();

        cint_data.with_cint_type(Cartesian, |cint_data| {
            let (out, shape) = cint_data.integrate(intor, aosym, shls_slice).into();
            assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-11);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    #[rstest]
    #[case("int2e"              , "s1"  , None                                  , -4.074835677031094    , 0.0                   , [86, 86, 86, 86]      )] // no-deriv
    #[case("int2e_ip1"          , "s1"  , None                                  , 65.47856903615192     , 3.951167514490777     , [86, 86, 86, 86, 3]   )] // deriv
    #[case("int2e"              , "s1"  , [[0, 12], [8, 18], [6, 15], [3, 10]]  , 6.7689010236720835    , 0.0                   , [64, 52, 58, 42]      )] // no-deriv | shl
    #[case("int2e_cg_ssa10ssp2" , "s1"  , [[0, 12], [8, 18], [6, 15], [3, 10]]  , 169.1237983415611     , -52.91956958665179    , [64, 52, 58, 42, 3]   )] // deriv | shl
    fn test_4c_spinor(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl Into<ShlsSlice>,
        #[case] ref_fp_re: f64,
        #[case] ref_fp_im: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_h2o_def2_tzvp();

        cint_data.with_cint_type(Spinor, |cint_data| {
            let (out, shape) = cint_data.integrate_spinor(intor, aosym, shls_slice).into();
            let fp = cint_fp(&out);
            assert_relative_eq!(ref_fp_re, fp.re(), epsilon = 1e-11);
            assert_relative_eq!(ref_fp_im, fp.im(), epsilon = 1e-11);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    #[rstest]
    // s2kl
    #[case("int2e"              , "s2kl", None                                  , -135.42091888440837   , [43, 43, 946]         )] // no-deriv
    #[case("int2e_ip1"          , "s2kl", None                                  , 71.84356356407886     , [43, 43, 946, 3]      )] // deriv
    #[case("int2e"              , "s2kl", [[0, 12], [8, 18], [6, 15], [6, 15]]  , -11.318569065551744   , [32, 26, 435]         )] // no-deriv | shl
    #[case("int2e_ip1"          , "s2kl", [[0, 12], [8, 18], [6, 15], [6, 15]]  , 30.740874132068775    , [32, 26, 435, 3]      )] // deriv
    fn test_4c_other_sym_sph(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl Into<ShlsSlice>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_h2o_def2_tzvp();

        cint_data.with_cint_type(Spheric, |cint_data| {
            let (out, shape) = cint_data.integrate(intor, aosym, shls_slice).into();
            assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-11);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    /* #endregion */
}
