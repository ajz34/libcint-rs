use approx::assert_relative_eq;
use libcint::prelude::*;
use rstest::rstest;

use CIntType::{Cartesian, Spheric};

#[cfg(test)]
mod test {
    use super::*;

    /* #region test usual case */

    #[rstest]
    #[case("int2c2e"        , "s1"  , []                  , -935.2694687735325    , [366, 366]    )] // no-deriv
    #[case("int1e_igovlp"   , "s1"  , []                  , -106.16961857524265   , [366, 366, 3] )] // deriv
    #[case("int1e_kin"      , "s1"  , [[0, 12], [8, 18]]  , -0.35004852378214735  , [41, 34]      )] // no-deriv | shl
    #[case("int1e_ipnuc"    , "s1"  , [[0, 12], [8, 18]]  , 45.968244650324706    , [41, 34, 3]   )] // deriv | shl
    #[case("int2c2e"        , "s2ij", []                  , 1809.6733664015762    , [67161]       )] // no-deriv | s2ij
    #[case("int1e_igovlp"   , "s2ij", []                  , -40.08094621810594    , [67161, 3]    )] // deriv | s2ij
    #[case("int1e_kin"      , "s2ij", [[8, 18], [8, 18]]  , 0.5049496096560406    , [595]         )] // no-deriv | shl | s2ij
    #[case("int1e_ipnuc"    , "s2ij", [[8, 18], [8, 18]]  , 4.499859812123924     , [595, 3]      )] // deriv | shl | s2ij
    fn test_2c_sph(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl AsRef<[[usize; 2]]>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_sb2me4_cc_pvtz();

        cint_data.with_cint_type(Spheric, |cint_data| {
            let (out, shape) = cint_data.integrate(intor, aosym, shls_slice).into();
            assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-10);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    #[rstest]
    #[case("int2c2e"        , "s1"  , []                  , -935.2694687735325    , [366, 366]    )] // no-deriv
    #[case("int1e_igovlp"   , "s1"  , []                  , -106.16961857524265   , [366, 366, 3] )] // deriv
    #[case("int1e_kin"      , "s1"  , [[0, 12], [8, 18]]  , -0.35004852378214735  , [41, 34]      )] // no-deriv | shl
    #[case("int1e_ipnuc"    , "s1"  , [[0, 12], [8, 18]]  , 45.968244650324706    , [41, 34, 3]   )] // deriv | shl
    #[case("int2c2e"        , "s2ij", []                  , 1809.6733664015762    , [67161]       )] // no-deriv | s2ij
    #[case("int1e_igovlp"   , "s2ij", []                  , -40.08094621810594    , [67161, 3]    )] // deriv | s2ij
    #[case("int1e_kin"      , "s2ij", [[8, 18], [8, 18]]  , 0.5049496096560406    , [595]         )] // no-deriv | shl | s2ij
    #[case("int1e_ipnuc"    , "s2ij", [[8, 18], [8, 18]]  , 4.499859812123924     , [595, 3]      )] // deriv | shl | s2ij
    fn test_2c_sph_ecp_merged(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl AsRef<[[usize; 2]]>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_sb2me4_cc_pvtz().merge_ecpbas();

        cint_data.with_cint_type(Spheric, |cint_data| {
            let (out, shape) = cint_data.integrate(intor, aosym, shls_slice).into();
            assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-10);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    /* #endregion */

    /* #region ECP */

    #[rstest]
    #[case("ECPscalar"         , "s1"  , []                  , 0.7534281412150909    , [366, 366]    )] // no-deriv
    #[case("ECPscalar_ignuc"   , "s1"  , []                  , -1.3426139111770714   , [366, 366, 3] )] // deriv
    #[case("ECPscalar"         , "s1"  , [[0, 12], [8, 18]]  , 1.4186154309495258    , [41, 34]      )] // no-deriv | shl
    #[case("ECPscalar_iprinv"  , "s1"  , [[0, 12], [8, 18]]  , -11.738302192911995   , [41, 34, 3]   )] // deriv | shl
    #[case("ECPscalar"         , "s2ij", []                  , -1.8388155709304113   , [67161]       )] // no-deriv | s2ij
    #[case("ECPscalar_iprinvip", "s2ij", []                  , -68.89774165957373    , [67161, 9]    )] // deriv | s2ij
    #[case("ECPscalar"         , "s2ij", [[8, 18], [8, 18]]  , 15.651229742132026    , [595]         )] // no-deriv | shl | s2ij
    #[case("ECPscalar_iprinvip", "s2ij", [[8, 18], [8, 18]]  , 0.5549200392488471    , [595, 9]      )] // deriv | shl | s2ij
    fn test_ecp_sph(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl AsRef<[[usize; 2]]>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_sb2me4_cc_pvtz();

        cint_data.with_cint_type(Spheric, |cint_data| {
            let (out, shape) = cint_data.integrate(intor, aosym, shls_slice).into();
            assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-10);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    #[rstest]
    #[case("ECPscalar"         , "s1"  , []                  , 1.6418489688115985    , [410, 410]    )] // no-deriv
    #[case("ECPscalar_ignuc"   , "s1"  , []                  , -16.389862409392386   , [410, 410, 3] )] // deriv
    #[case("ECPscalar"         , "s1"  , [[0, 12], [8, 18]]  , 1.3845811428580408    , [47, 39]      )] // no-deriv | shl
    #[case("ECPscalar_iprinv"  , "s1"  , [[0, 12], [8, 18]]  , -7.983349456393638    , [47, 39, 3]   )] // deriv | shl
    #[case("ECPscalar"         , "s2ij", []                  , -9.966466894809031    , [84255]       )] // no-deriv | s2ij
    #[case("ECPscalar_iprinvip", "s2ij", []                  , 132.938348145421      , [84255, 9]    )] // deriv | s2ij
    #[case("ECPscalar"         , "s2ij", [[8, 18], [8, 18]]  , -20.31359441160549    , [780]         )] // no-deriv | shl | s2ij
    #[case("ECPscalar_iprinvip", "s2ij", [[8, 18], [8, 18]]  , 1.9270442333836884    , [780, 9]      )] // deriv | shl | s2ij
    fn test_ecp_cart(
        #[case] intor: &str,
        #[case] aosym: &str,
        #[case] shls_slice: impl AsRef<[[usize; 2]]>,
        #[case] ref_fp: f64,
        #[case] ref_shape: impl AsRef<[usize]>,
    ) {
        let mut cint_data = init_sb2me4_cc_pvtz();

        cint_data.with_cint_type(Cartesian, |cint_data| {
            let (out, shape) = cint_data.integrate(intor, aosym, shls_slice).into();
            assert_relative_eq!(ref_fp, cint_fp(&out), epsilon = 1e-10);
            assert_eq!(shape, ref_shape.as_ref());
        });
    }

    /* #endregion */
}
