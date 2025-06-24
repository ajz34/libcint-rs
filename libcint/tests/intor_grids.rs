use approx::assert_relative_eq;
use itertools::Itertools;
use libcint::prelude::*;
use rstest::rstest;

#[cfg(test)]
mod test {
    use super::*;

    #[rstest]
    #[case("int1e_grids"   , -0.9862020295475049, [420, 43, 43]   )]
    #[case("int1e_grids_ip", 1.344336345578797  , [420, 43, 43, 3])]
    fn test_col_major(#[case] intor: &str, #[case] ref_fp: f64, #[case] ref_shape: impl AsRef<[usize]>) {
        let mut cint_data = init_h2o_def2_tzvp();

        let grids = (0..420).map(|x| [-2.3 + 0.012 * x as f64, 1.5 - 0.007 * x as f64, 0.2 + 0.015 * x as f64]).collect_vec();
        let (out, shape) = cint_data.with_grids(&grids, |data| data.integrate(intor, None, None).into());
        assert_relative_eq!(cint_fp(&out), ref_fp, epsilon = 1e-12);
        assert_eq!(shape, ref_shape.as_ref());
    }

    #[rstest]
    #[case("int1e_grids"   , 3.916348193223265 , [420, 43, 43]   )]
    #[case("int1e_grids_ip", -7.757919630556912, [3, 420, 43, 43])]
    fn test_row_major(#[case] intor: &str, #[case] ref_fp: f64, #[case] ref_shape: impl AsRef<[usize]>) {
        let mut cint_data = init_h2o_def2_tzvp();

        let grids = (0..420).map(|x| [-2.3 + 0.012 * x as f64, 1.5 - 0.007 * x as f64, 0.2 + 0.015 * x as f64]).collect_vec();
        let (out, shape) = cint_data.with_grids(&grids, |data| data.integrate_row_major(intor, None, None).into());
        assert_relative_eq!(cint_fp(&out), ref_fp, epsilon = 1e-12);
        assert_eq!(shape, ref_shape.as_ref());
    }

    #[test]
    #[should_panic]
    fn test_panic() {
        let cint_data = init_h2o_def2_tzvp();
        let (_out, _shape) = cint_data.integrate("int1e_grids", None, None).into();
    }
}
