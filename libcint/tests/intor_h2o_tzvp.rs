use libcint::prelude::*;

#[cfg(test)]
mod test {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn test_2c_s1() {
        let cint_data = init_h2o_def2_tzvp();

        // usual case
        let (out, shape) = cint_data.integrate("int2c2e", "s1", []).into();
        assert_relative_eq!(353.92327889575165, cint_fp(&out), epsilon = 1e-12);
        assert_eq!(shape, [43, 43]);
    }
}
