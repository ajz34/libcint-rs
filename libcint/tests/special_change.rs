// Reference source: libcint/scripts/pyscf_special_change.ipynb

use approx::assert_relative_eq;
use libcint::prelude::*;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_concat_molecules() {
        let mol1 = init_h2o_def2_tzvp();
        let mol2 = init_sb2me4_cc_pvtz();

        let mol = &mol1 + &mol1;
        assert_relative_eq!(cint_fp(&mol.atm), 220.2819072746202, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&mol.bas), 149.65011901575278, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&mol.env), -7210.728997361069, max_relative = 1e-10);

        let mol = &mol1 + &mol2;
        assert_relative_eq!(cint_fp(&mol.atm), -68.54444636475233, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&mol.bas), -58.06413507172557, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&mol.env), 26407.797117689446, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&mol.ecpbas), 164.89860996992238, max_relative = 1e-10);

        let mol = &mol2 + &mol1;
        assert_relative_eq!(cint_fp(&mol.atm), -714.6455874682379, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&mol.bas), -99.88753072795852, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&mol.env), -32210.1073858801, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&mol.ecpbas), 110.04541227168156, max_relative = 1e-10);

        let mol = &mol2 + &mol2;
        assert_relative_eq!(cint_fp(&mol.atm), 1355.7032887239495, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&mol.bas), -589.9932899367334, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&mol.env), 1519.3261601043766, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&mol.ecpbas), 779.9379788456364, max_relative = 1e-10);
    }

    #[test]
    fn test_fakemol_for_charges() {
        let coords = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [2.0, 1.5, 2.3], [3.0, 4.0, 5.0]];
        let exponent = 0.5;
        let fakemol = fakemol_for_charges(&coords, exponent);

        assert_relative_eq!(cint_fp(&fakemol.atm), 80.41184821237731, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&fakemol.bas), 3.9145826029379323, max_relative = 1e-10);
        assert_relative_eq!(cint_fp(&fakemol.env), 1.4694722804235696, max_relative = 1e-10);

        let mol = init_h2o_def2_tzvp();
        let (out, shape) = CInt::integrate_cross("int1e_ovlp", [&mol, &fakemol], "s1", None).into();
        assert_relative_eq!(cint_fingerprint(&out), 0.3218956097258599, max_relative = 1e-10);
        assert_eq!(shape, [43, 4]);
    }
}
