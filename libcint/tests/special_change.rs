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
    fn test_concat_frac_charge() {
        use libcint::ffi::cint_ffi;
        use std::ffi::c_int;
        const NUC_MOD_OF: usize = cint_ffi::NUC_MOD_OF as usize;
        const PTR_FRAC_CHARGE: usize = cint_ffi::PTR_FRAC_CHARGE as usize;

        let mol1 = init_h2o_def2_tzvp();
        let mut mol2 = init_h2o_def2_tzvp();
        mol2.atm[0][NUC_MOD_OF] = cint_ffi::FRAC_CHARGE_NUC as c_int;
        mol2.atm[0][PTR_FRAC_CHARGE] = mol2.env.len() as c_int;
        mol2.env.push(-0.834);
        assert_eq!(mol2.atom_charge(0), -0.834);

        let mol = &mol1 + &mol2;
        assert_eq!(mol.atom_charges()[..3], mol1.atom_charges()[..]);
        assert_eq!(mol.atom_charges()[3..], mol2.atom_charges()[..]);
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

    #[test]
    fn test_stale_ecpbas_offset() {
        use libcint::ffi::cecp_ffi;
        const AS_ECPBAS_OFFSET: usize = cecp_ffi::AS_ECPBAS_OFFSET as usize;
        const AS_NECPBAS: usize = cecp_ffi::AS_NECPBAS as usize;

        let mol = init_sb2me4_cc_pvtz();
        let merged = mol.merge_ecpbas();
        assert!(!mol.is_ecp_merged());
        assert!(merged.is_ecp_merged());

        // PySCF's `Mole._env` after an ECP integral: offset slots set, `_bas`
        // unmerged
        let mut stale = mol.clone();
        stale.env[AS_ECPBAS_OFFSET] = mol.bas.len() as f64;
        stale.env[AS_NECPBAS] = mol.ecpbas.len() as f64;
        assert!(!stale.is_ecp_merged());
        assert_eq!(stale.nbas(), mol.nbas());
        assert_eq!(stale.decopule_ecpbas().ecpbas, mol.ecpbas);
        assert_eq!(stale.merge_ecpbas(), merged);
        let (out, shape) = stale.integrate("int1e_ovlp", None, None).into();
        let (out_ref, shape_ref) = mol.integrate("int1e_ovlp", None, None).into();
        assert_eq!(shape, shape_ref);
        assert_eq!(out, out_ref);
    }
}
