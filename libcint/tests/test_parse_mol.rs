//! Tests for CIntMol molecule object generation.
//!
//! Ported from PySCF's pyscf/gto/test/test_mole.py
//! - test_ghost: ghost atom handling
//! - test_zmat: z-matrix conversion
//! - test_atom_types: labeled atoms with per-atom basis
//!
//! Note: Most tests are in libcint/src/parse/atom.rs, basis.rs, mol.rs
//! This file only contains integration tests that need the full molecule build.

#[cfg(feature = "bse")]
mod test_mol_integration {
    use approx::assert_relative_eq;
    use libcint::parse::atom::{parse_atom_string, parse_zmatrix, Unit, ANG_TO_BOHR};
    use libcint::parse::mol::{CIntMolInput, M_f, M};

    /// Test z-matrix parsing.
    /// Ported from PySCF test_mole.py::test_zmat
    #[test]
    fn test_zmatrix_parsing() {
        // Simple water molecule in z-matrix format
        let atoms = parse_zmatrix(
            "O
             H 1 0.94
             H 1 0.94 2 104.5",
            Unit::Angstrom,
        )
        .unwrap();

        assert_eq!(atoms.len(), 3);
        assert_eq!(atoms[0].symbol, "O");
        assert_eq!(atoms[1].symbol, "H");
        assert_eq!(atoms[2].symbol, "H");

        // First atom at origin
        assert_relative_eq!(atoms[0].coords[0], 0.0, epsilon = 1e-10);
        assert_relative_eq!(atoms[0].coords[1], 0.0, epsilon = 1e-10);
        assert_relative_eq!(atoms[0].coords[2], 0.0, epsilon = 1e-10);

        // Second atom along some axis (bond length)
        let bond_bohr = 0.94 * ANG_TO_BOHR;
        let dist1 = (atoms[1].coords[0].powi(2) + atoms[1].coords[1].powi(2) + atoms[1].coords[2].powi(2)).sqrt();
        assert_relative_eq!(dist1, bond_bohr, epsilon = 0.01);

        // Third atom with angle
        let dist2 = (atoms[2].coords[0].powi(2) + atoms[2].coords[1].powi(2) + atoms[2].coords[2].powi(2)).sqrt();
        assert_relative_eq!(dist2, bond_bohr, epsilon = 0.01);
    }

    /// Test labeled atoms with per-atom basis.
    /// Ported from PySCF test_mole.py::test_atom_types
    #[test]
    fn test_labeled_atoms() {
        // Parse labeled atoms: H0, H1, H, H3
        let atoms = parse_atom_string(
            "H0  0 0 0
             H1  0 0 0
             H   0 0 0
             H3  0 0 0",
            Unit::Angstrom,
        )
        .unwrap();

        assert_eq!(atoms.len(), 4);
        // All should be H atoms
        assert_eq!(atoms[0].symbol, "H");
        assert_eq!(atoms[1].symbol, "H");
        assert_eq!(atoms[2].symbol, "H");
        assert_eq!(atoms[3].symbol, "H");
        // Labels should be preserved
        assert_eq!(atoms[0].label, "H0");
        assert_eq!(atoms[1].label, "H1");
        assert_eq!(atoms[2].label, "H");
        assert_eq!(atoms[3].label, "H3");
    }

    /// Test M() convenience function.
    #[test]
    fn test_m_convenience() {
        let mol = M("H 0 0 0; H 0 0 0.74", "sto-3g");
        assert_eq!(mol.atoms.len(), 2);
        assert_eq!(mol.cint.bas.len(), 2);
    }

    /// Test M_f() fallible version.
    #[test]
    fn test_m_f_convenience() {
        let mol = M_f("H 0 0 0; H 0 0 0.74", "sto-3g").unwrap();
        assert_eq!(mol.atoms.len(), 2);
        assert_eq!(mol.cint.bas.len(), 2);
    }

    /// Test build_f() fallible version.
    #[test]
    fn test_build_f() {
        let mol = CIntMolInput::new().atom_str("O 0 0 0; H 0 0.94 0; H 0.94 0 0").basis_str("sto-3g").angstrom().build_f().unwrap();

        assert_eq!(mol.atoms.len(), 3);
        assert_eq!(mol.cint.atm.len(), 3);
    }
}
