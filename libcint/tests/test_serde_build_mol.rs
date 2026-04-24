#![cfg(feature = "bse")]

use libcint::parse::mole::CIntMol;

fn cint_fingerprint(data: &[f64]) -> f64 {
    let n = data.len();
    if n == 0 {
        return 0.0;
    }
    let mut sum = 0.0;
    for (i, &v) in data.iter().enumerate() {
        sum += (i as f64 + 1.0) * v.abs();
    }
    sum / n as f64
}

// Test basic JSON parsing
#[test]
fn test_json_basic() {
    let json = r#"{
        "atom": "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987",
        "basis": "STO-3G",
        "unit": "angstrom"
    }"#;

    let mol = CIntMol::from_json(json);

    let nao = mol.cint.nao();
    assert_eq!(nao, 7);

    let (int1e_ovlp, shape) = mol.cint.integrate("int1e_ovlp", None, None).into();
    assert_eq!(shape[0], 7);
    assert!(int1e_ovlp.iter().any(|&v| v.abs() > 0.0));
}

// Test basic TOML parsing
#[test]
fn test_toml_basic() {
    let toml = r#"
atom = "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987"
basis = "STO-3G"
unit = "angstrom"
"#;

    let mol = CIntMol::from_toml(toml);

    let nao = mol.cint.nao();
    assert_eq!(nao, 7);
}

// Test TOML with multi-line atom string
#[test]
fn test_toml_multiline_atom() {
    let toml = r#"
atom = """
O 0.000000 0.000000 0.000000
H 0.000000 0.000000 0.957200
H 0.000000 0.926627 -0.239987
"""
basis = "STO-3G"
"#;

    let mol = CIntMol::from_toml(toml);
    let nao = mol.cint.nao();
    assert_eq!(nao, 7);
}

// Test TOML with all optional fields
#[test]
fn test_toml_all_fields() {
    let toml = r#"
atom = "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987"
basis = "STO-3G"
ecp = ""
unit = "angstrom"
cart = false
ghost_ecp = false
allow_empty_basis = false
"#;

    let mol = CIntMol::from_toml(toml);
    let nao = mol.cint.nao();
    assert_eq!(nao, 7);
}

// Test TOML with Bohr unit
#[test]
fn test_toml_bohr_unit() {
    let toml = r#"
atom = "O 0 0 0; H 0 0 1.8; H 0 1.7 -0.4"
basis = "STO-3G"
unit = "bohr"
"#;

    let mol = CIntMol::from_toml(toml);
    let nao = mol.cint.nao();
    assert_eq!(nao, 7);
}

// Test TOML with custom basis
#[test]
fn test_toml_custom_basis() {
    let toml = r#"
atom = "O 0 0 0; H 0 0 0.9572"
basis = "custom"

[basis-custom]
O = "STO-3G"
default = "6-31G"
"#;

    let mol = CIntMol::from_toml(toml);

    let nao = mol.cint.nao();
    assert!(nao > 0);

    let (int1e_ovlp, _shape) = mol.cint.integrate("int1e_ovlp", None, None).into();
    assert!(int1e_ovlp.iter().any(|&v| v.abs() > 0.0));
}

// Test TOML with custom ECP
#[test]
fn test_toml_custom_ecp() {
    let toml = r#"
atom = "Sb 0 0 0; H 0 0 1.5"
basis = "def2-SVP"
ecp = "custom"

[ecp-custom]
Sb = "def2-SVP"
"#;

    let mol = CIntMol::from_toml(toml);

    let nao = mol.cint.nao();
    assert!(nao > 0);
}

// Test TOML with both custom basis and ECP
#[test]
fn test_toml_custom_basis_and_ecp() {
    let toml = r#"
atom = "Sb 0 0 0; H 0 0 1.5"
basis = "custom"
ecp = "custom"

[basis-custom]
Sb = "def2-SVP"
H = "def2-SVP"
default = "STO-3G"

[ecp-custom]
Sb = "def2-SVP"
"#;

    let mol = CIntMol::from_toml(toml);
    let nao = mol.cint.nao();
    assert!(nao > 0);
}

// Test JSON with inline dict basis
#[test]
fn test_json_dict_basis() {
    let json = r#"{
        "atom": "O 0 0 0; H 0 0 0.9572",
        "basis": {
            "O": "STO-3G",
            "default": "6-31G"
        },
        "unit": "angstrom"
    }"#;

    let mol = CIntMol::from_json(json);
    let nao = mol.cint.nao();
    assert!(nao > 0);
}

// Test JSON with inline dict ECP
#[test]
fn test_json_dict_ecp() {
    let json = r#"{
        "atom": "Sb 0 0 0; H 0 0 1.5",
        "basis": "def2-SVP",
        "ecp": {
            "Sb": "def2-SVP"
        },
        "unit": "angstrom"
    }"#;

    let mol = CIntMol::from_json(json);
    let nao = mol.cint.nao();
    assert!(nao > 0);
}

// Test that JSON with "custom" basis fails
#[test]
fn test_json_custom_basis_error() {
    let json = r#"{
        "atom": "O 0 0 0; H 0 0 0.9572",
        "basis": "custom"
    }"#;

    let result = CIntMol::from_json_f(json);
    assert!(result.is_err());
}

// Test that TOML with custom basis but no custom table fails
#[test]
fn test_toml_custom_basis_missing_table_error() {
    let toml = r#"
atom = "O 0 0 0; H 0 0 0.9572"
basis = "custom"
"#;

    let result = CIntMol::from_toml_f(toml);
    assert!(result.is_err());
}

// Test cartesian basis
#[test]
fn test_toml_cartesian() {
    let toml = r#"
atom = "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987"
basis = "6-31G"
unit = "angstrom"
cart = true
"#;

    let mol = CIntMol::from_toml(toml);
    let nao = mol.cint.nao();
    assert_eq!(nao, 13);
}

// Test z-matrix format
#[test]
fn test_toml_zmatrix() {
    let toml = r#"
atom = "O; H 1 0.94; H 1 0.94 2 104.5"
basis = "def2-TZVP"
unit = "angstrom"
"#;

    let mol = CIntMol::from_toml(toml);
    let nao = mol.cint.nao();
    assert_eq!(nao, 43);
}

// Test ghost atoms
#[test]
fn test_toml_ghost_atoms() {
    let toml = r#"
atom = "O 0 0 0; GHOST-H 0 0 1.0; H 0 0 1.0"
basis = "STO-3G"
unit = "angstrom"
allow_empty_basis = true
"#;

    let mol = CIntMol::from_toml(toml);
    let nao = mol.cint.nao();
    assert!(nao > 0);
}

// Test labeled atoms with custom basis
#[test]
fn test_toml_labeled_atoms_custom_basis() {
    let toml = r#"
atom = "O 0 0 0; H1 0 0 0.9572; H2 0 0.9266 -0.239987"
basis = "custom"

[basis-custom]
O = "STO-3G"
H1 = "6-31G"
H2 = "def2-SVP"
"#;

    let mol = CIntMol::from_toml(toml);
    let nao = mol.cint.nao();
    assert!(nao > 0);
}

// Test large molecule with ECP
#[test]
fn test_toml_large_ecp() {
    let toml = r#"
atom = """
Sb -1.33937843 0.44597852 -1.27279684
Sb1 1.33937843 -0.44597852 -1.27279684
C -1.40429524 1.10441871 0.83468205
"""
basis = "def2-SVP"
ecp = "def2-SVP"
unit = "angstrom"
allow_empty_basis = true
"#;

    let mol = CIntMol::from_toml(toml);
    let nao = mol.cint.nao();
    assert!(nao > 0);
}

// Test fingerprint consistency between JSON and TOML
#[test]
fn test_json_toml_fingerprint_match() {
    let json = r#"{
        "atom": "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987",
        "basis": "STO-3G",
        "unit": "angstrom"
    }"#;

    let toml = r#"
atom = "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987"
basis = "STO-3G"
unit = "angstrom"
"#;

    let mol_json = CIntMol::from_json(json);
    let mol_toml = CIntMol::from_toml(toml);

    let ovlp_json = mol_json.cint.integrate("int1e_ovlp", None, None).out.unwrap();
    let ovlp_toml = mol_toml.cint.integrate("int1e_ovlp", None, None).out.unwrap();

    let fp_json = cint_fingerprint(ovlp_json.as_slice());
    let fp_toml = cint_fingerprint(ovlp_toml.as_slice());

    assert!((fp_json - fp_toml).abs() < 1e-10);
}

// Test unit case insensitivity
#[test]
fn test_toml_unit_case_insensitive() {
    let toml_upper = r#"
atom = "O 0 0 0; H 0 0 1.8; H 0 1.7 -0.4"
basis = "STO-3G"
unit = "BOHR"
"#;

    let toml_lower = r#"
atom = "O 0 0 0; H 0 0 1.8; H 0 1.7 -0.4"
basis = "STO-3G"
unit = "bohr"
"#;

    let mol_upper = CIntMol::from_toml(toml_upper);
    let mol_lower = CIntMol::from_toml(toml_lower);

    assert_eq!(mol_upper.cint.nao(), mol_lower.cint.nao());
}