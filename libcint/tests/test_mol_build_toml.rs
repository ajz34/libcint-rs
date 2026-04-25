#![cfg(feature = "bse")]

use libcint::prelude::*;

fn check_cint_mol(toml_input: &str, refs: (usize, usize, f64, f64)) {
    let mol = CIntMol::from_toml(toml_input);
    let (ref_nao, ref_nbas, ref_kin_fp, ref_nuc_fp) = refs;

    let nao = mol.cint.nao();
    let nbas = mol.cint.nbas();
    let kin_fp = cint_fingerprint(mol.cint.integrate("int1e_kin", None, None).out.as_ref().unwrap());
    let nuc_fp = cint_fingerprint(mol.cint.integrate("int1e_nuc", None, None).out.as_ref().unwrap());

    assert_eq!(nao, ref_nao);
    assert_eq!(nbas, ref_nbas);
    assert!(
        (kin_fp - ref_kin_fp).abs() < 1e-6 || (kin_fp / ref_kin_fp - 1.0).abs() < 1e-4,
        "kinetic fingerprint mismatch: got {kin_fp}, expected {ref_kin_fp}"
    );
    assert!(
        (nuc_fp - ref_nuc_fp).abs() < 1e-6 || (nuc_fp / ref_nuc_fp - 1.0).abs() < 1e-4,
        "nuclear fingerprint mismatch: got {nuc_fp}, expected {ref_nuc_fp}"
    );
}

fn check_cint_mol_ecp(toml_input: &str, refs: (usize, usize, f64, f64, f64)) {
    let mol = CIntMol::from_toml(toml_input);
    let (ref_nao, ref_nbas, ref_kin_fp, ref_nuc_fp, ref_ecp_fp) = refs;

    let nao = mol.cint.nao();
    let nbas = mol.cint.nbas();
    let kin_fp = cint_fingerprint(mol.cint.integrate("int1e_kin", None, None).out.as_ref().unwrap());
    let nuc_fp = cint_fingerprint(mol.cint.integrate("int1e_nuc", None, None).out.as_ref().unwrap());
    let ecp_fp = cint_fingerprint(mol.cint.integrate_row_major("ECPscalar_iprinv", None, None).out.as_ref().unwrap());

    assert_eq!(nao, ref_nao);
    assert_eq!(nbas, ref_nbas);
    assert!(
        (kin_fp - ref_kin_fp).abs() < 1e-6 || (kin_fp / ref_kin_fp - 1.0).abs() < 1e-4,
        "kinetic fingerprint mismatch: got {kin_fp}, expected {ref_kin_fp}"
    );
    assert!(
        (nuc_fp - ref_nuc_fp).abs() < 1e-6 || (nuc_fp / ref_nuc_fp - 1.0).abs() < 1e-4,
        "nuclear fingerprint mismatch: got {nuc_fp}, expected {ref_nuc_fp}"
    );
    assert!(
        (ecp_fp - ref_ecp_fp).abs() < 1e-6 || (ecp_fp / ref_ecp_fp - 1.0).abs() < 1e-4,
        "ECP fingerprint mismatch: got {ecp_fp}, expected {ref_ecp_fp}"
    );
}

// simple case with sp shell
#[test]
fn simple_h2o_xyz_631g() {
    let toml_input = r#"
atom = """
O  0.0  0.0  1.0
H  1.2  0.5  1.3
H -0.7  1.2  0.1
"""

basis = "6-31G"
"#;
    let reference = (13, 9, 29.963738021815416, -74.71909470205469);
    check_cint_mol(toml_input, reference);
}

// zmat for three atoms
#[test]
fn simple_h2o_zmat_def2tzvp() {
    let toml_input = r#"
atom = "O; H 1 0.94; H 1 0.94 2 104.5"
basis = "def2-TZVP"
"#;
    let reference = (43, 19, 240.82414038242027, -455.36011189734484);
    check_cint_mol(toml_input, reference);
}

// large contraction in s orbital, cartesian
#[test]
fn simple_h2o_xyzcart_avqz() {
    let toml_input = r#"
atom = """
O  0.0  0.0  1.0
H  1.2  0.5  1.3
H -0.7  1.2  0.1
"""

basis = "aug-cc-pVQZ"
cart = true
"#;
    let reference = (215, 26, 117.86354585363418, -82.98711934741058);
    check_cint_mol(toml_input, reference);
}

// hybrid basis, zmat with dihedral
#[test]
fn difbas_ch4_zmat() {
    let toml_input = r#"
atom = """
C
H1 1 1.1
H  1 1.1 2 109.47
H2 1 1.1 2 109.47 3 120
H1 1 1.1 2 109.47 4 120
"""

basis = "custom"
cart = true

[basis-custom]
H1 = "aug-cc-pVQZ"
H2 = "def2-TZVP"
default = "6-31G"
"#;
    let reference = (127, 27, 26.883841020740512, -100.6569835983665);
    check_cint_mol(toml_input, reference);
}

// hybrid basis, zmat with dihedral, inline toml dict
#[test]
fn difbas_ch4_zmat_inline_dict() {
    let toml_input = r#"
atom = """
C
H1 1 1.1
H  1 1.1 2 109.47
H2 1 1.1 2 109.47 3 120
H1 1 1.1 2 109.47 4 120
"""

basis = { H1 = "aug-cc-pVQZ", H2 = "def2-TZVP", default = "6-31G" }
cart = true
"#;
    let reference = (127, 27, 26.883841020740512, -100.6569835983665);
    check_cint_mol(toml_input, reference);
}

// empty basis for some atoms
#[test]
fn empty_basis() {
    let toml_input = r#"
atom = """
C
H1 1 1.1
H  1 1.1 2 109.47
H2 1 1.1 2 109.47 3 120
H1 1 1.1 2 109.47 4 120
X-H2 1 1.1 2 50 4 60
"""

basis = "custom"
allow_empty_basis = true

[basis-custom]
H1 = "aug-cc-pVQZ"
H2 = "def2-TZVP"
C = "6-31G"
H = ""
"#;
    let reference = (113, 29, 25.5437545117713, -218.53479925755568);
    check_cint_mol(toml_input, reference);
}

// ECP tests

#[test]
fn simple_ecp() {
    let toml_input = r#"
atom = """
Sb        -1.33937843      0.44597852     -1.27279684
Sb1        1.33937843     -0.44597852     -1.27279684
C1        -1.40429524      1.10441871      0.83468205
C         -2.16210130     -1.56132398     -0.84717555
C2         2.16210130      1.56132398     -0.84717555
C          1.40429524     -1.10441871      0.83468205
H1        -0.69918639      1.91987631      1.00872018
H1        -1.16111477      0.29030616      1.51873028
H1        -2.40124532      1.47235562      1.08516843
H1        -2.02002046     -2.22909286     -1.69887295
H         -1.69052287     -2.01612927      0.02577778
H-1       -3.23450854     -1.49489801     -0.65423339
H-1        2.02002046      2.22909286     -1.69887295
H-1        3.23450854      1.49489801     -0.65423339
H          1.69052287      2.01612927      0.02577778
H          0.69918639     -1.91987631      1.00872018
H@2        2.40124532     -1.47235562      1.08516843
H@2        1.16111477     -0.29030616      1.51873028
"""

basis = "def2-TZVP"
ecp = "def2-TZVP"
"#;
    let reference = (296, 124, 46.36693654819834, -85.3713246949453, -382.5298577424865);
    check_cint_mol_ecp(toml_input, reference);
}

#[test]
fn hybrid_ecp() {
    let toml_input = r#"
atom = """
Sb        -1.33937843      0.44597852     -1.27279684
Sb1        1.33937843     -0.44597852     -1.27279684
C1        -1.40429524      1.10441871      0.83468205
C         -2.16210130     -1.56132398     -0.84717555
C2         2.16210130      1.56132398     -0.84717555
C          1.40429524     -1.10441871      0.83468205
H1        -0.69918639      1.91987631      1.00872018
H1        -1.16111477      0.29030616      1.51873028
H1        -2.40124532      1.47235562      1.08516843
H1        -2.02002046     -2.22909286     -1.69887295
H         -1.69052287     -2.01612927      0.02577778
H-1       -3.23450854     -1.49489801     -0.65423339
H-1        2.02002046      2.22909286     -1.69887295
H-1        3.23450854      1.49489801     -0.65423339
H          1.69052287      2.01612927      0.02577778
H          0.69918639     -1.91987631      1.00872018
H@2        2.40124532     -1.47235562      1.08516843
H@2        1.16111477     -0.29030616      1.51873028
"""

basis = "def2-TZVP"
ecp = "lanl2dz"
"#;
    let reference = (296, 124, 46.36693654819834, -101.63816733519229, -356.64926629582527);
    check_cint_mol_ecp(toml_input, reference);
}

#[test]
fn complicated_ecp() {
    let toml_input = r#"
atom = """
Sb        -1.33937843      0.44597852     -1.27279684
Sb1        1.33937843     -0.44597852     -1.27279684
C1        -1.40429524      1.10441871      0.83468205
C         -2.16210130     -1.56132398     -0.84717555
C2         2.16210130      1.56132398     -0.84717555
C          1.40429524     -1.10441871      0.83468205
H1        -0.69918639      1.91987631      1.00872018
H1        -1.16111477      0.29030616      1.51873028
H1        -2.40124532      1.47235562      1.08516843
H1        -2.02002046     -2.22909286     -1.69887295
H         -1.69052287     -2.01612927      0.02577778
H-1       -3.23450854     -1.49489801     -0.65423339
H-1        2.02002046      2.22909286     -1.69887295
H-1        3.23450854      1.49489801     -0.65423339
H          1.69052287      2.01612927      0.02577778
H          0.69918639     -1.91987631      1.00872018
H@2        2.40124532     -1.47235562      1.08516843
H@2        1.16111477     -0.29030616      1.51873028
"""

basis = "custom"
ecp = "custom"

[basis-custom]
H1 = "2ZaPa-NR"
"H@2" = "AHGBSP1-5"
C1 = "Ahlrichs VDZ"
C = "pcS-3"
Sb1 = "dyall-v2z"
default = "def2-SVP"

[ecp-custom]
C2 = "pSBKJC"
Sb1 = "lanl2dz"
"#;
    let reference = (499, 149, 56669463.51097445, -33390.20950360246, 241.4214981674055);
    check_cint_mol_ecp(toml_input, reference);
}

#[test]
fn ghost_heavy() {
    let toml_input = r#"
atom = """
Sb        -1.33937843      0.44597852     -1.27279684
Sb1        1.33937843     -0.44597852     -1.27279684
C1        -1.40429524      1.10441871      0.83468205
C         -2.16210130     -1.56132398     -0.84717555
C2         2.16210130      1.56132398     -0.84717555
C          1.40429524     -1.10441871      0.83468205
H1        -0.69918639      1.91987631      1.00872018
H1        -1.16111477      0.29030616      1.51873028
H1        -2.40124532      1.47235562      1.08516843
H1        -2.02002046     -2.22909286     -1.69887295
H         -1.69052287     -2.01612927      0.02577778
H-1       -3.23450854     -1.49489801     -0.65423339
H-1        2.02002046      2.22909286     -1.69887295
H-1        3.23450854      1.49489801     -0.65423339
H          1.69052287      2.01612927      0.02577778
H          0.69918639     -1.91987631      1.00872018
H@2        2.40124532     -1.47235562      1.08516843
H@2        1.16111477     -0.29030616      1.51873028
GHOST-Sb  -2.33937843      1.44597852     -1.27279684
GHOST-Sb1  2.33937843     -1.44597852     -1.27279684
"""

basis = "custom"
ecp = "custom"

[basis-custom]
H1 = "2ZaPa-NR"
"H@2" = "AHGBSP1-5"
C1 = "Ahlrichs VDZ"
C = "pcS-3"
Sb1 = "dyall-v2z"
default = "def2-SVP"

[ecp-custom]
C2 = "pSBKJC"
Sb1 = "lanl2dz"
"#;
    let reference = (646, 206, 64220631.18409411, 26146.394518924033, -41.935106631928406);
    check_cint_mol_ecp(toml_input, reference);
}
