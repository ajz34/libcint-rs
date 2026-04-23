use libcint::parse::basis::BasisSpec;
use libcint::prelude::*;

macro_rules! test_mol {
    (case: $case:ident; xyz: $xyz:expr; basis: $basis:expr, cart: $cart:expr, reference: ($ref_nao:expr, $ref_nbas:expr, $ref_kin_fp:expr, $ref_nuc_fp:expr)) => {
        #[test]
        fn $case() {
            let case_name = stringify!($case);
            check_cint_mol(case_name, $xyz, $basis, $cart, $ref_nao, $ref_nbas, $ref_kin_fp, $ref_nuc_fp);
        }
    };
}

fn check_cint_mol(
    case_name: &str,
    xyz: &str,
    basis: impl Into<BasisSpec>,
    cart: bool,
    ref_nao: usize,
    ref_nbas: usize,
    ref_kin_fp: f64,
    ref_nuc_fp: f64,
) {
    println!("Testing case: {case_name}");
    let mol_input = CIntMolInputBuilder::default().atom(xyz).basis(basis).cart(cart).allow_empty_basis(true).build().unwrap();
    let mol = mol_input.create_mol();

    let nao = mol.cint.nao();
    let nbas = mol.cint.nbas();
    let kin_fp = cint_fingerprint(mol.cint.integrate("int1e_kin", None, None).out.as_ref().unwrap());
    let nuc_fp = cint_fingerprint(mol.cint.integrate("int1e_nuc", None, None).out.as_ref().unwrap());

    assert_eq!(nao, ref_nao);
    assert_eq!(nbas, ref_nbas);
    assert!((kin_fp - ref_kin_fp).abs() < 1e-6, "kinetic fingerprint mismatch: got {kin_fp}, expected {ref_kin_fp}");
    assert!((nuc_fp - ref_nuc_fp).abs() < 1e-6, "nuclear fingerprint mismatch: got {nuc_fp}, expected {ref_nuc_fp}");
}

// simple case with sp shell
test_mol! {
    case: simple_h2o_xyz_631g;
    xyz: r"
        O  0.0  0.0  1.0
        H  1.2  0.5  1.3
        H -0.7  1.2  0.1
    ";
    basis: "6-31G",
    cart: false,
    reference: (13, 9, 29.963738021815416, -74.71909470205469)
}

// zmat for three atoms
test_mol! {
    case: simple_h2o_zmat_def2tzvp;
    xyz: r"
        O; H 1 0.94; H 1 0.94 2 104.5
    ";
    basis: "def2-TZVP",
    cart: false,
    reference: (43, 19, 240.82414038242027, -455.36011189734484)
}

// large contraction in s orbtial, cartesian
test_mol! {
    case: simple_h2o_xyzcart_avqz;
    xyz: r"
        O  0.0  0.0  1.0
        H  1.2  0.5  1.3
        H -0.7  1.2  0.1
    ";
    basis: "aug-cc-pVQZ",
    cart: true,
    reference: (215, 26, 117.86354585363418, -82.98711934741058)
}

// hybrid basis, zmat with dihedral
test_mol! {
    case: difbas_ch4_zmat;
    xyz: r"
        C
        H1 1 1.1
        H  1 1.1 2 109.47
        H2 1 1.1 2 109.47 3 120
        H1 1 1.1 2 109.47 4 120
    ";
    basis: vec![
        ("H1", "aug-cc-pVQZ"),
        ("H2", "def2-TZVP"),
        ("default", "6-31G"),
    ],
    cart: true,
    reference: (127, 27, 26.883841020740512, -100.6569835983665)
}

// hybrid basis, zmat with dihedral, hybrid formats
test_mol! {
    case: difbas_ch4_fmthybrid;
    xyz: r"
        C
        H1 1 1.1
        H  1 1.1 2 109.47
        H2 1 1.1 2 109.47 3 120
        H1 1 1.1 2 109.47 4 120
    ";
    basis: vec![
        ("H1", bse::get_formatted_basis("aug-cc-pVQZ", "nwchem", "")),
        ("H2", bse::get_formatted_basis("def2-TZVP", "cp2k", "")),
        ("default", bse::get_formatted_basis("6-31G", "g94", "")),
    ],
    cart: true,
    reference: (127, 27, 26.883841020740512, -100.6569835983665)
}

// hybrid basis, zmat with dihedral, hybrid formats
test_mol! {
    case: ghost_ch4_fmthybrid;
    xyz: r"
        C
        H1 1 1.1
        H  1 1.1 2 109.47
        H2 1 1.1 2 109.47 3 120
        H1 1 1.1 2 109.47 4 120
        GHOST-H2 1 1.1 2 50 4 60
    ";
    basis: vec![
        ("H1", bse::get_formatted_basis("aug-cc-pVQZ", "cp2k", "")),
        ("H2", bse::get_formatted_basis("def2-TZVP", "gamess_us", "")),
        ("default", bse::get_formatted_basis("6-31G", "crystal", "")),
    ],
    cart: false,
    reference: (115, 31, 22.917103194160873, -59.51548379533634)
}

// hybrid basis, zmat with dihedral, hybrid formats
test_mol! {
    case: ghost_ch4_uncontract_general;
    xyz: r"
        C
        H1 1 1.1
        H  1 1.1 2 109.47
        H2 1 1.1 2 109.47 3 120
        H1 1 1.1 2 109.47 4 120
        GHOST-H2 1 1.1 2 50 4 60
    ";
    basis: vec![
        ("H1", bse::get_formatted_basis("aug-cc-pVQZ", "g94", "")),
        ("H2", bse::get_formatted_basis("def2-TZVP", "turbomole", "")),
        ("default", bse::get_formatted_basis("6-31G", "cp2k", "")),
    ],
    cart: false,
    reference: (115, 43, 22.917103194160873, -59.51548379533635)
}

// hybrid basis, zmat with dihedral, hybrid formats
test_mol! {
    case: x_ch4_fmthybrid;
    xyz: r"
        C
        H1 1 1.1
        H  1 1.1 2 109.47
        H2 1 1.1 2 109.47 3 120
        H1 1 1.1 2 109.47 4 120
        X-H2 1 1.1 2 50 4 60
    ";
    basis: vec![
        ("H1", bse::get_formatted_basis("aug-cc-pVQZ", "nwchem", "")),
        ("H2", bse::get_formatted_basis("def2-TZVP", "nwchem", "")),
        ("DEFAULT", bse::get_formatted_basis("6-31G", "nwchem", "")),
    ],
    cart: true,
    reference: (133, 31, 6.762701948984031, 6.691179124574312)
}

// hybrid basis, zmat with dihedral, hybrid formats
test_mol! {
    case: empty_basis;
    xyz: r"
        C
        H1 1 1.1
        H  1 1.1 2 109.47
        H2 1 1.1 2 109.47 3 120
        H1 1 1.1 2 109.47 4 120
        X-H2 1 1.1 2 50 4 60
    ";
    basis: vec![
        ("H1", "aug-cc-pVQZ"),
        ("H2", "def2-TZVP"),
        ("C", "6-31G"),
        ("H", ""),
    ],
    cart: false,
    reference: (113, 29, 25.5437545117713, -218.53479925755568)
}

#[test]
fn playground() {
    use bse::prelude::*;

    let token = bse::get_formatted_basis("aug-cc-pvqz", "molcas", r#"elements = "O, H""#);
    let parsed = bse::read_formatted_basis_str(&token, "molcas");
    let token = bse::write_formatted_basis_str(&BseBasis::from_minimal(parsed), "nwchem", None);
    println!("Parsed basis:\n{:}", token);

    println!("======");

    let token = bse::get_formatted_basis("aug-cc-pvqz", "g94", r#"elements = "O, H""#);
    let parsed = bse::read_formatted_basis_str(&token, "g94");
    let token = bse::write_formatted_basis_str(&BseBasis::from_minimal(parsed), "nwchem", None);
    println!("Parsed basis:\n{:}", token);
}
