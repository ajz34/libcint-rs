use libcint::{parse::basis::BasisSpec, prelude::*};

macro_rules! test_mol {
    (case: $case:ident; xyz: $xyz:expr; basis: $basis:expr, cart: $cart:expr, reference: ($ref_nao:expr, $ref_nbas:expr, $ref_kin_fp:expr)) => {
        #[test]
        fn $case() {
            let case_name = stringify!($case);
            check_cint_mol(case_name, $xyz, $basis, $cart, $ref_nao, $ref_nbas, $ref_kin_fp);
        }
    };
    (case: $case:ident; xyz: $xyz:expr; basis: $basis:expr, reference: ($ref_nao:expr, $ref_nbas:expr, $ref_kin_fp:expr)) => {
        #[test]
        fn $case() {
            let case_name = stringify!($case);
            check_cint_mol(case_name, $xyz, $basis, false, $ref_nao, $ref_nbas, $ref_kin_fp);
        }
    };
}

fn check_cint_mol(case_name: &str, xyz: &str, basis: impl Into<BasisSpec>, cart: bool, ref_nao: usize, ref_nbas: usize, ref_kin_fp: f64) {
    println!("Testing case: {case_name}");
    let mol_input = CIntMolInputBuilder::default().atom(xyz).basis(basis).cart(cart).build().unwrap();
    let mol = mol_input.create_mol();

    let nao = mol.cint.nao();
    let nbas = mol.cint.nbas();
    let kin_fp = cint_fingerprint(mol.cint.integrate("int1e_kin", None, None).out.as_ref().unwrap());

    assert_eq!(nao, ref_nao);
    assert_eq!(nbas, ref_nbas);
    assert!((kin_fp - ref_kin_fp).abs() < 1e-6, "kinetic fingerprint mismatch: got {kin_fp}, expected {ref_kin_fp}");
}

test_mol! {
    case: simple_h2o_xyz_631g;
    xyz: r"
        O  0.0  0.0  1.0
        H  1.2  0.5  1.3
        H -0.7  1.2  0.1
    ";
    basis: "6-31G",
    reference: (13, 9, 29.963738021815416)
}

test_mol! {
    case: simple_h2o_zmat_def2tzvp;
    xyz: r"
        O; H 1 0.94; H 1 0.94 2 104.5
    ";
    basis: "def2-TZVP",
    reference: (43, 19, 240.82414038242027)
}

test_mol! {
    case: difbas_h2o_xyzcart_avqz;
    xyz: r"
        O  0.0  0.0  1.0
        H  1.2  0.5  1.3
        H -0.7  1.2  0.1
    ";
    basis: "aug-cc-pVQZ",
    cart: true,
    reference: (215, 26, 117.86354585363418)
}

#[test]
fn playground() {
    use bse::prelude::*;

    let arg = BseGetBasisArgsBuilder::default().elements("O".to_string()).build().unwrap();
    // let basis = bse::get_formatted_basis("aug-cc-pVQZ", "g94", arg);
    let basis = bse::get_basis("aug-cc-pVQZ", arg);
    println!("Basis data:\n{:#?}", basis);
}
