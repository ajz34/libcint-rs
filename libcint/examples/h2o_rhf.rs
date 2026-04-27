//! H2O RHF calculation example using libcint.
//!
//! This example demonstrates a simple RHF calculation on H2O/def2-TZVP.

use libcint::prelude::*;
use rstsr::prelude::*;

pub type Tsr = Tensor<f64, DeviceCpu>;
pub type TsrView<'a> = TensorView<'a, f64, DeviceCpu>;

/// Obtain integrals (in row-major, same to PySCF but reverse of libcint)
pub fn intor_row_major(cint: &CInt, intor: &str) -> Tsr {
    // use up all rayon available threads for tensor operations
    let device = DeviceCpu::default();
    // intor, "s1", full_shls_slice
    let (out, shape) = cint.integrate_row_major(intor, None, None).into();
    // row-major by transposition of col-major shape
    rt::asarray((out, shape.c(), &device))
}

fn main() {
    let device = DeviceCpu::default();
    // Assuming H2O/def2-TZVP data for `CInt` has been prepared

    let mol_input = r#"
        atom = "O; H 1 0.94; H 1 0.94 2 104.5"
        basis = "def2-TZVP"
    "#;
    let cint_mol = CIntMol::from_toml(mol_input);
    let cint = cint_mol.cint;

    /* #region rhf total energy algorithms */
    let atom_coords = {
        let coords = cint.atom_coords();
        let coords = coords.into_iter().flatten().collect::<Vec<f64>>();
        rt::asarray((coords, &device)).into_shape((-1, 3))
    };
    let atom_charges = rt::asarray((cint.atom_charges(), &device));
    let mut dist = rt::sci::cdist((atom_coords.view(), atom_coords.view()));
    dist.diagonal_mut(None).fill(f64::INFINITY);
    let eng_nuc = 0.5 * (&atom_charges * atom_charges.i((.., None)) / dist).sum();
    println!("Nuclear repulsion energy: {eng_nuc}");

    let hcore = intor_row_major(&cint, "int1e_kin") + intor_row_major(&cint, "int1e_nuc");
    let ovlp = intor_row_major(&cint, "int1e_ovlp");
    let int2e = intor_row_major(&cint, "int2e");
    let nocc = 5; // hardcoded for H2O, 5 occupied orbitals

    let mut dm = ovlp.zeros_like();
    for _idx_iter in 0..40 {
        // hardcoded SCF iterations
        let fock = &hcore + ((1.0_f64 * &int2e - 0.5_f64 * int2e.swapaxes(1, 2)) * &dm).sum_axes([-1, -2]);
        let (_mo_energy, mo_coeff) = rt::linalg::eigh((&fock, &ovlp)).into();
        dm = 2.0_f64 * mo_coeff.i((.., ..nocc)) % mo_coeff.i((.., ..nocc)).t();
    }
    let eng_scratch = &hcore + ((0.5_f64 * &int2e - 0.25_f64 * int2e.swapaxes(1, 2)) * &dm).sum_axes([-1, -2]);
    let eng_elec = (dm * &eng_scratch).sum();
    println!("Total elec energy: {eng_elec}");
    println!("Total RHF energy: {}", eng_nuc + eng_elec);
    /* #endregion */
    assert!((eng_nuc + eng_elec - -76.05945519696209).abs() < 1e-8)
}
