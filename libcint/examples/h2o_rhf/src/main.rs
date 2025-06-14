use core::f64;
use libcint::prelude::*;
use rstsr::prelude::*;

pub type Tsr = Tensor<f64, DeviceBLAS>;
pub type TsrView<'a> = TensorView<'a, f64, DeviceBLAS>;

/// Obtain integrals (in row-major, same to PySCF but reverse of libcint)
pub fn intor_row_major(cint_data: &CInt, intor: &str) -> Tsr {
    let device = DeviceBLAS::default(); // use up all rayon available threads for tensor operations
    let (out, shape) = cint_data.integrate(intor, None, None).into(); // intor, "s1", full_shls_slice
    rt::asarray((out, shape.f(), &device)).into_reverse_axes() // row-major by transposition of col-major shape
}

fn main() {
    let device = DeviceBLAS::default();
    let cint_data = init_h2o_def2_tzvp();

    /* #region rhf total energy algorithms */
    let atom_coords = {
        let coords = cint_data.atom_coords();
        let coords = coords.into_iter().flatten().collect::<Vec<f64>>();
        rt::asarray((coords, &device)).into_shape((-1, 3))
    };
    let atom_charges = rt::asarray((cint_data.atom_charges(), &device));
    let mut dist = rt::sci::cdist((atom_coords.view(), atom_coords.view()));
    dist.diagonal_mut(None).fill(f64::INFINITY);
    let eng_nuc = 0.5 * (&atom_charges * atom_charges.i((.., None)) / dist).sum();
    println!("Nuclear repulsion energy: {eng_nuc}");

    let hcore = intor_row_major(&cint_data, "int1e_kin") + intor_row_major(&cint_data, "int1e_nuc");
    let ovlp = intor_row_major(&cint_data, "int1e_ovlp");
    let int2e = intor_row_major(&cint_data, "int2e");
    let nocc = 5; // hardcoded for H2O, 5 occupied orbitals

    let mut dm = ovlp.zeros_like();
    for _ in 0..40 {
        // hardcoded SCF iterations
        let fock = &hcore + ((1.0_f64 * &int2e - 0.5_f64 * int2e.swapaxes(1, 2)) * &dm).sum_axes([-1, -2]);
        let (_, c) = rt::linalg::eigh((&fock, &ovlp)).into();
        dm = 2.0_f64 * c.i((.., ..nocc)) % c.i((.., ..nocc)).t();
    }
    let eng_scratch = &hcore + ((0.5_f64 * &int2e - 0.25_f64 * int2e.swapaxes(1, 2)) * &dm).sum_axes([-1, -2]);
    let eng_elec = (dm * &eng_scratch).sum();
    println!("Total elec energy: {eng_elec}");
    println!("Total RHF energy: {}", eng_nuc + eng_elec);
    /* #endregion */
}
