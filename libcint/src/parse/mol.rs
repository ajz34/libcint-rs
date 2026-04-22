//! CIntMol molecule object generation.
//!
//! This module provides `CIntMol` struct similar to PySCF's `gto.Mole`:
//! - Parse atom coordinates from string, list, or z-matrix format
//! - Resolve basis sets from BSE names, per-element dict, or custom format
//! - Generate `CInt` instance for integral calculations
//!
//! # Example
//!
//! ```text
//! use libcint::parse::mol::CIntMolInput;
//!
//! let mol = CIntMolInput::new()
//!     .atom_str("O; H 1 0.94; H 1 0.94 2 104.5")
//!     .basis_str("def2-TZVP")
//!     .angstrom()
//!     .build()
//!     .unwrap();
//!
//! // Now you can use mol.cint for integral calculations
//! let overlap = mol.cint.integrate("int1e_ovlp", None, None).unwrap();
//! ```

use crate::parse::atom::{parse_atom_string, parse_zmatrix, AtomInfo, Unit};
use crate::parse::basis::{resolve_basis, BasisInput, BasisSpec};
use crate::prelude::*;
use crate::util;

/// Builder for creating CIntMol molecule object.
///
/// Similar to PySCF's `gto.Mole`, this struct collects input parameters
/// and generates a `CInt` instance through the `build()` method.
#[non_exhaustive]
#[derive(Debug, Clone, Builder, Default)]
#[builder(pattern = "owned", build_fn(error = "CIntError"), default)]
pub struct CIntMolInput {
    /// Atom coordinates (string).
    #[builder(setter(into))]
    pub atom: String,
    /// Basis set specification.
    #[builder(setter(into))]
    pub basis: BasisSpec,
    /// ECP specification (optional).
    #[builder(setter(into), default = "BasisSpec::Uniform(BasisInput::None)")]
    pub ecp: BasisSpec,
    /// Unit of coordinates (Angstrom or Bohr).
    #[builder(setter(into), default = "Unit::Angstrom")]
    pub unit: Unit,
    /// Use cartesian basis (default: spherical).
    #[builder(default = "false")]
    pub cart: bool,
    /// Whether to assign ECP to ghost atoms (default: false).
    #[builder(default = "false")]
    pub ghost_ecp: bool,
}

/// Built molecule object containing parsed atoms and CInt instance.
#[derive(Debug, Clone)]
pub struct CIntMol {
    /// Parsed atom information.
    pub atoms: Vec<AtomInfo>,
    /// Resolved basis elements for each atom.
    ///
    /// This is label-basis mapping.
    pub basis_elements: IndexMap<String, BseBasisElement>,
    /// Map of atom labels to their corresponding basis element indices.
    pub atom_basis_map: Vec<String>,
    /// Generated CInt instance for integral calculations.
    pub cint: CInt,
}

impl CIntMolInput {
    /// Build the CIntMol object (fallible version).
    pub fn create_mol_f(self) -> Result<CIntMol, CIntError> {
        // 1. Parse atoms
        let atoms = self.parse_atoms()?;

        if atoms.is_empty() {
            return cint_raise!(ParseError, "No atoms provided");
        }

        // 2. Resolve basis for each atom
        let (basis_elements, atom_basis_map) = resolve_basis(&atoms, &self.basis, &self.ecp, self.ghost_ecp)?;

        // 3. Generate CInt arrays
        const PTR_ENV_START: usize = crate::ffi::cint_ffi::PTR_ENV_START as usize;
        let pre_env = vec![0.0; PTR_ENV_START]; // reserve env slots
        let cint_type = if self.cart { CIntType::Cartesian } else { CIntType::Spheric };
        let cint = make_env(&atoms, &basis_elements, &atom_basis_map, cint_type, pre_env)?;

        // 4. Add ECP data if any atom has ECP
        let has_ecp = basis_elements.values().any(|b| b.ecp_electrons.is_some() && b.ecp_potentials.is_some());
        let cint = if has_ecp { make_ecp_env(cint, &atoms, &basis_elements, &atom_basis_map)? } else { cint };

        Ok(CIntMol { atoms, basis_elements, atom_basis_map, cint })
    }

    /// Build the CIntMol object (infallible version, panics on error).
    pub fn create_mol(self) -> CIntMol {
        self.create_mol_f().cint_unwrap()
    }

    /// Parse atoms from input format.
    fn parse_atoms(&self) -> Result<Vec<AtomInfo>, CIntError> {
        match parse_atom_string(&self.atom, self.unit) {
            Ok(atoms) => Ok(atoms),
            Err(e_atom) => match parse_zmatrix(&self.atom, self.unit) {
                Ok(atoms) => Ok(atoms),
                Err(e_zmat) => {
                    // concat error messages from both attempts
                    let msg =
                        format!("Failed to parse atoms:\n- Atom string error: {:?}\n- Z-matrix error: {:?}", e_atom.kind, e_zmat.kind);
                    let trace_msg = format!(
                        "Backtrace for atom parsing:\n- Atom string error:\n{:?}\n- Z-matrix error:\n{:?}",
                        e_atom.backtrace, e_zmat.backtrace
                    );
                    cint_raise!(ParseError, "{msg}\n======\n{trace_msg}")
                },
            },
        }
    }
}

/// Add atm and env per atom.
///
/// atm fields: (charge, ptr_coord, nuc_mod, ptr_zeta, ptr_frac_charge, 0)
///
/// This function currently do not handle the following issues that PySCF do:
/// - zeta
/// - nuc_mod (assume NUC_POINT)
fn make_atm_env(atom: &AtomInfo, ptr: usize) -> ([i32; 6], [f64; 4]) {
    const CHARGE_OF: usize = crate::ffi::cint_ffi::CHARGE_OF as usize;
    const PTR_COORD: usize = crate::ffi::cint_ffi::PTR_COORD as usize;
    const POINT_NUC: i32 = crate::ffi::cint_ffi::POINT_NUC as i32;
    const NUC_MOD_OF: usize = crate::ffi::cint_ffi::NUC_MOD_OF as usize;
    const PTR_ZETA: usize = crate::ffi::cint_ffi::PTR_ZETA as usize;

    let nuc_charge = atom.charge;
    let mut atm = [0; 6];
    let mut env = [0.0; 4];

    atm[CHARGE_OF] = nuc_charge;
    atm[PTR_COORD] = ptr as i32;
    atm[NUC_MOD_OF] = POINT_NUC;
    atm[PTR_ZETA] = (ptr + 3) as i32;

    env[0..3].copy_from_slice(&atom.coords);

    (atm, env)
}

/// Add bas and env per atom.
///
/// bas fields: (atom_id, l, nprim, nctr, kappa, ptr_exp, ptr_coeff, 0)
///
/// This implementation will also absorb normalization into GTO contraction
/// coefficients, to match PySCF's convention exactly.
///
/// This function currently do not handle the following issues that PySCF do:
/// - kappa for relativistic basis (assume non-relativistic)
fn make_bas_env(basis: &BseBasisElement, atom_id: usize, ptr: usize) -> (Vec<[i32; 8]>, Vec<f64>) {
    let mut bas = Vec::new();
    let mut env = Vec::new();

    basis.electron_shells.as_ref().unwrap().iter().for_each(|shell| {
        shell.angular_momentum.iter().for_each(|&l| {
            let exponents: Vec<f64> = shell.exponents.iter().map(|s| s.parse::<f64>().unwrap()).collect();
            let coefficients: Vec<f64> = shell.coefficients.iter().flatten().map(|s| s.parse::<f64>().unwrap()).collect();

            let nprim = exponents.len() as i32;
            let nctr = coefficients.len() as i32 / nprim;

            // first sort exponents, and the corresponding coefficients, by exponent value
            // (decending)
            let mut sort_indices: Vec<usize> = (0..exponents.len()).collect();
            sort_indices.sort_by(|&i, &j| exponents[j].partial_cmp(&exponents[i]).unwrap());
            let exponents: Vec<f64> = sort_indices.iter().map(|&i| exponents[i]).collect();
            let coefficients: Vec<f64> = sort_indices
                .iter()
                .flat_map(|&i| {
                    let start = i * nctr as usize;
                    let end = start + nctr as usize;
                    coefficients[start..end].to_vec()
                })
                .collect();

            // Normalize coefficients to match PySCF convention
            let normalized_coeffs = normalize_shell(l, &exponents, &coefficients);

            // construct bas and env entries
            let ptr_exp = (ptr + env.len()) as i32;
            env.extend(exponents);
            let ptr_coeff = (ptr + env.len()) as i32;
            env.extend(normalized_coeffs);

            bas.push([atom_id as i32, l, nprim, nctr, 0, ptr_exp, ptr_coeff, 0]);
        });
    });

    (bas, env)
}

fn make_env(
    atoms: &[AtomInfo],
    basis_dict: &IndexMap<String, BseBasisElement>,
    atom_basis_map: &[String],
    cint_type: CIntType,
    pre_env: Vec<f64>,
) -> Result<CInt, CIntError> {
    const ATOM_OF: usize = crate::ffi::cint_ffi::ATOM_OF as usize;

    assert_eq!(atoms.len(), atom_basis_map.len(), "Number of atoms and atom_basis_map entries must match");

    let mut atm = vec![];
    let mut bas = vec![];
    let mut env = pre_env;

    // prepare atom and env
    for atom in atoms.iter() {
        let (atm_entry, env_entry) = make_atm_env(atom, env.len());
        atm.push(atm_entry);
        env.extend_from_slice(&env_entry);
    }

    // prepare bas and env
    // env directly write, bas will be paired to atom info
    let mut basdic = IndexMap::new();
    for (symb, basis_add) in basis_dict.iter() {
        let (bas_entries, env_entries) = make_bas_env(basis_add, 0, env.len());
        basdic.insert(symb.clone(), bas_entries);
        env.extend_from_slice(&env_entries);
    }

    for (ia, symb) in atom_basis_map.iter().enumerate() {
        let bas_entries = basdic.get(symb).ok_or_else(|| cint_error!(ParseError, "Basis entry must exist for atom {ia} symbol {symb}"))?;
        for bas_entry in bas_entries.iter() {
            let mut bas_entry = *bas_entry;
            bas_entry[ATOM_OF] = ia as i32;
            bas.push(bas_entry);
        }
    }

    let mut cint = CInt::new();
    cint.atm = atm;
    cint.bas = bas;
    cint.env = env;
    cint.cint_type = cint_type;
    Ok(cint)
}

fn make_ecp_env(
    mut cint: CInt,
    atoms: &[AtomInfo],
    basis_dict: &IndexMap<String, BseBasisElement>,
    atom_basis_map: &[String],
) -> Result<CInt, CIntError> {
    const CHARGE_OF: usize = crate::ffi::cint_ffi::CHARGE_OF as usize;
    const NUC_MOD_OF: usize = crate::ffi::cint_ffi::NUC_MOD_OF as usize;
    const ATOM_OF: usize = crate::ffi::cint_ffi::ATOM_OF as usize;
    const NUC_ECP: i32 = 4; // atoms with pseudo potential

    let mut ecpbas: Vec<[i32; 8]> = Vec::new();
    let ptr_env = cint.env.len();

    // Build ECP dictionary: label -> (nelec, ecpbas_entries)
    let mut ecpdic: IndexMap<String, (i32, Vec<[i32; 8]>)> = IndexMap::new();
    let mut ptr = ptr_env;

    for (symb, basis_add) in basis_dict.iter() {
        // Skip if no ECP potentials
        if basis_add.ecp_electrons.is_none() || basis_add.ecp_potentials.is_none() {
            continue;
        }

        let nelec = basis_add.ecp_electrons.unwrap();
        let potentials = basis_add.ecp_potentials.as_ref().unwrap();
        let mut ecp0: Vec<[i32; 8]> = Vec::new();

        for potential in potentials.iter() {
            // Get angular momentum (usually one value)
            for &l in potential.angular_momentum.iter() {
                // Parse exponents and coefficients
                let exponents: Vec<f64> = potential.gaussian_exponents.iter().map(|s| s.parse::<f64>().unwrap()).collect();
                let coefficients: Vec<f64> = potential.coefficients.iter().flatten().map(|s| s.parse::<f64>().unwrap()).collect();

                // Check for SO-ECP (if coefficients array has more than 1 inner vec)
                let has_so_ecp = potential.coefficients.len() > 1;

                // Sort by exponent (descending), similar to PySCF
                let mut sort_indices: Vec<usize> = (0..exponents.len()).collect();
                sort_indices.sort_by(|&i, &j| exponents[j].partial_cmp(&exponents[i]).unwrap());

                let sorted_exps: Vec<f64> = sort_indices.iter().map(|&i| exponents[i]).collect();
                let sorted_coeffs: Vec<f64> = sort_indices.iter().map(|&i| coefficients[i]).collect();

                let nexp = sorted_exps.len() as i32;

                // Determine radial power (r_exponents should be same for all terms in most
                // cases)
                let rorder = if potential.r_exponents.is_empty() { 0 } else { potential.r_exponents[0] };

                // Add exponents to env
                let ptr_exp = ptr as i32;
                cint.env.extend(&sorted_exps);
                ptr += sorted_exps.len();

                // Add coefficients to env
                let ptr_coeff = ptr as i32;
                cint.env.extend(&sorted_coeffs);
                ptr += sorted_coeffs.len();

                // Create ecpbas entry: [atom_id, l, nexp, rorder, so_type, ptr_exp, ptr_coeff,
                // 0] atom_id will be filled later when iterating through atoms
                // so_type = 0 for scalar ECP
                ecp0.push([0, l, nexp, rorder, 0, ptr_exp, ptr_coeff, 0]);

                // Handle SO-ECP: add second ecpbas entry for SO part
                if has_so_ecp {
                    // For SO-ECP, coefficients has multiple inner vecs
                    // The first inner vec is scalar, subsequent are SO
                    let so_coeffs: Vec<f64> = potential.coefficients.iter().skip(1).flatten().map(|s| s.parse::<f64>().unwrap()).collect();
                    let sorted_so_coeffs: Vec<f64> = sort_indices.iter().map(|&i| so_coeffs[i]).collect();

                    let ptr_so_coeff = ptr as i32;
                    cint.env.extend(&sorted_so_coeffs);
                    ptr += sorted_so_coeffs.len();

                    // so_type = 1 for SO-ECP
                    ecp0.push([0, l, nexp, rorder, 1, ptr_exp, ptr_so_coeff, 0]);
                }
            }
        }

        ecpdic.insert(symb.clone(), (nelec, ecp0));
    }

    // Apply ECP to atoms
    if !ecpdic.is_empty() {
        for (ia, symb) in atom_basis_map.iter().enumerate() {
            // Check if this atom's basis has ECP
            if let Some((nelec, ecp_entries)) = ecpdic.get(symb) {
                // Modify atm entry
                let original_charge = atoms[ia].charge;
                cint.atm[ia][CHARGE_OF] = original_charge - *nelec;
                cint.atm[ia][NUC_MOD_OF] = NUC_ECP;

                // Add ecpbas entries with correct atom_id
                for ecp_entry in ecp_entries.iter() {
                    let mut ecp_entry = *ecp_entry;
                    ecp_entry[ATOM_OF] = ia as i32;
                    ecpbas.push(ecp_entry);
                }
            }
        }
    }

    cint.ecpbas = ecpbas;
    Ok(cint)
}

/// Normalize shell coefficients to match PySCF convention exactly.
fn normalize_shell(l: i32, exponents: &[f64], coefficients: &[f64]) -> Vec<f64> {
    let nprim = exponents.len();
    let nctr = coefficients.len() / nprim;

    // Step 1: Primitive normalization using PySCF's formula
    // gto_norm(l, exp) = 1/sqrt(gaussian_int(l*2+2, 2*exp))
    let prim_norms: Vec<f64> = exponents.iter().map(|&exp| gto_norm(l, exp)).collect();

    // Apply primitive normalization to raw coefficients
    // cs = einsum('pi,p->pi', cs, gto_norm)
    // For coefficients in contraction-first layout: coeff[i*nprim + j] for contr i,
    // prim j
    let prim_normalized: Vec<f64> = (0..nprim * nctr)
        .map(|idx| {
            let i = idx % nctr; // contraction index
            let j = idx / nctr; // primitive index
            coefficients[i * nprim + j] * prim_norms[j]
        })
        .collect();

    // Step 2: Contraction normalization
    // ee[i,j] = gaussian_int(l*2+2, exp_i + exp_j)
    // s1 = 1/sqrt(einsum('pi,pq,qi->i', cs, ee, cs))
    // cs = einsum('pi,i->pi', cs, s1)
    let mut normalized: Vec<f64> = Vec::with_capacity(coefficients.len());

    for i in 0..nctr {
        // Calculate contraction integral: sum over j,k: cs[j,i] * ee[j,k] * cs[k,i]
        // Using contraction-first layout: prim_normalized[i*nprim + j]
        let mut contr_int = 0.0;
        for j in 0..nprim {
            for k in 0..nprim {
                let cj = prim_normalized[i * nprim + j];
                let ck = prim_normalized[i * nprim + k];
                let ee_jk = util::gaussian_int((l * 2 + 2) as f64, exponents[j] + exponents[k]);
                contr_int += cj * ck * ee_jk;
            }
        }
        let contr_norm_factor = 1.0 / contr_int.sqrt();

        // Normalize this contraction
        for j in 0..nprim {
            normalized.push(prim_normalized[i * nprim + j] * contr_norm_factor);
        }
    }

    normalized
}

/// Calculate GTO normalization factor following PySCF's formula.
/// gto_norm(l, exp) = 1/sqrt(gaussian_int(l*2+2, 2*exp))
fn gto_norm(l: i32, exp: f64) -> f64 {
    1.0 / util::gaussian_int((l * 2 + 2) as f64, 2.0 * exp).sqrt()
}
