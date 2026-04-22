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
use crate::parse::basis::{resolve_basis, BasisElement, BasisInput, BasisSpec};
use crate::prelude::*;

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
    pub basis_elements: BTreeMap<String, BasisElement>,
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
        let basis_elements = resolve_basis(&atoms, &self.basis, &self.ecp, self.ghost_ecp)?;

        // 3. Generate CInt arrays
        let cint_type = if self.cart { CIntType::Cartesian } else { CIntType::Spheric };
        let cint = generate_cint_arrays(&atoms, &basis_elements, cint_type)?;

        Ok(CIntMol { atoms, basis_elements, cint })
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

/// Generate CInt arrays from parsed atoms and basis data.
fn generate_cint_arrays(
    atoms: &[AtomInfo],
    basis_elements: &BTreeMap<String, BasisElement>,
    cint_type: CIntType,
) -> Result<CInt, CIntError> {
    let mut atm: Vec<[c_int; 6]> = Vec::new();
    let mut bas: Vec<[c_int; 8]> = Vec::new();
    let mut env: Vec<f64> = Vec::new();
    let mut ecpbas: Vec<[c_int; 8]> = Vec::new();

    // Initialize env with placeholder values
    // Index 0-19: reserved slots (set to 0)
    env.extend_from_slice(&[0.0; 20]);

    // Generate _atm and coordinates
    // Note: PySCF uses 4 slots per atom: x, y, z coordinates + zeta slot (for
    // nuclear model)
    for atom in atoms.iter() {
        // Get ECP electrons for this atom
        let ecp_core = basis_elements.get(&atom.label).and_then(|be| be.basis_data.ecp_electrons).unwrap_or(0);

        // Calculate effective charge (nuclear charge minus ECP core electrons)
        let effective_charge = atom.charge - ecp_core as f64;

        // Coordinates already in Bohr from AtomInfo
        // Layout: [x, y, z, zeta] - 4 slots per atom
        let coord_start = env.len() as c_int;
        env.push(atom.coords[0]);
        env.push(atom.coords[1]);
        env.push(atom.coords[2]);
        env.push(0.0); // zeta slot (unused for point nucleus)

        // ptr_zeta = coord_start + 3 (points to the slot after coordinates)
        let ptr_zeta = coord_start + 3;

        // atm slot: [charge, ptr_coord, nuc_mod, ptr_zeta, ptr_frac_charge, 0]
        atm.push([effective_charge as c_int, coord_start, 1, ptr_zeta, 0, 0]);
    }

    // Collect unique element symbols with their basis data
    // Process in order of appearance in atom list (not alphabetical)
    // PySCF stores basis data in the order elements appear in the basis dict
    // Use atom labels for lookup, but track unique elements by symbol
    let mut unique_elements: Vec<(String, String)> = Vec::new(); // (symbol, first_atom_label)
    for atom in atoms.iter() {
        if !unique_elements.iter().any(|(sym, _)| sym == &atom.symbol) {
            unique_elements.push((atom.symbol.clone(), atom.label.clone()));
        }
    }
    // Do NOT sort alphabetically - keep order of appearance in atom list

    // Pre-generate all basis data for unique elements
    let mut basis_cache: HashMap<String, Vec<(c_int, c_int, c_int, c_int)>> = HashMap::new();

    for (element_symbol, first_atom_label) in &unique_elements {
        let basis_elem =
            basis_elements.get(first_atom_label).ok_or_else(|| cint_error!(ParseError, "No basis for element {}", element_symbol))?;

        if basis_elem.basis_data.electron_shells.is_none() {
            continue;
        }

        let shells = basis_elem.basis_data.electron_shells.as_ref().unwrap();
        let mut indices: Vec<(c_int, c_int, c_int, c_int)> = Vec::new();

        for shell in shells {
            for &l in &shell.angular_momentum {
                let exponents: Vec<f64> = shell
                    .exponents
                    .iter()
                    .map(|s| s.parse::<f64>())
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| cint_error!(ParseError, "Failed to parse exponents: {}", e))?;

                let nprim = exponents.len() as c_int;

                let coeff_idx = shell.angular_momentum.iter().position(|&am| am == l).unwrap_or(0);
                let coefficients: Vec<f64> = shell
                    .coefficients
                    .get(coeff_idx)
                    .unwrap_or(&shell.coefficients[0])
                    .iter()
                    .map(|s| s.parse::<f64>())
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| cint_error!(ParseError, "Failed to parse coefficients: {}", e))?;

                let nctr = coefficients.len() / nprim as usize;
                let normalized_coeffs = normalize_shell(l, &exponents, &coefficients, cint_type);

                let exp_start = env.len() as c_int;
                env.extend(exponents.iter().cloned());

                let coeff_start = env.len() as c_int;
                env.extend(normalized_coeffs.iter().cloned());

                indices.push((exp_start, coeff_start, nprim, nctr as c_int));
            }
        }

        basis_cache.insert(element_symbol.clone(), indices);
    }

    // Generate _bas for each atom (using pre-generated basis data)
    for (idx, atom) in atoms.iter().enumerate() {
        let basis_elem = basis_elements.get(&atom.label).ok_or_else(|| cint_error!(ParseError, "No basis for atom {}", atom.label))?;

        if basis_elem.basis_data.electron_shells.is_none() {
            continue;
        }

        let shells = basis_elem.basis_data.electron_shells.as_ref().unwrap();
        let element_symbol = &atom.symbol;
        let cached_indices =
            basis_cache.get(element_symbol).ok_or_else(|| cint_error!(ParseError, "No cached basis for element {}", element_symbol))?;

        let mut shell_idx = 0;
        for shell in shells {
            for &l in &shell.angular_momentum {
                let (exp_start, coeff_start, nprim, nctr) = cached_indices[shell_idx];
                bas.push([idx as c_int, l as c_int, nprim, nctr, 0, exp_start, coeff_start, 0]);
                shell_idx += 1;
            }
        }
    }

    // Generate _ecpbas for atoms with ECP
    for (idx, atom) in atoms.iter().enumerate() {
        let ecp_core = basis_elements.get(&atom.label).and_then(|be| be.basis_data.ecp_electrons).unwrap_or(0);
        if ecp_core == 0 {
            continue;
        }

        let basis_elem = basis_elements.get(&atom.label).ok_or_else(|| cint_error!(ParseError, "No basis for atom {}", atom.label))?;

        if basis_elem.basis_data.ecp_potentials.is_none() {
            continue;
        }

        let potentials = basis_elem.basis_data.ecp_potentials.as_ref().unwrap();

        for potential in potentials {
            // Each potential can have multiple angular momentants
            for &l in &potential.angular_momentum {
                let exponents: Vec<f64> = potential
                    .gaussian_exponents
                    .iter()
                    .map(|s| s.parse::<f64>())
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| cint_error!(ParseError, "Failed to parse ECP exponents: {}", e))?;

                let r_exponents: Vec<i32> = potential.r_exponents.to_vec();
                let nexp = exponents.len() as c_int;

                // Get coefficients for this angular momentum
                let coeff_idx = potential.angular_momentum.iter().position(|&am| am == l).unwrap_or(0);
                let coefficients: Vec<f64> = potential
                    .coefficients
                    .get(coeff_idx)
                    .unwrap_or(&potential.coefficients[0])
                    .iter()
                    .map(|s| s.parse::<f64>())
                    .collect::<Result<Vec<_>, _>>()
                    .map_err(|e| cint_error!(ParseError, "Failed to parse ECP coefficients: {}", e))?;

                // Add exponents to env
                let exp_start = env.len() as c_int;
                env.extend(exponents.iter().cloned());

                // Add coefficients to env
                let coeff_start = env.len() as c_int;
                env.extend(coefficients.iter().cloned());

                // Max r_power for this potential
                let r_power = r_exponents.iter().max().copied().unwrap_or(0) as c_int;

                // ecpbas slot: [atom_id, l, nexp, r_power, so_type, ptr_exp, ptr_coeff, 0]
                ecpbas.push([idx as c_int, l as c_int, nexp, r_power, 0, exp_start, coeff_start, 0]);
            }
        }
    }

    let mut cint = CInt::new();
    cint.atm = atm;
    cint.bas = bas;
    cint.env = env;
    cint.ecpbas = ecpbas;
    cint.cint_type = cint_type;
    Ok(cint)
}

/// Normalize shell coefficients to match PySCF convention exactly.
fn normalize_shell(l: i32, exponents: &[f64], coefficients: &[f64], _cint_type: CIntType) -> Vec<f64> {
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
                let ee_jk = gaussian_int(l * 2 + 2, exponents[j] + exponents[k]);
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

/// Calculate gaussian_int(n, alpha) following PySCF's formula exactly.
/// gaussian_int(n, alpha) = gamma((n+1)/2) / (2 * alpha^((n+1)/2))
fn gaussian_int(n: i32, alpha: f64) -> f64 {
    let n1 = (n + 1) as f64 * 0.5;
    libm::tgamma(n1) / (2.0 * alpha.powf(n1))
}

/// Calculate GTO normalization factor following PySCF's formula.
/// gto_norm(l, exp) = 1/sqrt(gaussian_int(l*2+2, 2*exp))
fn gto_norm(l: i32, exp: f64) -> f64 {
    1.0 / gaussian_int(l * 2 + 2, 2.0 * exp).sqrt()
}
