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
use crate::parse::basis::{resolve_basis, resolve_ecp, BasisElement, BasisSpec, EcpSpec};
use crate::prelude::*;
use std::collections::HashMap;

/// Builder for creating CIntMol molecule object.
///
/// Similar to PySCF's `gto.Mole`, this struct collects input parameters
/// and generates a `CInt` instance through the `build()` method.
#[derive(Debug, Clone, Builder, Default)]
#[builder(pattern = "owned", build_fn(error = "CIntError"), default)]
pub struct CIntMolInput {
    /// Atom coordinates (string).
    pub atom: String,
    /// Basis set specification.
    pub basis: BasisSpec,
    /// ECP specification (default: Auto for heavy elements).
    pub ecp: EcpSpec,
    /// Unit of coordinates (Angstrom or Bohr).
    pub unit: Unit,
    /// Use cartesian basis (default: spherical).
    pub cart: bool,
}

/// Built molecule object containing parsed atoms and CInt instance.
#[derive(Debug, Clone)]
pub struct CIntMol {
    /// Parsed atom information.
    pub atoms: Vec<AtomInfo>,
    /// Resolved basis elements for each atom.
    pub basis_elements: HashMap<usize, BasisElement>,
    /// ECP electrons for each atom.
    pub ecp_electrons: HashMap<usize, i32>,
    /// Generated CInt instance for integral calculations.
    pub cint: CInt,
    /// Whether basis is cartesian.
    pub cart: bool,
}

impl CIntMolInput {
    /// Create a new CIntMolInput builder.
    pub fn new() -> Self {
        Self::default()
    }

    /// Set atom coordinates from string.
    pub fn atom_str(self, s: &str) -> Self {
        Self { atom: s.to_string(), ..self }
    }

    /// Set atom coordinates from list.
    pub fn atom_list(self, list: Vec<String>) -> Self {
        Self { atom: list.join("\n"), ..self }
    }

    /// Set basis set name (uniform for all atoms).
    pub fn basis_str(self, name: &str) -> Self {
        Self { basis: BasisSpec::uniform(name), ..self }
    }

    /// Set basis set from per-element map.
    pub fn basis_map(self, map: HashMap<String, String>) -> Self {
        Self { basis: BasisSpec::per_element(map), ..self }
    }

    /// Set coordinate unit to Angstrom.
    pub fn angstrom(self) -> Self {
        Self { unit: Unit::Angstrom, ..self }
    }

    /// Set coordinate unit to Bohr.
    pub fn bohr(self) -> Self {
        Self { unit: Unit::Bohr, ..self }
    }

    /// Enable cartesian basis.
    pub fn cartesian(self) -> Self {
        Self { cart: true, ..self }
    }

    /// Enable spherical basis (default).
    pub fn spherical(self) -> Self {
        Self { cart: false, ..self }
    }

    /// Build the CIntMol object (fallible version).
    pub fn build_f(self) -> Result<CIntMol, CIntError> {
        // 1. Parse atoms
        let atoms = self.parse_atoms()?;

        if atoms.is_empty() {
            return cint_raise!(ParseError, "No atoms provided");
        }

        // 2. Resolve basis for each atom
        let basis_elements = resolve_basis(&atoms, &self.basis)?;

        // 3. Resolve ECP
        let ecp_electrons = resolve_ecp(&atoms, &basis_elements, &self.ecp)?;

        // 4. Generate CInt arrays
        let cint_type = if self.cart { CIntType::Cartesian } else { CIntType::Spheric };
        let (atm, bas, env, ecpbas) = generate_cint_arrays(&atoms, &basis_elements, &ecp_electrons, cint_type)?;

        // 5. Create CInt instance
        let cint = CInt { atm, bas, env, ecpbas, cint_type };

        Ok(CIntMol { atoms, basis_elements, ecp_electrons, cint, cart: self.cart })
    }

    /// Build the CIntMol object (infallible version, panics on error).
    pub fn build(self) -> CIntMol {
        self.build_f().cint_unwrap()
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
                        "Backtrace for atom parsing:\n- Atom string error: {:?}\n- Z-matrix error: {:?}",
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
    basis_elements: &HashMap<usize, BasisElement>,
    ecp_electrons: &HashMap<usize, i32>,
    cint_type: CIntType,
) -> Result<(Vec<[c_int; 6]>, Vec<[c_int; 8]>, Vec<f64>, Vec<[c_int; 8]>), CIntError> {
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
    for (idx, atom) in atoms.iter().enumerate() {
        // Get ECP electrons for this atom
        let ecp_core = ecp_electrons.get(&idx).copied().unwrap_or(0);

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
    // Process in order of atomic number (H=1, He=2, ...) to match PySCF's ordering
    let mut unique_elements: Vec<(String, usize)> = Vec::new(); // (symbol, first_atom_idx)
    for (idx, atom) in atoms.iter().enumerate() {
        if !unique_elements.iter().any(|(sym, _)| sym == &atom.symbol) {
            unique_elements.push((atom.symbol.clone(), idx));
        }
    }
    // Sort alphabetically by element symbol (H < O) - PySCF stores basis data in
    // alphabetical order
    unique_elements.sort_by(|(sym_a, _), (sym_b, _)| sym_a.cmp(sym_b));

    // Pre-generate all basis data for unique elements
    // This ensures basis data is stored in element order (H first, then O)
    let mut basis_cache: HashMap<String, Vec<(c_int, c_int, c_int, c_int)>> = HashMap::new();

    for (element_symbol, first_atom_idx) in &unique_elements {
        let basis_elem =
            basis_elements.get(first_atom_idx).ok_or_else(|| cint_error!(ParseError, "No basis for element {}", element_symbol))?;

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
        let basis_elem = basis_elements.get(&idx).ok_or_else(|| cint_error!(ParseError, "No basis for atom {}", idx))?;

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
    for (idx, _atom) in atoms.iter().enumerate() {
        let ecp_core = ecp_electrons.get(&idx).copied().unwrap_or(0);
        if ecp_core == 0 {
            continue;
        }

        let basis_elem = basis_elements.get(&idx).ok_or_else(|| cint_error!(ParseError, "No basis for atom {}", idx))?;

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

    Ok((atm, bas, env, ecpbas))
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
    gamma(n1) / (2.0 * alpha.powf(n1))
}

/// Calculate GTO normalization factor following PySCF's formula.
/// gto_norm(l, exp) = 1/sqrt(gaussian_int(l*2+2, 2*exp))
fn gto_norm(l: i32, exp: f64) -> f64 {
    1.0 / gaussian_int(l * 2 + 2, 2.0 * exp).sqrt()
}

/// Calculate gamma function using recursion from known values.
/// gamma(x) for x = n + 0.5 where n is integer:
///   gamma(0.5) = sqrt(pi)
///   gamma(1.5) = 0.5 * sqrt(pi)
///   gamma(2.5) = 1.5 * gamma(1.5) = 1.5 * 0.5 * sqrt(pi)
///   gamma(3.5) = 2.5 * gamma(2.5) = 2.5 * 1.5 * 0.5 * sqrt(pi)
/// For integer x: gamma(x) = (x-1)!
fn gamma(x: f64) -> f64 {
    use std::f64::consts::PI;
    let sqrt_pi = PI.sqrt();

    if x == 0.5 {
        sqrt_pi
    } else if x == 1.0 {
        1.0
    } else if x == 1.5 {
        sqrt_pi / 2.0
    } else if x == 2.0 {
        1.0
    } else if x == 2.5 {
        1.5 * sqrt_pi / 2.0
    } else if x == 3.0 {
        2.0
    } else if x == 3.5 {
        2.5 * 1.5 * sqrt_pi / 2.0
    }
    // = 3.75 * sqrt_pi / 2
    else if x == 4.0 {
        6.0
    } else if x == 4.5 {
        3.5 * 2.5 * 1.5 * sqrt_pi / 2.0
    }
    // = 13.125 * sqrt_pi / 2
    else if x == 5.0 {
        24.0
    } else if x == 5.5 {
        4.5 * 3.5 * 2.5 * 1.5 * sqrt_pi / 2.0
    }
    // = 59.0625 * sqrt_pi / 2
    else if x == 6.0 {
        120.0
    } else if x == 6.5 {
        5.5 * 4.5 * 3.5 * 2.5 * 1.5 * sqrt_pi / 2.0
    } else if x == 7.0 {
        720.0
    } else if x > 1.0 {
        (x - 1.0) * gamma(x - 1.0)
    } else {
        sqrt_pi / x
    } // For x < 0.5, use gamma(x) * x = gamma(x+1)
}

/// Convenience function to create molecule from atom string and basis name
/// (fallible). Named `M_f` following the fallible convention.
#[allow(non_snake_case)]
pub fn M_f(atom: &str, basis: &str) -> Result<CIntMol, CIntError> {
    CIntMolInput::new().atom_str(atom).basis_str(basis).angstrom().build_f()
}

/// Convenience function to create molecule from atom string and basis name
/// (infallible). Named `M` to match PySCF's convention.
#[allow(non_snake_case)]
pub fn M(atom: &str, basis: &str) -> CIntMol {
    M_f(atom, basis).cint_unwrap()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_molecule_creation() {
        let mol = M("H 0 0 0; H 0 0 0.74", "sto-3g");
        assert_eq!(mol.atoms.len(), 2);
        assert_eq!(mol.atoms[0].symbol, "H");
        assert_eq!(mol.atoms[1].symbol, "H");
    }

    #[test]
    fn test_water_sto3g() {
        let mol = CIntMolInput::new().atom_str("O 0 0 0; H 0 0.94 0; H 0.94 0 0").basis_str("sto-3g").angstrom().build();

        assert_eq!(mol.atoms.len(), 3);
        assert_eq!(mol.cint.atm.len(), 3);
        assert!(!mol.cint.bas.is_empty());
    }

    #[test]
    fn test_ghost_atom() {
        let mol = CIntMolInput::new().atom_str("H 0 0 0; GHOST-O 0 0 1.2").basis_str("sto-3g").angstrom().build();

        assert_eq!(mol.atoms.len(), 2);
        assert!(mol.atoms[1].is_ghost);
        assert_eq!(mol.atoms[1].charge, 0.0);
        // Ghost atom should have 0 in atm charge
        assert_eq!(mol.cint.atm[1][0], 0);
    }

    #[test]
    fn test_zmatrix_water() {
        let mol = CIntMolInput::new().atom_str("O\nH 1 0.94\nH 1 0.94 2 104.5").basis_str("sto-3g").angstrom().build();

        assert_eq!(mol.atoms.len(), 3);
        assert_eq!(mol.atoms[0].symbol, "O");
    }
}
