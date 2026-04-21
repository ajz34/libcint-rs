//! Basis set resolution and handling.
//!
//! This module provides flexible basis set specification:
//! - Uniform basis: single BSE name for all atoms (e.g., "def2-TZVP")
//! - Per-element basis: dict by element symbol (e.g., {"H": "sto-3g", "O":
//!   "cc-pvdz"})
//! - Per-atom basis: dict by atom label (e.g., {"H1": "sto-3g", "H2":
//!   "cc-pvdz"})
//! - Custom basis: BseBasisElement format
//!
//! Supports ghost atoms and automatic ECP detection for heavy elements.

use std::collections::{HashMap, HashSet};

use crate::parse::atom::AtomInfo;
use crate::prelude::*;
use bse::lut::element_Z_from_sym;
use bse::prelude::*;

/// Specification for basis set assignment.
#[derive(Debug, Clone, PartialEq)]
pub enum BasisSpec {
    /// Single BSE name for all atoms (e.g., "def2-TZVP").
    Uniform(String),
    /// Per-element basis (e.g., {"H": "sto-3g", "O": "cc-pvdz"}).
    PerElement(HashMap<String, String>),
    /// Per-atom basis by label (e.g., {"H1": "sto-3g", "H2": "cc-pvdz"}).
    PerAtom(HashMap<String, String>),
    /// Custom basis in BseBasisElement format.
    Custom(HashMap<String, BseBasisElement>),
}

impl Default for BasisSpec {
    fn default() -> Self {
        BasisSpec::Uniform("sto-3g".to_string())
    }
}

impl BasisSpec {
    /// Create a uniform basis specification.
    pub fn uniform(name: &str) -> Self {
        BasisSpec::Uniform(name.to_string())
    }

    /// Create a per-element basis specification.
    pub fn per_element(specs: HashMap<String, String>) -> Self {
        BasisSpec::PerElement(specs)
    }

    /// Create a per-atom basis specification.
    pub fn per_atom(specs: HashMap<String, String>) -> Self {
        BasisSpec::PerAtom(specs)
    }

    /// Create from a single string (uniform basis).
    pub fn from_string(s: &str) -> Self {
        BasisSpec::Uniform(s.to_string())
    }

    /// Create from a HashMap of element symbols to basis names.
    pub fn from_map(map: HashMap<String, String>) -> Self {
        BasisSpec::PerElement(map)
    }
}

/// Resolved basis element for a specific atom.
#[derive(Debug, Clone, PartialEq)]
pub struct BasisElement {
    /// Element symbol (e.g., "H", "O").
    pub symbol: String,
    /// Atomic number.
    pub charge: f64,
    /// Basis set data from BSE.
    pub basis_data: BseBasisElement,
    /// Whether this atom has ECP.
    pub has_ecp: bool,
    /// Number of core electrons replaced by ECP.
    pub ecp_electrons: i32,
}

/// Resolve basis sets for a list of atoms.
///
/// # Arguments
/// * `atoms` - List of parsed atom information
/// * `basis_spec` - Basis set specification
///
/// # Returns
/// HashMap mapping atom index to resolved BasisElement.
pub fn resolve_basis(atoms: &[AtomInfo], basis_spec: &BasisSpec) -> Result<HashMap<usize, BasisElement>, CIntError> {
    let mut result: HashMap<usize, BasisElement> = HashMap::new();

    for (idx, atom) in atoms.iter().enumerate() {
        let basis_data = resolve_basis_for_atom(atom, idx, basis_spec)?;

        let has_ecp = basis_data.ecp_electrons.is_some() && basis_data.ecp_electrons.unwrap() > 0;
        let ecp_electrons = basis_data.ecp_electrons.unwrap_or(0);

        result.insert(idx, BasisElement { symbol: atom.symbol.clone(), charge: atom.charge, basis_data, has_ecp, ecp_electrons });
    }

    Ok(result)
}

/// Resolve basis set for a single atom.
fn resolve_basis_for_atom(atom: &AtomInfo, _idx: usize, spec: &BasisSpec) -> Result<BseBasisElement, CIntError> {
    // Ghost atoms can have basis if explicitly specified in PerAtom or PerElement
    // map Otherwise they get empty basis
    let skip_ghost_default = atom.is_ghost;

    match spec {
        BasisSpec::Uniform(name) => {
            // Ghost atoms get empty basis with uniform spec
            if skip_ghost_default {
                return Ok(BseBasisElement::default());
            }
            fetch_basis_element(name, &atom.symbol)
        },
        BasisSpec::PerElement(map) => {
            // Try label match first (for ghost-O, X_H, etc.)
            if let Some(name) = map.get(&atom.label) {
                fetch_basis_element(name, &atom.symbol)
            } else {
                // Try element symbol match
                let key = atom.symbol.to_uppercase();
                if let Some(name) = map.get(&key) {
                    // Ghost atoms get empty basis unless explicitly matched
                    if skip_ghost_default {
                        return Ok(BseBasisElement::default());
                    }
                    fetch_basis_element(name, &atom.symbol)
                } else if let Some(name) = map.get(&atom.symbol) {
                    if skip_ghost_default {
                        return Ok(BseBasisElement::default());
                    }
                    fetch_basis_element(name, &atom.symbol)
                } else {
                    cint_raise!(ParseError, "No basis specified for element '{}' in per-element basis map", atom.symbol)
                }
            }
        },
        BasisSpec::PerAtom(map) => {
            // Try exact label match first (H1, H2, ghost-O, etc.)
            if let Some(name) = map.get(&atom.label) {
                fetch_basis_element(name, &atom.symbol)
            } else if let Some(name) = map.get(&atom.symbol) {
                // Ghost atoms get empty basis unless explicitly matched
                if skip_ghost_default {
                    return Ok(BseBasisElement::default());
                }
                fetch_basis_element(name, &atom.symbol)
            } else {
                // Fallback to element symbol
                cint_raise!(ParseError, "No basis specified for atom '{}' (symbol '{}') in per-atom basis map", atom.label, atom.symbol)
            }
        },
        BasisSpec::Custom(map) => {
            // Try element symbol match
            let key = atom.symbol.to_uppercase();
            if let Some(elem) = map.get(&key) {
                Ok(elem.clone())
            } else if let Some(elem) = map.get(&atom.symbol) {
                Ok(elem.clone())
            } else {
                cint_raise!(ParseError, "No custom basis specified for element '{}'", atom.symbol)
            }
        },
    }
}

/// Fetch basis element data from BSE library.
fn fetch_basis_element(basis_name: &str, element_symbol: &str) -> Result<BseBasisElement, CIntError> {
    // Get atomic number from element symbol
    let z = element_Z_from_sym(element_symbol).ok_or_else(|| cint_error!(ParseError, "Unknown element symbol: {}", element_symbol))?;

    // Fetch basis set from BSE
    let args = BseGetBasisArgsBuilder::default()
        .elements(z.to_string())
        .build()
        .map_err(|e| cint_error!(ParseError, "Failed to build BSE args: {}", e))?;

    let basis = get_basis_f(basis_name, args)
        .map_err(|e| cint_error!(ParseError, "Failed to fetch basis '{}' for element '{}': {}", basis_name, element_symbol, e))?;

    // Extract element data
    let elem_key = z.to_string();
    basis
        .elements
        .get(&elem_key)
        .cloned()
        .ok_or_else(|| cint_error!(ParseError, "Element {} not found in basis '{}'", element_symbol, basis_name))
}

/// Check if a basis set provides ECP for a given element.
pub fn has_ecp_for_element(basis_name: &str, element_symbol: &str) -> Result<bool, CIntError> {
    let elem = fetch_basis_element(basis_name, element_symbol)?;
    Ok(elem.ecp_electrons.is_some() && elem.ecp_electrons.unwrap() > 0)
}

/// Get ECP electrons for a basis set and element.
pub fn get_ecp_electrons(basis_name: &str, element_symbol: &str) -> Result<i32, CIntError> {
    let elem = fetch_basis_element(basis_name, element_symbol)?;
    Ok(elem.ecp_electrons.unwrap_or(0))
}

/// ECP specification options.
#[derive(Default, Debug, Clone, PartialEq)]
pub enum EcpSpec {
    /// Use ECP automatically for heavy elements (default).
    #[default]
    Auto,
    /// Specify ECP basis by element.
    PerElement(HashMap<String, String>),
    /// Explicitly disable ECP for certain elements.
    Disabled(HashSet<String>),
    /// No ECP at all.
    None,
}

/// Threshold atomic number for automatic ECP detection.
/// Elements with Z > ECP_AUTO_THRESHOLD will use ECP if available.
pub const ECP_AUTO_THRESHOLD: i32 = 36; // Start from Rb (Z=37)

/// Resolve ECP for atoms based on specification.
///
/// # Arguments
/// * `atoms` - List of parsed atom information
/// * `basis_elements` - Already resolved basis elements
/// * `ecp_spec` - ECP specification
///
/// # Returns
/// HashMap mapping atom index to ECP electrons count.
pub fn resolve_ecp(
    atoms: &[AtomInfo],
    basis_elements: &HashMap<usize, BasisElement>,
    ecp_spec: &EcpSpec,
) -> Result<HashMap<usize, i32>, CIntError> {
    let mut result: HashMap<usize, i32> = HashMap::new();

    for (idx, atom) in atoms.iter().enumerate() {
        if atom.is_ghost {
            result.insert(idx, 0);
            continue;
        }

        let basis_elem = basis_elements.get(&idx).ok_or_else(|| cint_error!(ParseError, "Basis not resolved for atom {}", idx))?;

        let ecp_electrons = match ecp_spec {
            EcpSpec::Auto => {
                // Use ECP for heavy elements if available in basis
                if atom.charge > ECP_AUTO_THRESHOLD as f64 && basis_elem.has_ecp {
                    basis_elem.ecp_electrons
                } else {
                    0
                }
            },
            EcpSpec::PerElement(map) => {
                // Check if ECP specified for this element
                let key = atom.symbol.to_uppercase();
                if let Some(ecp_name) = map.get(&key) {
                    get_ecp_electrons(ecp_name, &atom.symbol)?
                } else if let Some(ecp_name) = map.get(&atom.symbol) {
                    get_ecp_electrons(ecp_name, &atom.symbol)?
                } else {
                    // Fall back to basis ECP if available
                    if basis_elem.has_ecp {
                        basis_elem.ecp_electrons
                    } else {
                        0
                    }
                }
            },
            EcpSpec::Disabled(disabled_set) => {
                // Check if element is in disabled set
                let key = atom.symbol.to_uppercase();
                if disabled_set.contains(&key) || disabled_set.contains(&atom.symbol) {
                    0
                } else {
                    // Use basis ECP if available
                    if basis_elem.has_ecp {
                        basis_elem.ecp_electrons
                    } else {
                        0
                    }
                }
            },
            EcpSpec::None => 0,
        };

        result.insert(idx, ecp_electrons);
    }

    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse::atom::{parse_atom_string, Unit};

    #[test]
    fn test_basis_spec_uniform() {
        let spec = BasisSpec::uniform("def2-TZVP");
        assert_eq!(spec, BasisSpec::Uniform("def2-TZVP".to_string()));
    }

    #[test]
    fn test_basis_spec_per_element() {
        let mut map = HashMap::new();
        map.insert("H".to_string(), "sto-3g".to_string());
        map.insert("O".to_string(), "cc-pvdz".to_string());
        let spec = BasisSpec::per_element(map);
        assert!(matches!(spec, BasisSpec::PerElement(_)));
    }

    #[test]
    fn test_resolve_basis_uniform() {
        let atoms = parse_atom_string("H 0 0 0; O 0 0 1.2", Unit::Angstrom).unwrap();
        let spec = BasisSpec::uniform("sto-3g");
        let basis = resolve_basis(&atoms, &spec).unwrap();

        assert_eq!(basis.len(), 2);
        assert!(basis.contains_key(&0));
        assert!(basis.contains_key(&1));
    }

    #[test]
    fn test_resolve_basis_ghost() {
        let atoms = parse_atom_string("H 0 0 0; GHOST-O 0 0 1.2", Unit::Angstrom).unwrap();
        let spec = BasisSpec::uniform("sto-3g");
        let basis = resolve_basis(&atoms, &spec).unwrap();

        // Ghost atom has empty basis
        assert!(basis.get(&1).unwrap().basis_data.electron_shells.is_none());
    }

    #[test]
    fn test_fetch_basis_element() {
        let elem = fetch_basis_element("sto-3g", "H").unwrap();
        assert!(elem.electron_shells.is_some());
        let shells = elem.electron_shells.unwrap();
        assert!(!shells.is_empty());
    }

    #[test]
    fn test_has_ecp_for_element() {
        // def2-TZVP has ECP for heavy elements like Au
        let has_ecp = has_ecp_for_element("def2-TZVP", "Au").unwrap();
        assert!(has_ecp);

        // def2-TZVP has no ECP for light elements like H
        let has_ecp = has_ecp_for_element("def2-TZVP", "H").unwrap();
        assert!(!has_ecp);
    }
}
