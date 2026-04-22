//! Basis set resolution and handling.

use super::atom;
use crate::parse::atom::AtomInfo;
use crate::prelude::*;
use bse::prelude::*;

/* #region BasisInput */

#[derive(Debug, Clone, PartialEq)]
pub enum BasisInput {
    None,
    String(String),
    Element(Box<BseBasisElement>),
    Basis(Box<BseBasis>),
}

impl From<String> for BasisInput {
    fn from(s: String) -> Self {
        BasisInput::String(s)
    }
}

impl From<&str> for BasisInput {
    fn from(s: &str) -> Self {
        BasisInput::String(s.to_string())
    }
}

impl From<BseBasisElement> for BasisInput {
    fn from(elem: BseBasisElement) -> Self {
        BasisInput::Element(Box::new(elem))
    }
}

impl From<BseBasis> for BasisInput {
    fn from(basis: BseBasis) -> Self {
        BasisInput::Basis(Box::new(basis))
    }
}

impl From<Option<String>> for BasisInput {
    fn from(opt: Option<String>) -> Self {
        match opt {
            Some(s) => BasisInput::String(s),
            None => BasisInput::None,
        }
    }
}

/* #endregion */

/* #region BasisSpec */

#[derive(Debug, Clone, PartialEq)]
pub enum BasisSpec {
    Uniform(BasisInput),
    Dict(BTreeMap<String, BasisInput>),
    List(Vec<BasisInput>),
}

impl Default for BasisSpec {
    fn default() -> Self {
        BasisSpec::Uniform(BasisInput::String("sto-3g".to_string()))
    }
}

#[allow(clippy::useless_conversion)]
#[duplicate_item(TY; [String]; [&str]; [Option<String>]; [BseBasisElement]; [BseBasis];)]
impl From<TY> for BasisSpec {
    fn from(s: TY) -> Self {
        BasisSpec::Uniform(s.into())
    }
}

#[duplicate_item(TY; [HashMap<String, T>]; [BTreeMap<String, T>];)]
impl<T> From<TY> for BasisSpec
where
    T: Into<BasisInput>,
{
    fn from(map: TY) -> Self {
        let dict = map.into_iter().map(|(k, v)| (k, v.into())).collect();
        BasisSpec::Dict(dict)
    }
}

impl<T> From<Vec<T>> for BasisSpec
where
    T: Into<BasisInput>,
{
    fn from(list: Vec<T>) -> Self {
        let vec = list.into_iter().map(|v| v.into()).collect();
        BasisSpec::List(vec)
    }
}

/* #endregion */

/// Resolved basis element for a specific atom.
#[derive(Debug, Clone, PartialEq)]
pub struct BasisElement {
    /// Element symbol (e.g., "H", "O").
    pub symbol: String,
    /// Atomic number.
    pub charge: f64,
    /// Basis set data of BSE format.
    pub basis_data: BseBasisElement,
}

/// Resolve basis sets for a list of atoms.
///
/// # Arguments
/// * `atoms` - List of parsed atom information
/// * `basis_spec` - Basis set specification
///
/// # Returns
/// HashMap mapping atom label to resolved BasisElement.
pub fn resolve_basis(
    atoms: &[AtomInfo],
    basis_spec: &BasisSpec,
    ecp_spec: &BasisSpec,
    ghost_ecp: bool,
) -> Result<BTreeMap<String, BasisElement>, CIntError> {
    let mut result = BTreeMap::new();

    for atom in atoms.iter() {
        // skip if label already processed
        if result.contains_key(&atom.label) {
            continue;
        }

        // resolve basis for this atom
        let mut basis_data = resolve_basis_for_atom(atom, basis_spec)?;

        // handle ecp basis
        if let Ok(ecp_data) = resolve_basis_for_atom(atom, ecp_spec) {
            basis_data.ecp_potentials = ecp_data.ecp_potentials;
            basis_data.ecp_electrons = ecp_data.ecp_electrons;
        };

        // ecp should not be assigned to ghost atoms, for usual cases
        if atom.is_ghost && !ghost_ecp {
            basis_data.ecp_electrons = None;
            basis_data.ecp_potentials = None;
        }

        result.insert(atom.label.clone(), BasisElement { symbol: atom.symbol.clone(), charge: atom.charge, basis_data });
    }

    Ok(result)
}

/// Resolve basis set for a single atom.
fn resolve_basis_for_atom(atom: &AtomInfo, spec: &BasisSpec) -> Result<BseBasisElement, CIntError> {
    match spec {
        BasisSpec::Uniform(input) => {
            // Ghost atoms get empty basis with uniform spec
            resolve_basis_input(input, &atom.symbol)
        },
        BasisSpec::Dict(map) => {
            // Try label match first (for ghost-O, X_H, etc.)
            if let Some(input) = map.get(&atom.label) {
                resolve_basis_input(input, &atom.symbol)
            } else {
                // Try element symbol match
                let key = atom.symbol.to_uppercase();
                if let Some(input) = map.get(&key) {
                    resolve_basis_input(input, &atom.symbol)
                } else if let Some(input) = map.get(&atom.symbol) {
                    resolve_basis_input(input, &atom.symbol)
                } else if let Some(input) = map.get("default") {
                    resolve_basis_input(input, &atom.symbol)
                } else {
                    cint_raise!(ParseError, "No basis specified for element '{}' in dict basis map", atom.symbol)
                }
            }
        },
        BasisSpec::List(list) => {
            // List is like Uniform for each item, but check for duplicate elements
            // For now, treat as uniform - the first matching input
            // Find matching input for this element
            for input in list {
                // Try to resolve and see if it matches this element
                if let Ok(elem) = resolve_basis_input(input, &atom.symbol) {
                    return Ok(elem);
                }
            }
            cint_raise!(ParseError, "No matching basis in list for element '{}'", atom.symbol)
        },
    }
}

/// Resolve basis from BasisInput.
fn resolve_basis_input(input: &BasisInput, element_symbol: &str) -> Result<BseBasisElement, CIntError> {
    match input {
        BasisInput::None => Ok(BseBasisElement::default()),
        BasisInput::String(str) => parse_basis_format(str, element_symbol),
        BasisInput::Element(elem) => Ok((**elem).clone()),
        BasisInput::Basis(basis) => {
            // Extract element from full basis
            let z = atom::symbol_to_charge(element_symbol)?;
            let elem_key = z.to_string();
            basis
                .elements
                .get(&elem_key)
                .cloned()
                .ok_or_else(|| cint_error!(ParseError, "Element {element_symbol} not found in provided basis"))
        },
    }
}

fn parse_basis_format(input: &str, element_symbol: &str) -> Result<BseBasisElement, CIntError> {
    let input = input.trim();
    // if input is nothing, return empty basis
    if input.is_empty() {
        return Ok(BseBasisElement::default());
    }
    // try fetch basis by name first
    if let Ok(elem) = fetch_basis_element(input, element_symbol) {
        return Ok(elem);
    }
    // try formats, only lines are larger than one
    if input.lines().count() > 1 {
        let formats = ["json", "nwchem", "gaussian94", "cp2k", "gamess_us", "turbomole", "molcas", "cfour"];
        for fmt in &formats {
            if let Ok(elem) = parse_basis_format_element(input, element_symbol, fmt) {
                return Ok(elem);
            }
        }
    }
    cint_raise!(ParseError, "Failed to parse basis token '{input}' for element '{element_symbol}' in any known format")
}

/// Fetch basis element data from BSE library.
fn fetch_basis_element(basis_name: &str, element_symbol: &str) -> Result<BseBasisElement, CIntError> {
    // Get atomic number from element symbol
    let z = atom::symbol_to_charge(element_symbol)?;

    // Fetch basis set from BSE
    let args = BseGetBasisArgsBuilder::default()
        .elements(z.to_string())
        .build()
        .map_err(|e| cint_error!(ParseError, "Failed to build BSE args: {e}"))?;

    let basis = get_basis_f(basis_name, args)
        .map_err(|e| cint_error!(ParseError, "Failed to fetch basis '{basis_name}' for element '{element_symbol}': {e}"))?;

    // Extract element data
    let elem_key = z.to_string();
    basis
        .elements
        .get(&elem_key)
        .cloned()
        .ok_or_else(|| cint_error!(ParseError, "Element {element_symbol} not found in basis '{basis_name}'"))
}

/// Parse a formatted basis string for a specific element.
fn parse_basis_format_element(basis_token: &str, element_symbol: &str, fmt: &str) -> Result<BseBasisElement, CIntError> {
    let z = atom::symbol_to_charge(element_symbol)?;

    // Try to parse as formatted basis string
    let basis = read_formatted_basis_str_f(basis_token, fmt)
        .map_err(|e| cint_error!(ParseError, "Failed to parse basis token '{basis_token}' as format '{fmt}': {e}"))?;

    // Extract element data
    let elem_key = z.to_string();
    basis.elements.get(&elem_key).cloned().ok_or_else(|| {
        cint_error!(ParseError, "Element {element_symbol} not found in parsed basis from token '{basis_token}' with format '{fmt}'")
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::parse::atom::{parse_atom_string, Unit};

    #[test]
    fn test_basis_spec_uniform() {
        let spec = BasisSpec::from("def2-TZVP");
        assert!(matches!(spec, BasisSpec::Uniform(BasisInput::String(_))));
    }

    #[test]
    fn test_basis_spec_dict() {
        let mut map = BTreeMap::new();
        map.insert("H".to_string(), BasisInput::String("sto-3g".to_string()));
        map.insert("O".to_string(), BasisInput::String("cc-pvdz".to_string()));
        let spec = BasisSpec::Dict(map);
        assert!(matches!(spec, BasisSpec::Dict(_)));
    }

    #[test]
    fn test_resolve_basis_uniform() {
        let atoms = parse_atom_string("H 0 0 0; O 0 0 1.2", Unit::Angstrom).unwrap();
        let spec = BasisSpec::from("sto-3g");
        let basis = resolve_basis(&atoms, &spec, &(None.into()), false).unwrap();

        assert_eq!(basis.len(), 2);
        // Keys are atom labels (H, O)
        assert!(basis.contains_key("H"));
        assert!(basis.contains_key("O"));
    }

    #[test]
    fn test_resolve_basis_ghost() {
        let atoms = parse_atom_string("H 0 0 0; GHOST-O 0 0 1.2", Unit::Angstrom).unwrap();
        let spec = BasisSpec::from("sto-3g");
        let basis = resolve_basis(&atoms, &spec, &(None.into()), false).unwrap();

        // Ghost atom has empty basis (keyed by label "GHOST-O")
        assert!(basis.get("GHOST-O").unwrap().basis_data.electron_shells.is_some());
    }

    #[test]
    fn test_fetch_basis_element() {
        let elem = fetch_basis_element("sto-3g", "H").unwrap();
        assert!(elem.electron_shells.is_some());
        let shells = elem.electron_shells.unwrap();
        assert!(!shells.is_empty());
    }
}
