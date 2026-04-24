//! Basis set resolution and handling.

use super::atom;
use crate::parse::atom::AtomInfo;
use crate::prelude::*;
use bse::manip::uncontract_spdf_in_element;

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

/// Specification for basis sets, which can be uniform (same for all atoms), a
/// dict mapping atom labels to basis, or a list of basis inputs to try in
/// order.
///
/// Please note for dict basis, the keys must be uppercase (`AU` instead of
/// `Au`) to be matched.
#[derive(Debug, Clone, PartialEq)]
pub enum BasisSpec {
    Uniform(BasisInput),
    Dict(IndexMap<String, BasisInput>),
    List(Vec<BasisInput>),
}

impl Default for BasisSpec {
    fn default() -> Self {
        BasisSpec::Uniform("".into())
    }
}

#[allow(clippy::useless_conversion)]
#[duplicate_item(TY; [String]; [&str]; [Option<String>]; [BseBasisElement]; [BseBasis];)]
impl From<TY> for BasisSpec {
    fn from(s: TY) -> Self {
        BasisSpec::Uniform(s.into())
    }
}

#[duplicate_item(TY;
    [HashMap<String, T>]; [HashMap<&str, T>];
    [BTreeMap<String, T>]; [BTreeMap<&str, T>];
    [IndexMap<String, T>]; [IndexMap<&str, T>];
    [Vec<(String, T)>]; [Vec<(&str, T)>];
)]
impl<T> From<TY> for BasisSpec
where
    T: Into<BasisInput>,
{
    fn from(map: TY) -> Self {
        let dict = map.into_iter().map(|(k, v)| (k.to_ascii_uppercase(), v.into())).collect();
        BasisSpec::Dict(dict)
    }
}

#[duplicate_item(TY; [String]; [&str])]
impl<T, const N: usize> From<[(TY, T); N]> for BasisSpec
where
    T: Into<BasisInput>,
{
    fn from(arr: [(TY, T); N]) -> Self {
        let dict = arr.into_iter().map(|(k, v)| (k.to_ascii_uppercase(), v.into())).collect();
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

/// Resolve basis sets for a list of atoms.
///
/// # Arguments
/// - `atoms`: List of parsed atom information
/// - `basis_spec`: Basis set specification
/// - `ecp_spec`: ECP specification (same format as basis_spec)
/// - `ghost_ecp`: Whether to assign ECP to ghost atoms (default: false)
///
/// # Returns
/// - A mapping from atom labels to resolved basis elements
pub fn resolve_basis(
    atoms: &[AtomInfo],
    basis_spec: &BasisSpec,
    ecp_spec: &BasisSpec,
    ghost_ecp: bool,
) -> Result<(IndexMap<String, BseBasisElement>, Vec<String>), CIntError> {
    // step 1: generate the dictionary and atom label list first
    let mut result = BTreeMap::new();
    let mut name_list: Vec<String> = Vec::new(); // list of parsed
    let mut name_map: BTreeMap<&str, String> = BTreeMap::new(); // atom.label -> parsed

    for atom in atoms.iter() {
        // skip if label already processed
        if let Some(parsed) = name_map.get(atom.label.as_str()) {
            name_list.push(parsed.clone());
            continue;
        }

        // resolve basis for this atom
        let (mut basis_data, mut parsed_name) = resolve_basis_for_atom(atom, basis_spec)?;

        // handle ecp basis
        if let Ok((ecp_data, parsed_name_ecp)) = resolve_basis_for_atom(atom, ecp_spec) {
            basis_data.ecp_potentials = ecp_data.ecp_potentials;
            basis_data.ecp_electrons = ecp_data.ecp_electrons;
            // if ecp basis have larger precedence, use the ecp basis name
            if parsed_name_ecp.len() > parsed_name.len() || parsed_name == "DEFAULT" {
                parsed_name = parsed_name_ecp;
            }
        };

        // ecp should not be assigned to ghost atoms, for usual cases
        // this is very special case, use the label name to represent the basis
        if atom.is_ghost && !ghost_ecp {
            basis_data.ecp_electrons = None;
            basis_data.ecp_potentials = None;
            parsed_name = atom.label.clone();
        }

        // pyscf usually use the uncontract version
        uncontract_spdf_in_element(&mut basis_data, 0);
        // follows bse's convention of nwchem, sort shells here
        if let Some(ref mut electron_shells) = basis_data.electron_shells {
            bse::sort::sort_shells(electron_shells);
        }
        if let Some(ref mut ecp_potentials) = basis_data.ecp_potentials {
            bse::sort::sort_potentials(ecp_potentials);
        }

        // check if the result have already have this basis (by parsed name), if so,
        // check if the basis data is the same, otherwise raise.
        if let Some(existing) = result.get(&parsed_name) {
            if existing != &basis_data {
                cint_raise!(ParseError, "The basis parsing seems to give two different results with the same entry {parsed_name}.")?
            }
        } else {
            result.insert(parsed_name.clone(), basis_data.clone());
        }
        name_list.push(parsed_name.clone());
        name_map.insert(atom.label.as_str(), parsed_name);
    }

    // step 2. reorder btree with input order
    // order to be expected: BasisSpec::Dict's keys, then atom labels
    let mut order_guide = vec![];
    if let BasisSpec::Dict(dict) = basis_spec {
        for key in dict.keys() {
            order_guide.push(key.clone());
        }
    }
    order_guide.extend(name_list.iter().cloned());

    // create an ordered result based on the order guide
    let mut ordered_result = IndexMap::new();
    for key in order_guide {
        if let Some(value) = result.remove(&key) {
            ordered_result.insert(key.clone(), value);
        }
    }

    // check the result is moved correctly
    if !result.is_empty() {
        cint_raise!(ParseError, "Some parsed basis entries are not included in the final result: {:?}, probably bug.", result.keys())?
    }

    Ok((ordered_result, name_list))
}

/// Resolve basis set for a single atom.
///
/// The returned `String` is the parsed string for representing this basis (may
/// be label, identifier or symbol of atom_info).
fn resolve_basis_for_atom(atom: &AtomInfo, spec: &BasisSpec) -> Result<(BseBasisElement, String), CIntError> {
    match spec {
        BasisSpec::Uniform(input) => {
            // Ghost atoms get empty basis with uniform spec
            Ok((resolve_basis_input(input, &atom.symbol)?, atom.symbol.clone()))
        },
        BasisSpec::List(list) => {
            // treat as uniform, find first match
            for input in list {
                // Try to resolve and see if it matches this element
                if let Ok(parsed) = resolve_basis_input(input, &atom.symbol) {
                    return Ok((parsed, atom.symbol.clone()));
                }
            }
            cint_raise!(ParseError, "No matching basis in list for element '{}'", atom.symbol)
        },
        BasisSpec::Dict(map) => {
            // match by label, then identifier, then symbol, then default, then error
            if let Some(input) = map.get(&atom.label) {
                Ok((resolve_basis_input(input, &atom.symbol)?, atom.label.clone()))
            } else if let Some(input) = map.get(&atom.identifier) {
                Ok((resolve_basis_input(input, &atom.symbol)?, atom.identifier.clone()))
            } else if let Some(input) = map.get(&atom.symbol) {
                Ok((resolve_basis_input(input, &atom.symbol)?, atom.symbol.clone()))
            } else if let Some(input) = map.get("DEFAULT") {
                Ok((resolve_basis_input(input, &atom.symbol)?, atom.symbol.clone()))
            } else {
                cint_raise!(ParseError, "No matching basis in dict for element '{}'", atom.symbol)
            }
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

/// Parse basis from a string input, which can be either a basis name or a
/// formatted basis string.
///
/// This function uses bse.
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
        let formats = ["json", "nwchem", "gaussian94", "cp2k", "gamess_us", "turbomole", "molcas", "cfour", "molpro", "crystal"];
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
        let mut map = IndexMap::new();
        map.insert("H".to_string(), BasisInput::String("sto-3g".to_string()));
        map.insert("O".to_string(), BasisInput::String("cc-pvdz".to_string()));
        let spec = BasisSpec::Dict(map);
        assert!(matches!(spec, BasisSpec::Dict(_)));
    }

    #[test]
    fn test_resolve_basis_uniform() {
        let atoms = parse_atom_string("H 0 0 0; O 0 0 1.2", Unit::Angstrom).unwrap();
        let spec = BasisSpec::from("sto-3g");
        let (basis, _) = resolve_basis(&atoms, &spec, &(None.into()), false).unwrap();

        assert_eq!(basis.len(), 2);
        // Keys are atom labels (H, O)
        assert!(basis.contains_key("H"));
        assert!(basis.contains_key("O"));
    }

    #[test]
    fn test_resolve_basis_ghost() {
        let atoms = parse_atom_string("H 0 0 0; GHOST-O 0 0 1.2", Unit::Angstrom).unwrap();
        let spec = BasisSpec::from("sto-3g");
        let (basis, _) = resolve_basis(&atoms, &spec, &(None.into()), false).unwrap();

        // Ghost atom has empty basis (keyed by label "GHOST-O")
        assert!(basis.get("GHOST-O").unwrap().electron_shells.is_some());
    }

    #[test]
    fn test_fetch_basis_element() {
        let elem = fetch_basis_element("sto-3g", "H").unwrap();
        assert!(elem.electron_shells.is_some());
        let shells = elem.electron_shells.unwrap();
        assert!(!shells.is_empty());
    }
}
