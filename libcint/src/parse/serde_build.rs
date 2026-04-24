//! Serde-based molecule building from JSON and TOML formats.
//!
//! This module provides `CIntMol::from_toml` method for building molecule
//! instances from TOML files with custom basis/ECP table support.

use crate::parse::basis::{BasisInput, BasisSpec};
use crate::parse::mole::{CIntMol, CIntMolInput};
use crate::prelude::*;

/// TOML input wrapper that captures custom basis/ECP tables.
///
/// Handles the special case where `basis = "custom"` requires a
/// `[basis-custom]` table. Uses `#[serde(flatten)]` to capture remaining
/// CIntMolInput fields.
#[derive(Debug, Clone, Default, Deserialize)]
struct TomlInputWithCustomTables {
    #[serde(flatten)]
    input: CIntMolInput,
    #[serde(default, rename = "basis-custom")]
    basis_custom: HashMap<String, String>,
    #[serde(default, rename = "ecp-custom")]
    ecp_custom: HashMap<String, String>,
}

/// Parse TOML input and build molecule.
///
/// Handles the special "custom" marker for basis/ECP that requires separate
/// tables.
pub fn parse_toml_input(toml_str: &str) -> Result<CIntMol, CIntError> {
    // First parse as raw toml::Value to check for "custom" markers
    let raw_toml: toml::Value = toml::from_str(toml_str).map_err(|e| cint_error!(ParseError, "Failed to parse TOML: {e}"))?;

    // Check if basis/ecp are "custom" strings
    let basis_is_custom = check_is_custom(&raw_toml, "basis");
    let ecp_is_custom = check_is_custom(&raw_toml, "ecp");

    // Parse into TomlInputWithCustomTables
    let raw: TomlInputWithCustomTables = toml::from_str(toml_str).map_err(|e| cint_error!(ParseError, "Failed to parse TOML: {e}"))?;

    // Resolve custom basis if needed
    let basis = if basis_is_custom {
        if raw.basis_custom.is_empty() {
            return cint_raise!(ParseError, "basis = \"custom\" specified but no [basis-custom] table found");
        }
        let dict: IndexMap<String, BasisInput> =
            raw.basis_custom.into_iter().map(|(k, v)| (k.to_ascii_uppercase(), BasisInput::String(v))).collect();
        BasisSpec::Dict(dict)
    } else {
        raw.input.basis
    };

    // Resolve custom ECP if needed
    let ecp = if ecp_is_custom {
        if raw.ecp_custom.is_empty() {
            return cint_raise!(ParseError, "ecp = \"custom\" specified but no [ecp-custom] table found");
        }
        let dict: IndexMap<String, BasisInput> =
            raw.ecp_custom.into_iter().map(|(k, v)| (k.to_ascii_uppercase(), BasisInput::String(v))).collect();
        BasisSpec::Dict(dict)
    } else {
        raw.input.ecp
    };

    // Build CIntMolInput with resolved basis/ecp
    let input = CIntMolInput {
        atom: raw.input.atom,
        basis,
        ecp,
        unit: raw.input.unit,
        cart: raw.input.cart,
        ghost_ecp: raw.input.ghost_ecp,
        allow_empty_basis: raw.input.allow_empty_basis,
    };

    input.create_mol_f()
}

/// Check if a field in TOML is the "custom" string marker.
fn check_is_custom(toml: &toml::Value, field: &str) -> bool {
    if let toml::Value::Table(table) = toml {
        if let Some(toml::Value::String(s)) = table.get(field) {
            return s == "custom";
        }
    }
    false
}

#[cfg(test)]
mod tests {
    use crate::parse::atom::Unit;
    use crate::parse::basis::{BasisInput, BasisSpec};
    use crate::parse::mole::{CIntMol, CIntMolInput};

    #[test]
    fn test_cint_mol_input_json() {
        let json = r#"{
            "atom": "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987",
            "basis": "STO-3G",
            "unit": "angstrom",
            "cart": false
        }"#;

        let input: CIntMolInput = serde_json::from_str(json).unwrap();
        assert_eq!(input.atom, "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987");
        assert!(matches!(input.basis, BasisSpec::Uniform(BasisInput::String(_))));
        assert_eq!(input.unit, Unit::Angstrom);
        assert!(!input.cart);
    }

    #[test]
    fn test_cint_mol_input_toml() {
        let toml = r#"
            atom = "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987"
            basis = "STO-3G"
            unit = "angstrom"
            cart = false
        "#;

        let input: CIntMolInput = toml::from_str(toml).unwrap();
        assert_eq!(input.atom, "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987");
        assert!(matches!(input.basis, BasisSpec::Uniform(BasisInput::String(_))));
        assert_eq!(input.unit, Unit::Angstrom);
        assert!(!input.cart);
    }

    #[test]
    fn test_toml_custom_basis() {
        let toml = r#"
            atom = "O 0 0 0; H 0 0 0.9572"
            basis = "custom"

            [basis-custom]
            O = "STO-3G"
            default = "6-31G"
        "#;

        let mol = CIntMol::from_toml(toml);
        assert!(mol.cint.nao() > 0);
    }

    #[test]
    fn test_toml_list_basis() {
        let toml = r#"
            atom = "O 0 0 0; H 0 0 0.9572"
            basis = ["STO-3G", "6-31G"]
        "#;

        let input: CIntMolInput = toml::from_str(toml).unwrap();
        assert!(matches!(input.basis, BasisSpec::List(_)));
    }

    #[test]
    fn test_basis_spec_serialize() {
        let spec = BasisSpec::Uniform(BasisInput::String("STO-3G".to_string()));
        let json = serde_json::to_string(&spec).unwrap();
        assert_eq!(json, "\"STO-3G\"");
    }

    #[test]
    fn test_cint_mol_input_serialize() {
        let input = CIntMolInput {
            atom: "O 0 0 0".to_string(),
            basis: BasisSpec::Uniform(BasisInput::String("STO-3G".to_string())),
            ecp: BasisSpec::Uniform(BasisInput::None),
            unit: Unit::Angstrom,
            cart: false,
            ghost_ecp: false,
            allow_empty_basis: false,
        };
        let json = serde_json::to_string(&input).unwrap();
        assert!(json.contains("\"atom\""));
        assert!(json.contains("\"basis\""));
    }
}
