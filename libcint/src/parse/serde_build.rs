//! Serde-based molecule building from JSON and TOML formats.
//!
//! This module provides `CIntMol::from_json` and `CIntMol::from_toml` methods
//! for building molecule instances from structured input files.
//!
//! # TOML Format Example
//!
//! ```toml
//! atom = """
//! O 0.000000 0.000000 0.000000
//! H 0.000000 0.000000 0.957200
//! H 0.000000 0.926627 -0.239987
//! """
//! basis = "STO-3G"
//! ecp = "LANL2DZ"
//! unit = "Angstrom"
//! cart = false
//! ghost_ecp = false
//! allow_empty_basis = false
//! ```
//!
//! # Custom Basis/ECP Format
//!
//! For custom basis or ECP sets with per-element specifications:
//!
//! ```toml
//! atom = """
//! O 0.000000 0.000000 0.000000
//! H 0.000000 0.000000 0.957200
//! """
//! basis = "custom"
//!
//! [basis-custom]
//! O = "STO-3G"
//! default = "6-31G"
//! ```

use crate::parse::atom::Unit;
use crate::parse::basis::{BasisInput, BasisSpec};
use crate::parse::mole::{CIntMol, CIntMolInput};
use crate::prelude::*;
use indexmap::IndexMap;
use serde::{Deserialize, Deserializer, Serialize, Serializer};
use std::collections::HashMap;

/// Serde-friendly input structure for molecule building.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SerdeMolInput {
    /// Atom coordinates string (required).
    pub atom: String,
    /// Basis set specification.
    #[serde(default)]
    pub basis: SerdeBasisSpec,
    /// ECP specification.
    #[serde(default)]
    pub ecp: SerdeBasisSpec,
    /// Unit of coordinates (optional, default: Angstrom).
    #[serde(default)]
    pub unit: Unit,
    /// Use cartesian basis (optional, default: false).
    #[serde(default)]
    pub cart: bool,
    /// Whether to assign ECP to ghost atoms.
    #[serde(default)]
    pub ghost_ecp: bool,
    /// Allow empty basis.
    #[serde(default)]
    pub allow_empty_basis: bool,
}

/// Serde-friendly wrapper for BasisSpec.
///
/// Can deserialize from:
/// - A simple string: `"STO-3G"` -> Uniform basis
/// - A dictionary: `{"O": "STO-3G", "default": "6-31G"}` -> Dict basis
/// - The value "custom" (requires separate `[basis-custom]` table)
#[derive(Debug, Clone, PartialEq, Default)]
pub enum SerdeBasisSpec {
    #[default]
    None,
    Uniform(String),
    Dict(IndexMap<String, String>),
    Custom,
}

impl Serialize for SerdeBasisSpec {
    fn serialize<S>(&self, serializer: S) -> Result<S::Ok, S::Error>
    where
        S: Serializer,
    {
        match self {
            SerdeBasisSpec::None => serializer.serialize_str(""),
            SerdeBasisSpec::Uniform(s) => serializer.serialize_str(s),
            SerdeBasisSpec::Dict(map) => {
                use serde::ser::SerializeMap;
                let mut map_ser = serializer.serialize_map(Some(map.len()))?;
                for (k, v) in map.iter() {
                    map_ser.serialize_entry(k, v)?;
                }
                map_ser.end()
            },
            SerdeBasisSpec::Custom => serializer.serialize_str("custom"),
        }
    }
}

impl<'de> Deserialize<'de> for SerdeBasisSpec {
    fn deserialize<D>(deserializer: D) -> Result<Self, D::Error>
    where
        D: Deserializer<'de>,
    {
        // Try to deserialize as a string first, if fails try as a map
        // This is a common pattern for "string or map" deserialization

        let value = serde_json::Value::deserialize(deserializer)?;
        match value {
            serde_json::Value::String(s) => {
                if s == "custom" {
                    Ok(SerdeBasisSpec::Custom)
                } else if s.is_empty() {
                    Ok(SerdeBasisSpec::None)
                } else {
                    Ok(SerdeBasisSpec::Uniform(s))
                }
            },
            serde_json::Value::Object(map) => {
                let result: IndexMap<String, String> = map
                    .into_iter()
                    .map(|(k, v)| {
                        let v_str = match v {
                            serde_json::Value::String(s) => s,
                            _ => v.to_string(),
                        };
                        (k.to_ascii_uppercase(), v_str)
                    })
                    .collect();
                Ok(SerdeBasisSpec::Dict(result))
            },
            _ => Err(serde::de::Error::custom("expected string or map")),
        }
    }
}

impl From<SerdeBasisSpec> for BasisSpec {
    fn from(spec: SerdeBasisSpec) -> Self {
        match spec {
            SerdeBasisSpec::None => BasisSpec::Uniform(BasisInput::None),
            SerdeBasisSpec::Uniform(s) => BasisSpec::Uniform(BasisInput::String(s)),
            SerdeBasisSpec::Dict(map) => {
                let dict: IndexMap<String, BasisInput> = map.into_iter().map(|(k, v)| (k, BasisInput::String(v))).collect();
                BasisSpec::Dict(dict)
            },
            SerdeBasisSpec::Custom => BasisSpec::Uniform(BasisInput::None),
        }
    }
}

/// Raw TOML input that captures custom basis/ECP tables.
#[derive(Debug, Clone, Default, Deserialize)]
struct RawTomlInput {
    atom: String,
    #[serde(default)]
    basis: SerdeBasisSpec,
    #[serde(default)]
    ecp: SerdeBasisSpec,
    #[serde(default)]
    unit: Unit,
    #[serde(default)]
    cart: bool,
    #[serde(default)]
    ghost_ecp: bool,
    #[serde(default)]
    allow_empty_basis: bool,
    #[serde(default, rename = "basis-custom")]
    basis_custom: HashMap<String, String>,
    #[serde(default, rename = "ecp-custom")]
    ecp_custom: HashMap<String, String>,
}

impl SerdeMolInput {
    /// Convert to CIntMolInput for molecule building.
    pub fn to_mol_input(self) -> CIntMolInput {
        CIntMolInput {
            atom: self.atom,
            basis: self.basis.into(),
            ecp: self.ecp.into(),
            unit: self.unit,
            cart: self.cart,
            ghost_ecp: self.ghost_ecp,
            allow_empty_basis: self.allow_empty_basis,
        }
    }

    /// Build molecule from this input.
    pub fn build(self) -> Result<CIntMol, CIntError> {
        self.to_mol_input().create_mol_f()
    }
}

/// Parse TOML input and build molecule.
pub fn parse_toml_input(toml_str: &str) -> Result<CIntMol, CIntError> {
    let raw: RawTomlInput = toml::from_str(toml_str)
        .map_err(|e| cint_error!(ParseError, "Failed to parse TOML: {e}"))?;

    // Resolve custom basis if needed
    let basis = if matches!(raw.basis, SerdeBasisSpec::Custom) {
        if raw.basis_custom.is_empty() {
            return cint_raise!(ParseError, "basis = \"custom\" specified but no [basis-custom] table found");
        }
        let map: IndexMap<String, String> = raw
            .basis_custom
            .into_iter()
            .map(|(k, v)| (k.to_ascii_uppercase(), v))
            .collect();
        SerdeBasisSpec::Dict(map)
    } else {
        raw.basis
    };

    // Resolve custom ECP if needed
    let ecp = if matches!(raw.ecp, SerdeBasisSpec::Custom) {
        if raw.ecp_custom.is_empty() {
            return cint_raise!(ParseError, "ecp = \"custom\" specified but no [ecp-custom] table found");
        }
        let map: IndexMap<String, String> = raw
            .ecp_custom
            .into_iter()
            .map(|(k, v)| (k.to_ascii_uppercase(), v))
            .collect();
        SerdeBasisSpec::Dict(map)
    } else {
        raw.ecp
    };

    let input = SerdeMolInput {
        atom: raw.atom,
        basis,
        ecp,
        unit: raw.unit,
        cart: raw.cart,
        ghost_ecp: raw.ghost_ecp,
        allow_empty_basis: raw.allow_empty_basis,
    };

    input.build()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_serde_unit() {
        let json = "\"angstrom\"";
        let unit: Unit = serde_json::from_str(json).unwrap();
        assert_eq!(unit, Unit::Angstrom);

        let json = "\"bohr\"";
        let unit: Unit = serde_json::from_str(json).unwrap();
        assert_eq!(unit, Unit::Bohr);

        // Case insensitive
        let json = "\"BOHR\"";
        let unit: Unit = serde_json::from_str(json).unwrap();
        assert_eq!(unit, Unit::Bohr);
    }

    #[test]
    fn test_serde_basis_spec_uniform() {
        let json = "\"STO-3G\"";
        let spec: SerdeBasisSpec = serde_json::from_str(json).unwrap();
        assert_eq!(spec, SerdeBasisSpec::Uniform("STO-3G".to_string()));
    }

    #[test]
    fn test_serde_basis_spec_custom() {
        let json = "\"custom\"";
        let spec: SerdeBasisSpec = serde_json::from_str(json).unwrap();
        assert_eq!(spec, SerdeBasisSpec::Custom);
    }

    #[test]
    fn test_serde_basis_spec_dict() {
        let json = "{\"O\": \"STO-3G\", \"default\": \"6-31G\"}";
        let spec: SerdeBasisSpec = serde_json::from_str(json).unwrap();
        let map = match spec {
            SerdeBasisSpec::Dict(map) => map,
            _ => panic!("Expected Dict"),
        };
        assert_eq!(map.get("O"), Some(&"STO-3G".to_string()));
        assert_eq!(map.get("DEFAULT"), Some(&"6-31G".to_string()));
    }

    #[test]
    fn test_serde_mol_input_json() {
        let json = r#"{
            "atom": "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987",
            "basis": "STO-3G",
            "unit": "angstrom",
            "cart": false
        }"#;

        let input: SerdeMolInput = serde_json::from_str(json).unwrap();
        assert_eq!(input.atom, "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987");
        assert_eq!(input.basis, SerdeBasisSpec::Uniform("STO-3G".to_string()));
        assert_eq!(input.unit, Unit::Angstrom);
        assert!(!input.cart);
    }

    #[test]
    fn test_serde_mol_input_toml() {
        let toml = r#"
            atom = "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987"
            basis = "STO-3G"
            unit = "angstrom"
            cart = false
        "#;

        let input: SerdeMolInput = toml::from_str(toml).unwrap();
        assert_eq!(input.atom, "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987");
        assert_eq!(input.basis, SerdeBasisSpec::Uniform("STO-3G".to_string()));
        assert_eq!(input.unit, Unit::Angstrom);
        assert!(!input.cart);
    }

    #[test]
    fn test_serde_mol_input_toml_custom_basis() {
        let toml = r#"
            atom = "O 0 0 0; H 0 0 0.9572"
            basis = "custom"

            [basis-custom]
            O = "STO-3G"
            default = "6-31G"
        "#;

        let raw: RawTomlInput = toml::from_str(toml).unwrap();
        assert!(matches!(raw.basis, SerdeBasisSpec::Custom));
        assert_eq!(raw.basis_custom.get("O"), Some(&"STO-3G".to_string()));
        assert_eq!(raw.basis_custom.get("default"), Some(&"6-31G".to_string()));
    }
}