//! Atom coordinate parsing.
//!
//! This module provides functions to parse atom coordinates from various
//! formats:
//! - String format: `"H 0 0 0; O 0 0 1.2"`
//! - List format: `["H 0 0 0", "O 0 0 1.2"]`
//! - Z-matrix format: `"O; H 1 0.94; H 1 0.94 2 104.5"`
//!
//! Supports ghost atoms (`GHOST-X`, `X-` prefix) and labeled atoms (`H1`,
//! `H2`).

use crate::prelude::*;
use bse::lut::{element_Z_from_sym, element_sym_from_Z_with_normalize};

/// Bohr to Angstrom conversion constant.
///
/// <https://physics.nist.gov/cgi-bin/cuu/Value?bohrrada0>
pub const BOHR_TO_ANG: f64 = 0.529177210544;

/// Angstrom to Bohr conversion constant.
pub const ANG_TO_BOHR: f64 = 1.0 / BOHR_TO_ANG;

/// Unit of length for coordinates.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum Unit {
    #[default]
    Angstrom,
    Bohr,
}

/// Information about a parsed atom.
#[derive(Debug, Clone, PartialEq)]
pub struct AtomInfo {
    /// Raw symbol/label from input (e.g., "H", "H1", "GHOST-O", "O@2").
    pub label: String,
    /// Element symbol without label suffix (e.g., "H", "O").
    pub symbol: String,
    /// Atomic number (Z). 0 for ghost atoms.
    pub charge: f64,
    /// Whether this is a ghost atom.
    pub is_ghost: bool,
    /// Cartesian coordinates in Bohr.
    pub coords: [f64; 3],
}

impl AtomInfo {
    /// Create a new AtomInfo from symbol and coordinates.
    pub fn new(symbol: &str, coords: [f64; 3], unit: Unit) -> Self {
        let (parsed_symbol, is_ghost) = parse_atom_symbol(symbol);
        let charge = if is_ghost {
            0.0
        } else {
            element_charge(&parsed_symbol).unwrap_or_else(|| panic!("Unknown element: {}", parsed_symbol)) as f64
        };
        let coords_bohr =
            if unit == Unit::Angstrom { [coords[0] * ANG_TO_BOHR, coords[1] * ANG_TO_BOHR, coords[2] * ANG_TO_BOHR] } else { coords };
        AtomInfo { label: symbol.to_string(), symbol: parsed_symbol, charge, is_ghost, coords: coords_bohr }
    }
}

/// Parse atom symbol, extracting element symbol and ghost status.
///
/// Handles formats:
/// - `H`, `O`, `C` → ("H", false)
/// - `H1`, `H2`, `O@1` → ("H", false)
/// - `GHOST-H`, `GHOST_O` → ("H", true)
/// - `X-H`, `X_O` → ("H", true)
/// - `H@`, `O@` (trailing @) → ("H", true)
fn parse_atom_symbol(s: &str) -> (String, bool) {
    let s_trimmed = s.trim();
    let s_upper = s_trimmed.to_uppercase();

    // Check for ghost atom markers: GHOST-H, GHOST_O, ghost.H, etc.
    if let Some(rest_upper) = s_upper.strip_prefix("GHOST") {
        let rest = rest_upper.trim_start_matches(['-', '_', '.']);
        let (symbol, _) = parse_atom_symbol(rest);
        return (symbol, true);
    }

    // Check for X-H, X_O ghost markers
    if s_upper.strip_prefix("X-").is_some() || s_upper.strip_prefix("X_").is_some() {
        let rest = s_trimmed[2..].trim();
        let (symbol, _) = parse_atom_symbol(rest);
        return (symbol, true);
    }

    // Check for trailing @ (ghost marker): H@, O@
    if let Some(rest) = s_trimmed.strip_suffix('@') {
        // Avoid double @@ (not ghost marker)
        if !rest.ends_with('@') {
            let (symbol, _) = parse_atom_symbol(rest);
            return (symbol, true);
        }
    }

    // Extract pure element symbol: strip digits and @ suffixes
    // H1 → H, H2 → H, O@1 → O, Au1 → Au
    let symbol = s_trimmed
        .trim_end_matches(|c: char| c.is_ascii_digit())
        .trim_end_matches('@')
        .trim_end_matches(|c: char| c.is_ascii_digit())
        .to_string();

    (symbol, false)
}

/// Get atomic number from element symbol using bse lookup.
pub fn element_charge(symbol: &str) -> Option<i32> {
    element_Z_from_sym(symbol)
}

/// Get element symbol from atomic number using bse lookup.
fn charge_to_element(z: i32) -> Option<String> {
    element_sym_from_Z_with_normalize(z)
}

/// Parse atom coordinates from a string.
///
/// Formats supported:
/// - `"H 0 0 0; O 0 0 1.2"` (semicolon separated)
/// - `"H 0 0 0\nO 0 0 1.2"` (newline separated)
/// - `"H,0,0,0; O,0,0,1.2"` (comma separated)
///
/// # Arguments
/// * `s` - String containing atom coordinates
/// * `unit` - Unit of input coordinates (Angstrom or Bohr)
///
/// # Returns
/// Vector of AtomInfo structs.
pub fn parse_atom_string(s: &str, unit: Unit) -> Result<Vec<AtomInfo>, CIntError> {
    // Replace semicolons with newlines for unified processing
    let s = s.replace(';', "\n").replace(',', " ");

    let atoms: Vec<AtomInfo> = s
        .lines()
        .filter(|line| {
            let trimmed = line.trim();
            !trimmed.is_empty() && !trimmed.starts_with('#')
        })
        .map(|line| parse_atom_line(line.trim(), unit))
        .collect::<Result<Vec<_>, _>>()?;

    Ok(atoms)
}

/// Parse a single atom line.
///
/// Format: `<symbol> <x> <y> <z>` or `<symbol> <x> <y> <z> <unit>`
///
/// Supports numeric atomic numbers: `1 0 0 0` → H atom.
fn parse_atom_line(line: &str, default_unit: Unit) -> Result<AtomInfo, CIntError> {
    let parts: Vec<&str> = line.split_whitespace().collect();

    if parts.len() < 4 {
        return cint_raise!(ParseError, "Invalid atom line: '{}'. Expected at least 4 parts.", line);
    }

    // Parse symbol (could be numeric charge)
    let symbol = parts[0];
    let (parsed_symbol, is_ghost) = if let Ok(z) = symbol.parse::<i32>() {
        // Numeric charge: convert to symbol
        let sym = charge_to_symbol(z)?;
        (sym, false)
    } else {
        parse_atom_symbol(symbol)
    };

    // Parse coordinates with explicit error conversion
    let x: f64 = parts[1]
        .parse()
        .map_err(|e: std::num::ParseFloatError| cint_error!(ParseError, "Failed to parse x coordinate: '{}': {}", parts[1], e))?;
    let y: f64 = parts[2]
        .parse()
        .map_err(|e: std::num::ParseFloatError| cint_error!(ParseError, "Failed to parse y coordinate: '{}': {}", parts[2], e))?;
    let z: f64 = parts[3]
        .parse()
        .map_err(|e: std::num::ParseFloatError| cint_error!(ParseError, "Failed to parse z coordinate: '{}': {}", parts[3], e))?;

    // Check for unit override in 5th position
    let unit = if parts.len() > 4 {
        match parts[4].to_uppercase().as_str() {
            "ANG" | "ANGSTROM" | "A" => Unit::Angstrom,
            "BOHR" | "AU" | "B" => Unit::Bohr,
            _ => default_unit,
        }
    } else {
        default_unit
    };

    // Convert to Bohr if necessary
    let coords_bohr = if unit == Unit::Angstrom { [x * ANG_TO_BOHR, y * ANG_TO_BOHR, z * ANG_TO_BOHR] } else { [x, y, z] };

    let charge =
        if is_ghost { 0.0 } else { element_charge(&parsed_symbol).unwrap_or_else(|| panic!("Unknown element: {}", parsed_symbol)) as f64 };

    Ok(AtomInfo { label: symbol.to_string(), symbol: parsed_symbol, charge, is_ghost, coords: coords_bohr })
}

/// Convert atomic number to element symbol.
fn charge_to_symbol(z: i32) -> Result<String, CIntError> {
    charge_to_element(z).ok_or_else(|| cint_error!(ParseError, "Invalid atomic number: {}", z))
}

/// Parse z-matrix format coordinates.
///
/// Format:
/// ```text
/// Atom1
/// Atom2 Atom1 bond_length
/// Atom3 Atom1 bond_length Atom2 angle_degrees
/// Atom4 Atom1 bond_length Atom2 angle_degrees Atom3 dihedral_degrees
/// ```
///
/// The first atom is placed at origin. The second atom is placed along z-axis.
/// The third atom is in the xz-plane. The fourth and subsequent atoms use
/// full z-matrix specification.
pub fn parse_zmatrix(s: &str, unit: Unit) -> Result<Vec<AtomInfo>, CIntError> {
    // Replace semicolons with newlines for unified processing
    let s = s.replace(';', "\n");

    let lines: Vec<&str> = s
        .lines()
        .filter(|line| {
            let trimmed = line.trim();
            !trimmed.is_empty() && !trimmed.starts_with('#')
        })
        .collect();

    if lines.is_empty() {
        return Ok(vec![]);
    }

    let mut atoms: Vec<AtomInfo> = vec![];
    let mut atom_order: Vec<String> = vec![];

    for (idx, line) in lines.iter().enumerate() {
        let parts: Vec<&str> = line.split_whitespace().collect();

        if parts.is_empty() {
            continue;
        }

        let symbol = parts[0];
        let (parsed_symbol, is_ghost) = parse_atom_symbol(symbol);

        let coords = match idx {
            0 => {
                // First atom at origin
                [0.0, 0.0, 0.0]
            },
            1 => {
                // Second atom along x-axis (PySCF convention)
                if parts.len() < 3 {
                    return cint_raise!(ParseError, "Invalid z-matrix line: '{}'. Need bond length.", line);
                }
                let bond: f64 = parse_float(parts[2], "bond length")?;
                let bond_bohr = if unit == Unit::Angstrom { bond * ANG_TO_BOHR } else { bond };
                [bond_bohr, 0.0, 0.0]
            },
            2 => {
                // Third atom in xz-plane (PySCF convention)
                // Bond distance is from atom 1 (r1 = O at origin), angle is from direction of
                // atom 2 (a1 = H1) r1 is the atom at origin, a1 is the
                // reference direction atom
                if parts.len() < 5 {
                    return cint_raise!(ParseError, "Invalid z-matrix line: '{}'. Need bond and angle.", line);
                }
                let r1_name = parts[1];
                let r1_idx = find_atom_index(&atom_order, r1_name)?;
                let r1_pos = atoms[r1_idx].coords;

                let a1_name = parts[3];
                let _a1_idx = find_atom_index(&atom_order, a1_name)?;
                let _a1_pos = atoms[_a1_idx].coords;

                let bond: f64 = parse_float(parts[2], "bond length")?;
                let bond_bohr = if unit == Unit::Angstrom { bond * ANG_TO_BOHR } else { bond };

                let angle_deg: f64 = parse_float(parts[4], "angle")?;
                let angle_rad = angle_deg * std::f64::consts::PI / 180.0;

                // PySCF convention for third atom:
                // - r1 (atom 1) is at origin
                // - a1 (atom 2) is along positive x-axis
                // - New atom is placed at bond distance from origin, at angle from x-axis in
                //   xz-plane
                // Angle is measured from the direction of a1 (positive x-axis)
                let cos_a = angle_rad.cos();
                let sin_a = angle_rad.sin();

                // Place at bond distance from r1 (origin), at angle from x-axis
                let x = r1_pos[0] + bond_bohr * cos_a;
                let z = bond_bohr * sin_a;
                [x, 0.0, z]
            },
            _ => {
                // Full z-matrix specification
                if parts.len() < 7 {
                    return cint_raise!(ParseError, "Invalid z-matrix line: '{}'. Need bond, angle, and dihedral.", line);
                }
                let r1_name = parts[1];
                let r1_idx = find_atom_index(&atom_order, r1_name)?;
                let r1_pos = atoms[r1_idx].coords;

                let a1_name = parts[3];
                let a1_idx = find_atom_index(&atom_order, a1_name)?;
                let a1_pos = atoms[a1_idx].coords;

                let d1_name = parts[5];
                let d1_idx = find_atom_index(&atom_order, d1_name)?;
                let d1_pos = atoms[d1_idx].coords;

                let bond: f64 = parse_float(parts[2], "bond length")?;
                let bond_bohr = if unit == Unit::Angstrom { bond * ANG_TO_BOHR } else { bond };

                let angle_deg: f64 = parse_float(parts[4], "angle")?;
                let angle_rad = angle_deg * std::f64::consts::PI / 180.0;

                let dihedral_deg: f64 = parse_float(parts[6], "dihedral")?;
                let dihedral_rad = dihedral_deg * std::f64::consts::PI / 180.0;

                zmatrix_to_cartesian(r1_pos, a1_pos, d1_pos, bond_bohr, angle_rad, dihedral_rad)
            },
        };

        let charge = if is_ghost {
            0.0
        } else {
            element_charge(&parsed_symbol).unwrap_or_else(|| panic!("Unknown element: {}", parsed_symbol)) as f64
        };

        atoms.push(AtomInfo { label: symbol.to_string(), symbol: parsed_symbol, charge, is_ghost, coords });
        atom_order.push(symbol.to_string());
    }

    Ok(atoms)
}

/// Parse a float value with error context.
fn parse_float(s: &str, context: &str) -> Result<f64, CIntError> {
    s.parse().map_err(|e: std::num::ParseFloatError| cint_error!(ParseError, "Failed to parse {} '{}': {}", context, s, e))
}

/// Find atom index by name in z-matrix.
fn find_atom_index(order: &[String], name: &str) -> Result<usize, CIntError> {
    // Try exact match first
    for (i, atom_name) in order.iter().enumerate() {
        if atom_name == name {
            return Ok(i);
        }
    }
    // Try numeric index (1-based)
    if let Ok(idx) = name.parse::<usize>() {
        if idx > 0 && idx <= order.len() {
            return Ok(idx - 1);
        }
    }
    cint_raise!(ParseError, "Atom '{}' not found in z-matrix", name)
}

/// Convert z-matrix parameters to Cartesian coordinates.
fn zmatrix_to_cartesian(
    r1: [f64; 3], // atom being bonded to
    a1: [f64; 3], // atom for angle reference
    d1: [f64; 3], // atom for dihedral reference
    bond: f64,
    angle: f64,
    dihedral: f64,
) -> [f64; 3] {
    // Vector from a1 to r1
    let v1 = [r1[0] - a1[0], r1[1] - a1[1], r1[2] - a1[2]];
    let v1_norm = norm3(v1);
    let v1_unit = [v1[0] / v1_norm, v1[1] / v1_norm, v1[2] / v1_norm];

    // Vector from d1 to a1
    let v2 = [a1[0] - d1[0], a1[1] - d1[1], a1[2] - d1[2]];
    let v2_norm = norm3(v2);
    let v2_unit = [v2[0] / v2_norm, v2[1] / v2_norm, v2[2] / v2_norm];

    // Normal vector (perpendicular to plane containing a1, r1, d1)
    let n = cross3(v1_unit, v2_unit);
    let n_norm = norm3(n);
    if n_norm < 1e-10 {
        // Collinear atoms, use arbitrary perpendicular
        let n_alt = cross3(v1_unit, [1.0, 0.0, 0.0]);
        let n_alt_norm = norm3(n_alt);
        if n_alt_norm < 1e-10 {
            let n_alt2 = cross3(v1_unit, [0.0, 1.0, 0.0]);
            let n_alt2_norm = norm3(n_alt2);
            [n_alt2[0] / n_alt2_norm, n_alt2[1] / n_alt2_norm, n_alt2[2] / n_alt2_norm]
        } else {
            [n_alt[0] / n_alt_norm, n_alt[1] / n_alt_norm, n_alt[2] / n_alt_norm]
        }
    } else {
        [n[0] / n_norm, n[1] / n_norm, n[2] / n_norm]
    };

    // Create local coordinate system
    // x-axis: perpendicular to v1 in the plane
    let x_local = cross3(n, v1_unit);
    let x_norm = norm3(x_local);
    let x_unit = [x_local[0] / x_norm, x_local[1] / x_norm, x_local[2] / x_norm];

    // Position in local coordinates (bond along -z, angle from z-axis)
    let sin_a = angle.sin();
    let cos_a = angle.cos();
    let sin_d = dihedral.sin();
    let cos_d = dihedral.cos();

    let x = bond * sin_a * cos_d;
    let y = bond * sin_a * sin_d;
    let z = -bond * cos_a;

    // Transform to global coordinates
    [
        r1[0] + x * x_unit[0] + y * n[0] + z * v1_unit[0],
        r1[1] + x * x_unit[1] + y * n[1] + z * v1_unit[1],
        r1[2] + x * x_unit[2] + y * n[2] + z * v1_unit[2],
    ]
}

/// Compute norm of 3D vector.
fn norm3(v: [f64; 3]) -> f64 {
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt()
}

/// Compute cross product of two 3D vectors.
fn cross3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]]
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_parse_atom_symbol() {
        assert_eq!(parse_atom_symbol("H"), ("H".to_string(), false));
        assert_eq!(parse_atom_symbol("H1"), ("H".to_string(), false));
        assert_eq!(parse_atom_symbol("O@2"), ("O".to_string(), false));
        assert_eq!(parse_atom_symbol("GHOST-H"), ("H".to_string(), true));
        assert_eq!(parse_atom_symbol("GHOST_O"), ("O".to_string(), true));
        assert_eq!(parse_atom_symbol("X-H"), ("H".to_string(), true));
        assert_eq!(parse_atom_symbol("X_O"), ("O".to_string(), true));
        assert_eq!(parse_atom_symbol("H@"), ("H".to_string(), true));
        assert_eq!(parse_atom_symbol("Au1"), ("Au".to_string(), false));
    }

    #[test]
    fn test_element_charge() {
        assert_eq!(element_charge("H"), Some(1));
        assert_eq!(element_charge("O"), Some(8));
        assert_eq!(element_charge("Au"), Some(79));
        assert_eq!(element_charge("X"), None);
    }

    #[test]
    fn test_parse_atom_string() {
        let atoms = parse_atom_string("H 0 0 0; O 0 0 1.2", Unit::Angstrom).unwrap();
        assert_eq!(atoms.len(), 2);
        assert_eq!(atoms[0].symbol, "H");
        assert_eq!(atoms[0].charge, 1.0);
        assert_relative_eq!(atoms[0].coords[2], 0.0);

        assert_eq!(atoms[1].symbol, "O");
        assert_eq!(atoms[1].charge, 8.0);
        assert_relative_eq!(atoms[1].coords[2], 1.2 * ANG_TO_BOHR);
    }

    #[test]
    fn test_parse_atom_string_ghost() {
        let atoms = parse_atom_string("H 0 0 0; GHOST-O 0 0 1.2", Unit::Angstrom).unwrap();
        assert_eq!(atoms.len(), 2);
        assert!(!atoms[0].is_ghost);
        assert!(atoms[1].is_ghost);
        assert_eq!(atoms[1].charge, 0.0);
    }

    #[test]
    fn test_parse_atom_ghost_formats() {
        // Test single ghost atom: GHOST-O
        let atoms1 = parse_atom_string("GHOST-O 0 0 1.2", Unit::Angstrom).unwrap();
        assert_eq!(atoms1.len(), 1);
        assert!(atoms1[0].is_ghost);
        assert_eq!(atoms1[0].charge, 0.0);

        // Test ghost atom with real atom (two atoms)
        let atoms2 = parse_atom_string("GHOST-O 0 0 1.2; O 0 0 0", Unit::Angstrom).unwrap();
        assert_eq!(atoms2.len(), 2);
        assert!(atoms2[0].is_ghost);
        assert!(!atoms2[1].is_ghost);
        assert_eq!(atoms2[0].charge, 0.0);
        assert_eq!(atoms2[1].charge, 8.0);

        // Test X_H format
        let atoms3 = parse_atom_string("X_H 0 0 1; H 0 0 0", Unit::Angstrom).unwrap();
        assert_eq!(atoms3.len(), 2);
        assert!(atoms3[0].is_ghost);
        assert!(!atoms3[1].is_ghost);
        assert_eq!(atoms3[0].charge, 0.0);
        assert_eq!(atoms3[1].charge, 1.0);

        // Test ghost.H format (lowercase with dot separator)
        let atoms4 = parse_atom_string("ghost.H 0 0 1; H 0 0 0", Unit::Angstrom).unwrap();
        assert_eq!(atoms4.len(), 2);
        assert!(atoms4[0].is_ghost);
        assert!(!atoms4[1].is_ghost);
        assert_eq!(atoms4[0].charge, 0.0);
        assert_eq!(atoms4[1].charge, 1.0);

        // Test multiple ghost atoms
        let atoms5 = parse_atom_string("GHOST-O 0 0 0; X_H 0 0 1; ghost.H 0 0 2; O 0 0 3; H 0 0 4; H 0 0 5", Unit::Angstrom).unwrap();
        assert_eq!(atoms5.len(), 6);
        assert!(atoms5[0].is_ghost); // GHOST-O
        assert!(atoms5[1].is_ghost); // X_H
        assert!(atoms5[2].is_ghost); // ghost.H
        assert!(!atoms5[3].is_ghost); // O
        assert!(!atoms5[4].is_ghost); // H
        assert!(!atoms5[5].is_ghost); // H
    }

    #[test]
    fn test_parse_atom_ghost_dot_format() {
        // Test ghost.H format (with dot) - standalone test
        let atoms = parse_atom_string("ghost.H 0 0 1", Unit::Angstrom).unwrap();
        assert_eq!(atoms.len(), 1);
        assert!(atoms[0].is_ghost);
        assert_eq!(atoms[0].charge, 0.0);
        assert_eq!(atoms[0].symbol, "H");
    }

    #[test]
    fn test_parse_atom_string_labeled() {
        let atoms = parse_atom_string("H1 0 0 0; H2 0 0 1", Unit::Bohr).unwrap();
        assert_eq!(atoms.len(), 2);
        assert_eq!(atoms[0].label, "H1");
        assert_eq!(atoms[0].symbol, "H");
        assert_eq!(atoms[1].label, "H2");
        assert_eq!(atoms[1].symbol, "H");
    }

    #[test]
    fn test_parse_zmatrix_water() {
        // Water molecule: O at origin, H at bond distance, H with angle
        // PySCF convention: H1 along x-axis, H2 at angle from x-axis in xz-plane
        let atoms = parse_zmatrix("O\nH 1 0.94\nH 1 0.94 2 104.5", Unit::Angstrom).unwrap();
        assert_eq!(atoms.len(), 3);
        assert_eq!(atoms[0].symbol, "O");

        // First H along x-axis (PySCF convention)
        assert_relative_eq!(atoms[1].coords[0], 0.94 * ANG_TO_BOHR, epsilon = 1e-6);
        assert_relative_eq!(atoms[1].coords[1], 0.0, epsilon = 1e-10);
        assert_relative_eq!(atoms[1].coords[2], 0.0, epsilon = 1e-10);

        // Second H at angle 104.5° from x-axis in xz-plane
        // Bond distance from O (origin), not from H1
        let bond_bohr = 0.94 * ANG_TO_BOHR;
        let angle_rad = 104.5 * std::f64::consts::PI / 180.0;
        let expected_x = bond_bohr * angle_rad.cos(); // negative (second quadrant)
        let expected_z = bond_bohr * angle_rad.sin(); // positive
        assert_relative_eq!(atoms[2].coords[0], expected_x, epsilon = 1e-6);
        assert_relative_eq!(atoms[2].coords[1], 0.0, epsilon = 1e-10);
        assert_relative_eq!(atoms[2].coords[2], expected_z, epsilon = 1e-6);
    }

    #[test]
    fn test_charge_to_symbol() {
        assert_eq!(charge_to_symbol(1).unwrap(), "H");
        assert_eq!(charge_to_symbol(8).unwrap(), "O");
        assert_eq!(charge_to_symbol(79).unwrap(), "Au");
        assert!(charge_to_symbol(0).is_err());
        // Z=118 (Oganesson) is the highest known element, Z=119 may return
        // placeholder in bse
    }

    #[test]
    fn test_unit_conversion() {
        let coords_ang = [1.0, 2.0, 3.0];
        let atom = AtomInfo::new("H", coords_ang, Unit::Angstrom);
        assert_relative_eq!(atom.coords[0], coords_ang[0] * ANG_TO_BOHR);
        assert_relative_eq!(atom.coords[1], coords_ang[1] * ANG_TO_BOHR);
        assert_relative_eq!(atom.coords[2], coords_ang[2] * ANG_TO_BOHR);

        let atom_bohr = AtomInfo::new("H", coords_ang, Unit::Bohr);
        assert_relative_eq!(atom_bohr.coords[0], 1.0);
    }
}
