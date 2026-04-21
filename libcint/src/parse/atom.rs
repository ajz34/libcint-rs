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
    pub fn new(symbol: &str, coords: [f64; 3], unit: Unit) -> Result<Self, CIntError> {
        let (parsed_symbol, is_ghost) = parse_atom_symbol(symbol);
        let charge = if is_ghost { 0.0 } else { symbol_to_charge(&parsed_symbol)? as f64 };
        let coords_bohr =
            if unit == Unit::Angstrom { [coords[0] * ANG_TO_BOHR, coords[1] * ANG_TO_BOHR, coords[2] * ANG_TO_BOHR] } else { coords };
        Ok(AtomInfo { label: symbol.to_string(), symbol: parsed_symbol, charge, is_ghost, coords: coords_bohr })
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
pub fn symbol_to_charge(symbol: &str) -> Result<i32, CIntError> {
    element_Z_from_sym(symbol).ok_or_else(|| cint_error!(ParseError, "Unknown element: {symbol}"))
}

/// Convert atomic number to element symbol.
fn charge_to_symbol(z: i32) -> Result<String, CIntError> {
    element_sym_from_Z_with_normalize(z).ok_or_else(|| cint_error!(ParseError, "Invalid atomic number: {z}"))
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
        .collect::<Result<Vec<AtomInfo>, CIntError>>()?;

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
        return cint_raise!(ParseError, "Invalid atom line: '{line}'. Expected at least 4 parts.");
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
    let x: f64 = parse_float(parts[1], "x coordinate")?;
    let y: f64 = parse_float(parts[2], "y coordinate")?;
    let z: f64 = parse_float(parts[3], "z coordinate")?;

    // Check for unit override in 5th position
    let unit = if parts.len() > 4 {
        match parts[4].to_uppercase().as_str() {
            "ANG" | "ANGSTROM" => Unit::Angstrom,
            "BOHR" | "AU" | "A.U." => Unit::Bohr,
            _ => default_unit,
        }
    } else {
        default_unit
    };

    // Convert to Bohr if necessary
    let coords_bohr = if unit == Unit::Angstrom { [x * ANG_TO_BOHR, y * ANG_TO_BOHR, z * ANG_TO_BOHR] } else { [x, y, z] };

    let charge = if is_ghost { 0.0 } else { symbol_to_charge(&parsed_symbol)? as f64 };

    Ok(AtomInfo { label: symbol.to_string(), symbol: parsed_symbol, charge, is_ghost, coords: coords_bohr })
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
/// The first atom is placed at origin. The second atom is placed along x-axis.
/// The third atom is in the xz-plane. The fourth and subsequent atoms use
/// full z-matrix specification.
pub fn parse_zmatrix(s: &str, unit: Unit) -> Result<Vec<AtomInfo>, CIntError> {
    // Parse to ZmatEntry structs and atom labels
    let (entries, labels) = parse_zmatrix_to_entries(s)?;

    if entries.is_empty() {
        return Ok(vec![]);
    }

    // Convert to Cartesian coordinates
    let coords = zmat_entries_to_cart(&entries, unit)?;

    // Build AtomInfo structs
    let atoms: Vec<AtomInfo> = labels
        .iter()
        .zip(coords.iter())
        .map(|(label, coord)| {
            let (parsed_symbol, is_ghost) = parse_atom_symbol(label);
            let charge = if is_ghost { 0.0 } else { symbol_to_charge(&parsed_symbol)? as f64 };
            Ok(AtomInfo { label: label.clone(), symbol: parsed_symbol, charge, is_ghost, coords: *coord })
        })
        .collect::<Result<Vec<AtomInfo>, CIntError>>()?;

    Ok(atoms)
}

/// Parse z-matrix string to ZmatEntry structs and atom labels.
///
/// Returns a tuple of (ZmatEntry vector, label vector).
fn parse_zmatrix_to_entries(s: &str) -> Result<(Vec<ZmatEntry>, Vec<String>), CIntError> {
    let s = s.replace(';', "\n");

    let lines: Vec<&str> = s
        .lines()
        .filter(|line| {
            let trimmed = line.trim();
            !trimmed.is_empty() && !trimmed.starts_with('#')
        })
        .collect();

    if lines.is_empty() {
        return Ok((vec![], vec![]));
    }

    let mut entries: Vec<ZmatEntry> = vec![];
    let mut labels: Vec<String> = vec![];

    for (idx, line) in lines.iter().enumerate() {
        let parts: Vec<&str> = line.split_whitespace().collect();

        if parts.is_empty() {
            continue;
        }

        let label = parts[0].to_string();
        labels.push(label.clone());

        let entry = match idx {
            0 => {
                // First atom at origin
                if parts.len() > 1 {
                    eprintln!("Warning: Extra parameters in z-matrix line '{line}' will be ignored for first atom.");
                }
                ZmatEntry::first()
            },
            1 => {
                // Second atom: bond only
                match parts.len() {
                    ..3 => return cint_raise!(ParseError, "Invalid z-matrix line: '{line}'. Need bond length."),
                    3 => (),
                    4.. => eprintln!("Warning: Extra parameters in z-matrix line '{line}' will be ignored for second atom."),
                };
                // find_atom_index returns 0-indexed, keep as-is
                let bond_ref = find_atom_index(&labels[..idx], parts[1])?;
                let bond: f64 = parse_float(parts[2], "bond length")?;
                ZmatEntry { bond_ref, bond, angle_ref: 0, angle: 0.0, dihedral_ref: 0, dihedral: 0.0 }
            },
            2 => {
                // Third atom: bond and angle
                match parts.len() {
                    ..5 => return cint_raise!(ParseError, "Invalid z-matrix line: '{line}'. Need bond and angle."),
                    5 => (),
                    6.. => eprintln!("Warning: Extra parameters in z-matrix line '{line}' will be ignored for third atom."),
                };
                let bond_ref = find_atom_index(&labels[..idx], parts[1])?;
                let bond: f64 = parse_float(parts[2], "bond length")?;
                let angle_ref = find_atom_index(&labels[..idx], parts[3])?;
                let angle: f64 = parse_float(parts[4], "angle")?;
                ZmatEntry { bond_ref, bond, angle_ref, angle, dihedral_ref: 0, dihedral: 0.0 }
            },
            _ => {
                // Fourth+ atom: full specification
                match parts.len() {
                    ..7 => return cint_raise!(ParseError, "Invalid z-matrix line: '{line}'. Need bond, angle, and dihedral."),
                    7 => (),
                    8.. => eprintln!("Warning: Extra parameters in z-matrix line '{line}' will be ignored for fourth+ atom."),
                };
                let bond_ref = find_atom_index(&labels[..idx], parts[1])?;
                let bond: f64 = parse_float(parts[2], "bond length")?;
                let angle_ref = find_atom_index(&labels[..idx], parts[3])?;
                let angle: f64 = parse_float(parts[4], "angle")?;
                let dihedral_ref = find_atom_index(&labels[..idx], parts[5])?;
                let dihedral: f64 = parse_float(parts[6], "dihedral")?;
                ZmatEntry { bond_ref, bond, angle_ref, angle, dihedral_ref, dihedral }
            },
        };

        entries.push(entry);
    }

    Ok((entries, labels))
}

/// Convert ZmatEntry structs to Cartesian coordinates.
///
/// Handles arbitrary reference indices (0-indexed, not limited to atoms 0,1,2).
fn zmat_entries_to_cart(entries: &[ZmatEntry], unit: Unit) -> Result<Vec<[f64; 3]>, CIntError> {
    const TOL: f64 = 1e-7;

    if entries.is_empty() {
        return Ok(vec![]);
    }

    let mut coords = vec![[0.0, 0.0, 0.0]];

    for (i, entry) in entries.iter().enumerate().skip(1) {
        let bond_bohr = if unit == Unit::Angstrom { entry.bond * ANG_TO_BOHR } else { entry.bond };

        let coord = if i == 1 {
            // Second atom: place along x-axis from referenced atom
            let ref_pos = coords[entry.bond_ref];
            [ref_pos[0] + bond_bohr, ref_pos[1], ref_pos[2]]
        } else if i == 2 {
            // Third atom: use rotation matrix approach (no dihedral)
            let bonda_pos = coords[entry.bond_ref];
            let anga_pos = coords[entry.angle_ref];
            let angle_rad = entry.angle * std::f64::consts::PI / 180.0;

            // PySCF algorithm for third atom
            let v1 = [anga_pos[0] - bonda_pos[0], anga_pos[1] - bonda_pos[1], anga_pos[2] - bonda_pos[2]];
            let v1_norm = norm3(v1);

            if v1_norm < TOL {
                return cint_raise!(ParseError, "Reference atoms for z-matrix are too close");
            }

            let vecn = if !((v1[0].abs() < TOL) && (v1[1].abs() < TOL)) { cross3(v1, [0.0, 0.0, 1.0]) } else { [0.0, 0.0, 1.0] };

            let rmat = rotation_mat(vecn, angle_rad);
            let c = mat_vec_mul(&rmat, v1);
            let scale = bond_bohr / v1_norm;
            [bonda_pos[0] + c[0] * scale, bonda_pos[1] + c[1] * scale, bonda_pos[2] + c[2] * scale]
        } else {
            // Fourth+ atom: use zmatrix_to_cartesian
            let r1_pos = coords[entry.bond_ref];
            let a1_pos = coords[entry.angle_ref];
            let d1_pos = coords[entry.dihedral_ref];
            let angle_rad = entry.angle * std::f64::consts::PI / 180.0;
            let dihedral_rad = entry.dihedral * std::f64::consts::PI / 180.0;

            zmatrix_to_cartesian(r1_pos, a1_pos, d1_pos, bond_bohr, angle_rad, dihedral_rad)?
        };

        coords.push(coord);
    }

    Ok(coords)
}

/// Parse a float value with error context.
fn parse_float(s: &str, context: &str) -> Result<f64, CIntError> {
    s.parse().map_err(|e| cint_error!(ParseError, "Failed to parse {context} '{s}': {e}"))
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
    cint_raise!(ParseError, "Atom '{name}' not found in z-matrix")
}

/// Convert z-matrix parameters to Cartesian coordinates.
///
/// Uses the same algorithm as PySCF's from_zmatrix function.
fn zmatrix_to_cartesian(
    r1: [f64; 3], // atom being bonded to (bonda in PySCF)
    a1: [f64; 3], // atom for angle reference (anga in PySCF)
    d1: [f64; 3], // atom for dihedral reference (diha in PySCF)
    bond: f64,
    angle: f64,
    dihedral: f64,
) -> Result<[f64; 3], CIntError> {
    const TOL: f64 = 1e-7;

    // v1 = a1 - r1 (vector from r1 to a1, following PySCF convention)
    let v1 = [a1[0] - r1[0], a1[1] - r1[1], a1[2] - r1[2]];
    let v1_norm = norm3(v1);
    if v1_norm < TOL {
        return cint_raise!(ParseError, "Reference atoms for z-matrix ({r1:?}, {a1:?}) are too close");
    }
    let v1_unit = [v1[0] / v1_norm, v1[1] / v1_norm, v1[2] / v1_norm];

    // Handle edge cases: angle = 0 or angle = pi
    if angle < TOL {
        // Angle is 0, place atom along v1 direction
        return Ok([r1[0] + bond * v1_unit[0], r1[1] + bond * v1_unit[1], r1[2] + bond * v1_unit[2]]);
    }
    if std::f64::consts::PI - angle < TOL {
        // Angle is pi, place atom opposite to v1 direction
        return Ok([r1[0] - bond * v1_unit[0], r1[1] - bond * v1_unit[1], r1[2] - bond * v1_unit[2]]);
    }

    // v2 = d1 - a1 (vector from a1 to d1)
    let v2 = [d1[0] - a1[0], d1[1] - a1[1], d1[2] - a1[2]];

    // vecn = cross(v2, -v1) - normal to the plane
    let neg_v1 = [-v1_unit[0], -v1_unit[1], -v1_unit[2]];
    let vecn = cross3(v2, neg_v1);
    let vecn_norm = norm3(vecn);

    if vecn_norm < TOL {
        // Collinear atoms (d1, a1, r1 are on a line)
        // Use alternative normal vector
        let vecn_alt = if !((v1_unit[0].abs() < TOL) && (v1_unit[1].abs() < TOL)) {
            // v1 is not along z-axis, use cross(v1, z_axis)
            cross3(v1_unit, [0.0, 0.0, 1.0])
        } else {
            // v1 is along z-axis, use cross(v1, x_axis)
            cross3(v1_unit, [1.0, 0.0, 0.0])
        };
        let vecn_alt_norm = norm3(vecn_alt);
        if vecn_alt_norm < TOL {
            return cint_raise!(ParseError, "Cannot determine rotation axis for z-matrix");
        }
        let vecn_unit = [vecn_alt[0] / vecn_alt_norm, vecn_alt[1] / vecn_alt_norm, vecn_alt[2] / vecn_alt_norm];
        let rmat = rotation_mat(vecn_unit, angle);
        let c = mat_vec_mul(&rmat, v1_unit);
        Ok([r1[0] + bond * c[0], r1[1] + bond * c[1], r1[2] + bond * c[2]])
    } else {
        // Normal case: apply rotation
        // First rotate vecn around v1 by -dihedral
        let vecn_unit = [vecn[0] / vecn_norm, vecn[1] / vecn_norm, vecn[2] / vecn_norm];
        let rmat_dih = rotation_mat(v1_unit, -dihedral);
        let vecn_rot = mat_vec_mul(&rmat_dih, vecn_unit);

        // Then rotate v1 around the rotated vecn by angle
        let rmat_ang = rotation_mat(vecn_rot, angle);
        let c = mat_vec_mul(&rmat_ang, v1_unit);

        Ok([r1[0] + bond * c[0], r1[1] + bond * c[1], r1[2] + bond * c[2]])
    }
}

/// Rodrigues' rotation formula: rotate by angle theta around axis vec.
/// Returns 3x3 rotation matrix.
fn rotation_mat(vec: [f64; 3], theta: f64) -> [[f64; 3]; 3] {
    let norm = norm3(vec);
    let k = [vec[0] / norm, vec[1] / norm, vec[2] / norm];

    let c = theta.cos();
    let s = theta.sin();

    // uu = k * k^T (outer product)
    let uu = [[k[0] * k[0], k[0] * k[1], k[0] * k[2]], [k[1] * k[0], k[1] * k[1], k[1] * k[2]], [k[2] * k[0], k[2] * k[1], k[2] * k[2]]];

    // ux = cross product matrix of k
    let ux = [[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]];

    // r = c * I + s * ux + (1-c) * uu
    [
        [c + (1.0 - c) * uu[0][0] + s * ux[0][0], (1.0 - c) * uu[0][1] + s * ux[0][1], (1.0 - c) * uu[0][2] + s * ux[0][2]],
        [(1.0 - c) * uu[1][0] + s * ux[1][0], c + (1.0 - c) * uu[1][1] + s * ux[1][1], (1.0 - c) * uu[1][2] + s * ux[1][2]],
        [(1.0 - c) * uu[2][0] + s * ux[2][0], (1.0 - c) * uu[2][1] + s * ux[2][1], c + (1.0 - c) * uu[2][2] + s * ux[2][2]],
    ]
}

/// Multiply 3x3 matrix by 3-vector.
fn mat_vec_mul(m: &[[f64; 3]; 3], v: [f64; 3]) -> [f64; 3] {
    [
        m[0][0] * v[0] + m[0][1] * v[1] + m[0][2] * v[2],
        m[1][0] * v[0] + m[1][1] * v[1] + m[1][2] * v[2],
        m[2][0] * v[0] + m[2][1] * v[1] + m[2][2] * v[2],
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

/// Compute dot product of two 3D vectors.
fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Z-matrix entry for a single atom.
///
/// Internal references use 0-based indexing. When parsing/printing strings,
/// convert to/from 1-based indexing (standard z-matrix format).
#[derive(Debug, Clone, PartialEq)]
pub struct ZmatEntry {
    /// Reference atom for bond distance (0-based index).
    pub bond_ref: usize,
    /// Bond distance in the specified unit.
    pub bond: f64,
    /// Reference atom for angle (0-based index).
    pub angle_ref: usize,
    /// Angle in degrees.
    pub angle: f64,
    /// Reference atom for dihedral (0-based index).
    pub dihedral_ref: usize,
    /// Dihedral angle in degrees.
    pub dihedral: f64,
}

impl ZmatEntry {
    /// Create entry for first atom (at origin).
    pub fn first() -> Self {
        ZmatEntry { bond_ref: 0, bond: 0.0, angle_ref: 0, angle: 0.0, dihedral_ref: 0, dihedral: 0.0 }
    }

    /// Create entry for second atom (bond distance only, referencing atom 0).
    pub fn second(bond: f64) -> Self {
        ZmatEntry { bond_ref: 0, bond, angle_ref: 0, angle: 0.0, dihedral_ref: 0, dihedral: 0.0 }
    }

    /// Create entry for third atom (bond and angle, referencing atoms 0 and 1).
    pub fn third(bond: f64, angle: f64) -> Self {
        ZmatEntry { bond_ref: 0, bond, angle_ref: 1, angle, dihedral_ref: 0, dihedral: 0.0 }
    }

    /// Create entry for fourth+ atom (full specification, 0-indexed refs).
    pub fn full(bond_ref: usize, bond: f64, angle_ref: usize, angle: f64, dihedral_ref: usize, dihedral: f64) -> Self {
        ZmatEntry { bond_ref, bond, angle_ref, angle, dihedral_ref, dihedral }
    }
}

/// Convert Cartesian coordinates to z-matrix format.
///
/// Returns a vector of ZmatEntry structs representing the internal coordinates.
/// The coordinates are computed relative to atoms 1, 2, 3 as reference atoms.
///
/// # Arguments
/// * `coords` - Vector of Cartesian coordinates in Bohr
///
/// # Returns
/// Vector of ZmatEntry structs.
///
/// # Roundtrip Limitation
/// The `cart2zmat` -> `zmat2cart` roundtrip only exactly recovers the original
/// coordinates when the input follows the PySCF convention:
/// - Atom 1 at origin (0, 0, 0)
/// - Atom 2 along positive x-axis (y=0, z=0)
/// - Atom 3 in the xz-plane (y=0)
///
/// For arbitrary coordinates, the roundtrip preserves internal
/// distances/angles/dihedrals but loses the global orientation of the first 3
/// atoms. This is expected behavior matching PySCF, as z-matrix format
/// inherently loses orientation information.
pub fn cart2zmat(coords: &[[f64; 3]]) -> Vec<ZmatEntry> {
    const TOL: f64 = 1e-7;

    if coords.is_empty() {
        return vec![];
    }

    let mut zmat = vec![ZmatEntry::first()];

    if coords.len() > 1 {
        // Second atom: bond distance from first
        let r1 = [coords[1][0] - coords[0][0], coords[1][1] - coords[0][1], coords[1][2] - coords[0][2]];
        let bond = norm3(r1);
        zmat.push(ZmatEntry::second(bond));
    }

    if coords.len() > 2 {
        // Third atom: bond and angle
        let r1 = [coords[1][0] - coords[0][0], coords[1][1] - coords[0][1], coords[1][2] - coords[0][2]];
        let nr1 = norm3(r1);

        let r2 = [coords[2][0] - coords[0][0], coords[2][1] - coords[0][1], coords[2][2] - coords[0][2]];
        let nr2 = norm3(r2);

        let cos_a = dot3(r1, r2) / (nr1 * nr2);
        let angle = cos_a.acos() * 180.0 / std::f64::consts::PI;

        zmat.push(ZmatEntry::third(nr2, angle));
    }

    if coords.len() > 3 {
        // Fourth+ atoms: full z-matrix specification
        // Reference atoms: o0=coords[0], o1=coords[1], o2=coords[2] (0-indexed)
        let o0 = coords[0];
        let o1 = coords[1];
        let mut o2 = coords[2];
        let mut p2: usize = 2; // 0-indexed reference to atom 2

        for (k, c) in coords.iter().enumerate().skip(3) {
            // Bond distance from o0 (atom 0)
            let r0 = [c[0] - o0[0], c[1] - o0[1], c[2] - o0[2]];
            let nr0 = norm3(r0);

            // Angle: from direction o1-o0 (atom 1 to atom 0)
            let r1 = [o1[0] - o0[0], o1[1] - o0[1], o1[2] - o0[2]];
            let nr1 = norm3(r1);
            let cos_a1 = dot3(r0, r1) / (nr0 * nr1);
            let a1 = cos_a1.acos() * 180.0 / std::f64::consts::PI;

            // Dihedral: cross products to determine angle
            let b0 = cross3(r0, r1);
            let nb0 = norm3(b0);

            if nb0 < TOL {
                // o0, o1, c are collinear
                zmat.push(ZmatEntry::full(0, nr0, 1, a1, p2, 0.0));
            } else {
                let b1 = cross3([o2[0] - o0[0], o2[1] - o0[1], o2[2] - o0[2]], r1);
                let nb1 = norm3(b1);

                if nb1 < TOL {
                    // o0, o1, o2 are collinear
                    zmat.push(ZmatEntry::full(0, nr0, 1, a1, p2, 0.0));
                    // Update reference for next atoms
                    o2 = *c;
                    p2 = 3 + k; // 0-indexed: atom index = 3 + k - 1
                } else {
                    // Compute dihedral angle
                    let cos_a2 = dot3(b1, b0) / (nb0 * nb1);
                    let cross_b = cross3(b1, b0);
                    let sign = dot3(cross_b, r1);

                    let a2 = if sign < 0.0 {
                        cos_a2.acos() * 180.0 / std::f64::consts::PI
                    } else {
                        -cos_a2.acos() * 180.0 / std::f64::consts::PI
                    };

                    zmat.push(ZmatEntry::full(0, nr0, 1, a1, p2, a2));
                }
            }
        }
    }

    zmat
}

/// Convert z-matrix entries back to Cartesian coordinates.
///
/// Uses PySCF's algorithm for consistent results.
///
/// # Arguments
/// * `zmat` - Vector of ZmatEntry structs (0-indexed references)
/// * `unit` - Unit for bond distances (Angstrom or Bohr)
///
/// # Returns
/// Vector of Cartesian coordinates in Bohr.
pub fn zmat2cart(zmat: &[ZmatEntry], unit: Unit) -> Result<Vec<[f64; 3]>, CIntError> {
    if zmat.is_empty() {
        return Ok(vec![]);
    }

    let mut coords = vec![[0.0, 0.0, 0.0]];

    for (i, entry) in zmat.iter().enumerate() {
        if i == 0 {
            continue; // First atom at origin
        }

        let bond_bohr = if unit == Unit::Angstrom { entry.bond * ANG_TO_BOHR } else { entry.bond };

        let coord = match i {
            1 => {
                // Second atom along x-axis
                [bond_bohr, 0.0, 0.0]
            },
            2 => {
                // Third atom using PySCF's rotation matrix approach
                // References are 0-indexed
                let bonda_pos = coords[entry.bond_ref];
                let anga_pos = coords[entry.angle_ref];
                let angle_rad = entry.angle * std::f64::consts::PI / 180.0;

                // v1 = anga_pos - bonda_pos
                let v1 = [anga_pos[0] - bonda_pos[0], anga_pos[1] - bonda_pos[1], anga_pos[2] - bonda_pos[2]];
                let v1_norm = norm3(v1);

                if v1_norm < 1e-7 {
                    return cint_raise!(ParseError, "Reference atoms for z-matrix are too close");
                }

                // Create perpendicular vector vecn
                let vecn = if !((v1[0].abs() < 1e-7) && (v1[1].abs() < 1e-7)) { cross3(v1, [0.0, 0.0, 1.0]) } else { [0.0, 0.0, 1.0] };

                let rmat = rotation_mat(vecn, angle_rad);
                let c = mat_vec_mul(&rmat, v1);
                let scale = bond_bohr / v1_norm;
                [bonda_pos[0] + c[0] * scale, bonda_pos[1] + c[1] * scale, bonda_pos[2] + c[2] * scale]
            },
            _ => {
                // Fourth+ atom: full z-matrix (0-indexed references)
                let r1_idx = entry.bond_ref;
                let a1_idx = entry.angle_ref;
                let d1_idx = entry.dihedral_ref;

                if r1_idx >= coords.len() || a1_idx >= coords.len() || d1_idx >= coords.len() {
                    return cint_raise!(ParseError, "Invalid reference index in z-matrix");
                }

                let angle_rad = entry.angle * std::f64::consts::PI / 180.0;
                let dihedral_rad = entry.dihedral * std::f64::consts::PI / 180.0;

                zmatrix_to_cartesian(coords[r1_idx], coords[a1_idx], coords[d1_idx], bond_bohr, angle_rad, dihedral_rad)?
            },
        };

        coords.push(coord);
    }

    Ok(coords)
}

#[cfg(test)]
mod tests_llm_assist {
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
        assert_eq!(symbol_to_charge("H").unwrap(), 1);
        assert_eq!(symbol_to_charge("O").unwrap(), 8);
        assert_eq!(symbol_to_charge("Au").unwrap(), 79);
        assert!(symbol_to_charge("X").is_err());
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
        assert_relative_eq!(atoms[1].coords[1], 0.0, epsilon = 1e-6);
        assert_relative_eq!(atoms[1].coords[2], 0.0, epsilon = 1e-6);

        // Second H at angle 104.5° from x-axis in xz-plane
        // Bond distance from O (origin), not from H1
        let bond_bohr = 0.94 * ANG_TO_BOHR;
        let angle_rad = 104.5 * std::f64::consts::PI / 180.0;
        let expected_x = bond_bohr * angle_rad.cos(); // negative (second quadrant)
        let expected_z = bond_bohr * angle_rad.sin(); // positive
        assert_relative_eq!(atoms[2].coords[0], expected_x, epsilon = 1e-6);
        assert_relative_eq!(atoms[2].coords[1], 0.0, epsilon = 1e-6);
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
        let atom = AtomInfo::new("H", coords_ang, Unit::Angstrom).unwrap();
        assert_relative_eq!(atom.coords[0], coords_ang[0] * ANG_TO_BOHR);
        assert_relative_eq!(atom.coords[1], coords_ang[1] * ANG_TO_BOHR);
        assert_relative_eq!(atom.coords[2], coords_ang[2] * ANG_TO_BOHR);

        let atom_bohr = AtomInfo::new("H", coords_ang, Unit::Bohr).unwrap();
        assert_relative_eq!(atom_bohr.coords[0], 1.0);
    }
}

#[cfg(test)]
mod test_zmat {
    use approx::assert_relative_eq;

    use super::*;

    #[test]
    fn test_parse_zmatrix() {
        // Test z-matrix parsing (from PySCF docstring)
        // PySCF interprets bond lengths directly (in Bohr by default)
        // Expected coordinates from PySCF (in Bohr):
        // Atom 1: [0, 0, 0]
        // Atom 2: [2.67247631, 0, 0]
        // Atom 3: [2.67247631, 0, 3.27310166]
        // Atom 4: [0.53449526, 0.30859098, 2.83668811]
        let token = r#"
        H
        H 1 2.67247631453057
        H 1 4.22555607338457 2 50.7684795164077
        H 1 2.90305235726773 2 79.3904651036893 3 6.20854462618583
        "#;
        // Use Unit::Bohr to match PySCF's output
        let atoms = parse_zmatrix(token, Unit::Bohr).unwrap();

        // Expected coordinates from PySCF (in Bohr)
        let expected_coords =
            [[0.0, 0.0, 0.0], [2.67247631, 0.0, 0.0], [2.67247631, 0.0, 3.27310166], [0.53449526, 0.30859098, 2.83668811]];

        for (i, atom) in atoms.iter().enumerate() {
            println!("Atom {}: symbol={}, coords={:?}", i + 1, atom.symbol, atom.coords);
            assert_relative_eq!(atom.coords[0], expected_coords[i][0], epsilon = 1e-6);
            assert_relative_eq!(atom.coords[1], expected_coords[i][1], epsilon = 1e-6);
            assert_relative_eq!(atom.coords[2], expected_coords[i][2], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_parse_zmatrix_collinear() {
        // Test collinear atoms: first three atoms on a line
        // This tests the n_norm < TOL case
        // Atom1 at origin, Atom2 along x-axis, Atom3 also along x-axis (collinear)
        // Atom4 needs special handling for dihedral
        let token = "H\nH 1 1.0\nH 1 2.0 2 0.0\nH 1 1.5 2 90.0 3 45.0";
        let atoms = parse_zmatrix(token, Unit::Angstrom).unwrap();

        let expected_coords_angstrom = [
            [0.0, 0.0, 0.0], // Atom1
            [1.0, 0.0, 0.0], // Atom2
            [2.0, 0.0, 0.0], // Atom3 (collinear)
            [0.0, 0.0, 1.5], // Atom4 (special case)
        ];

        for (i, atom) in atoms.iter().enumerate() {
            println!("Atom {}: symbol={}, coords={:?}", i + 1, atom.symbol, atom.coords);
            assert_relative_eq!(atom.coords[0], expected_coords_angstrom[i][0] * ANG_TO_BOHR, epsilon = 1e-6);
            assert_relative_eq!(atom.coords[1], expected_coords_angstrom[i][1] * ANG_TO_BOHR, epsilon = 1e-6);
            assert_relative_eq!(atom.coords[2], expected_coords_angstrom[i][2] * ANG_TO_BOHR, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_parse_zmatrix_linear_molecule() {
        // Test linear molecule: CO2 (collinear)
        // C at origin, O1 along x-axis, O2 opposite direction (angle 180)
        let token = "C\nO 1 1.16\nO 1 1.16 2 180.0";
        let atoms = parse_zmatrix(token, Unit::Angstrom).unwrap();

        let expected_coords_angstrom = [
            [0.0, 0.0, 0.0],   // C at origin
            [1.16, 0.0, 0.0],  // O1 along x-axis
            [-1.16, 0.0, 0.0], // O2 opposite (angle 180)
        ];

        for (i, atom) in atoms.iter().enumerate() {
            println!("Atom {}: symbol={}, coords={:?}", i + 1, atom.symbol, atom.coords);
            assert_relative_eq!(atom.coords[0], expected_coords_angstrom[i][0] * ANG_TO_BOHR, epsilon = 1e-6);
            assert_relative_eq!(atom.coords[1], expected_coords_angstrom[i][1] * ANG_TO_BOHR, epsilon = 1e-6);
            assert_relative_eq!(atom.coords[2], expected_coords_angstrom[i][2] * ANG_TO_BOHR, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_parse_zmatrix_angle_zero() {
        // Test angle = 0 edge case: all atoms collinear along x-axis
        let token = "H\nH 1 1.0\nH 1 2.0 2 0.0\nH 1 3.0 2 0.0 3 0.0";
        let atoms = parse_zmatrix(token, Unit::Angstrom).unwrap();

        let expected_coords_angstrom = [
            [0.0, 0.0, 0.0], // Atom1 at origin
            [1.0, 0.0, 0.0], // Atom2 at x=1
            [2.0, 0.0, 0.0], // Atom3 at x=2 (angle 0)
            [3.0, 0.0, 0.0], // Atom4 at x=3 (angle 0, dihedral 0)
        ];

        for (i, atom) in atoms.iter().enumerate() {
            println!("Atom {}: symbol={}, coords={:?}", i + 1, atom.symbol, atom.coords);
            assert_relative_eq!(atom.coords[0], expected_coords_angstrom[i][0] * ANG_TO_BOHR, epsilon = 1e-6);
            assert_relative_eq!(atom.coords[1], expected_coords_angstrom[i][1] * ANG_TO_BOHR, epsilon = 1e-6);
            assert_relative_eq!(atom.coords[2], expected_coords_angstrom[i][2] * ANG_TO_BOHR, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_parse_zmatrix_angle_pi() {
        // Test angle = 180 (pi) edge case: atom placed opposite direction
        let token = "H\nH 1 1.0\nH 1 1.0 2 180.0";
        let atoms = parse_zmatrix(token, Unit::Angstrom).unwrap();

        let expected_coords_angstrom = [
            [0.0, 0.0, 0.0],  // Atom1 at origin
            [1.0, 0.0, 0.0],  // Atom2 along x-axis
            [-1.0, 0.0, 0.0], // Atom3 opposite (angle 180)
        ];

        for (i, atom) in atoms.iter().enumerate() {
            println!("Atom {}: symbol={}, coords={:?}", i + 1, atom.symbol, atom.coords);
            assert_relative_eq!(atom.coords[0], expected_coords_angstrom[i][0] * ANG_TO_BOHR, epsilon = 1e-6);
            assert_relative_eq!(atom.coords[1], expected_coords_angstrom[i][1] * ANG_TO_BOHR, epsilon = 1e-6);
            assert_relative_eq!(atom.coords[2], expected_coords_angstrom[i][2] * ANG_TO_BOHR, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_parse_zmatrix_dihedral_zero() {
        // Test dihedral = 0: all atoms in xz-plane (y = 0)
        let token = "H\nH 1 1.0\nH 1 1.0 2 90.0\nH 1 1.0 2 90.0 3 0.0";
        let atoms = parse_zmatrix(token, Unit::Angstrom).unwrap();

        let expected_coords_angstrom = [
            [0.0, 0.0, 0.0], // Atom1
            [1.0, 0.0, 0.0], // Atom2
            [0.0, 0.0, 1.0], // Atom3
            [0.0, 0.0, 1.0], // Atom4
        ];

        for (i, atom) in atoms.iter().enumerate() {
            println!("Atom {}: symbol={}, coords={:?}", i + 1, atom.symbol, atom.coords);
            assert_relative_eq!(atom.coords[0], expected_coords_angstrom[i][0] * ANG_TO_BOHR, epsilon = 1e-6);
            assert_relative_eq!(atom.coords[1], expected_coords_angstrom[i][1] * ANG_TO_BOHR, epsilon = 1e-6);
            assert_relative_eq!(atom.coords[2], expected_coords_angstrom[i][2] * ANG_TO_BOHR, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_parse_zmatrix_dihedral_180() {
        // Test dihedral = 180: Atom4 on opposite side of xz-plane from Atom3
        let token = "H\nH 1 1.0\nH 1 1.0 2 90.0\nH 1 1.0 2 90.0 3 180.0";
        let atoms = parse_zmatrix(token, Unit::Angstrom).unwrap();

        let expected_coords_angstrom = [
            [0.0, 0.0, 0.0],  // Atom1
            [1.0, 0.0, 0.0],  // Atom2
            [0.0, 0.0, 1.0],  // Atom3
            [0.0, 0.0, -1.0], // Atom4
        ];

        for (i, atom) in atoms.iter().enumerate() {
            println!("Atom {}: symbol={}, coords={:?}", i + 1, atom.symbol, atom.coords);
            assert_relative_eq!(atom.coords[0], expected_coords_angstrom[i][0] * ANG_TO_BOHR, epsilon = 1e-6);
            assert_relative_eq!(atom.coords[1], expected_coords_angstrom[i][1] * ANG_TO_BOHR, epsilon = 1e-6);
            assert_relative_eq!(atom.coords[2], expected_coords_angstrom[i][2] * ANG_TO_BOHR, epsilon = 1e-6);
        }
    }

    #[test]
    fn test_parse_zmatrix_bohr_unit() {
        // Test with Bohr unit: no conversion needed
        let token = "H\nH 1 1.0\nH 1 1.0 2 90.0";
        let atoms = parse_zmatrix(token, Unit::Bohr).unwrap();

        let expected_coords_bohr = [
            [0.0, 0.0, 0.0], // Atom1 at origin
            [1.0, 0.0, 0.0], // Atom2 along x-axis
            [0.0, 0.0, 1.0], // Atom3 at angle 90
        ];

        for (i, atom) in atoms.iter().enumerate() {
            println!("Atom {}: symbol={}, coords={:?}", i + 1, atom.symbol, atom.coords);
            assert_relative_eq!(atom.coords[0], expected_coords_bohr[i][0], epsilon = 1e-6);
            assert_relative_eq!(atom.coords[1], expected_coords_bohr[i][1], epsilon = 1e-6);
            assert_relative_eq!(atom.coords[2], expected_coords_bohr[i][2], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_cart2zmat_basic() {
        // Test basic cart2zmat conversion
        // Coordinates from PySCF docstring example (in Bohr)
        let coords: Vec<[f64; 3]> = vec![
            [0.0, 1.889726124565, 0.0],
            [0.0, 0.0, -1.889726124565],
            [1.889726124565, -1.889726124565, 0.0],
            [1.889726124565, 0.0, 1.133835674739],
        ];

        let zmat = cart2zmat(&coords);
        assert_eq!(zmat.len(), 4);

        // First atom
        assert_eq!(zmat[0].bond_ref, 0);

        // Second atom: bond from first (0-indexed: atom 0)
        assert_eq!(zmat[1].bond_ref, 0);
        assert_relative_eq!(zmat[1].bond, 2.67247631453057, epsilon = 1e-6);

        // Third atom: bond and angle (0-indexed: bond_ref=0, angle_ref=1)
        assert_eq!(zmat[2].bond_ref, 0);
        assert_relative_eq!(zmat[2].bond, 4.22555607338457, epsilon = 1e-6);
        assert_eq!(zmat[2].angle_ref, 1);
        assert_relative_eq!(zmat[2].angle, 50.7684795164077, epsilon = 1e-6);

        // Fourth atom: full specification (0-indexed: bond_ref=0, angle_ref=1,
        // dihedral_ref=2)
        assert_eq!(zmat[3].bond_ref, 0);
        assert_relative_eq!(zmat[3].bond, 2.90305235726773, epsilon = 1e-6);
        assert_eq!(zmat[3].angle_ref, 1);
        assert_relative_eq!(zmat[3].angle, 79.3904651036893, epsilon = 1e-6);
        assert_eq!(zmat[3].dihedral_ref, 2);
        assert_relative_eq!(zmat[3].dihedral.abs(), 6.20854462618583, epsilon = 1e-5);
    }

    #[test]
    fn test_cart2zmat_roundtrip_non_pyscf_convention() {
        // Test that cart2zmat -> zmat2cart does NOT recover original coordinates
        // when atoms don't follow PySCF convention.
        // This is EXPECTED behavior matching PySCF - z-matrix loses orientation info.
        // PySCF convention requires: Atom 1 at origin, Atom 2 along x-axis, Atom 3 in
        // xz-plane (y=0)
        let coords: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [1.0, 1.0, 0.0], // y != 0, NOT in xz-plane
            [0.5, 0.5, 1.0],
        ];

        let zmat = cart2zmat(&coords);
        let recovered = zmat2cart(&zmat, Unit::Bohr).unwrap();

        println!("Original coords: {:?}", coords);
        println!("Recovered coords: {:?}", recovered);

        // First two atoms should match (origin and x-axis)
        assert_relative_eq!(coords[0][0], recovered[0][0], epsilon = 1e-6);
        assert_relative_eq!(coords[0][1], recovered[0][1], epsilon = 1e-6);
        assert_relative_eq!(coords[0][2], recovered[0][2], epsilon = 1e-6);
        assert_relative_eq!(coords[1][0], recovered[1][0], epsilon = 1e-6);
        assert_relative_eq!(coords[1][1], recovered[1][1], epsilon = 1e-6);
        assert_relative_eq!(coords[1][2], recovered[1][2], epsilon = 1e-6);

        // Third atom: distance and angle preserved, but y becomes 0 (PySCF convention)
        // PySCF places atom 3 in xz-plane, so recovered[3][1] = 0
        let orig_dist = norm3([coords[2][0] - coords[0][0], coords[2][1] - coords[0][1], coords[2][2] - coords[0][2]]);
        let rec_dist = norm3([recovered[2][0] - recovered[0][0], recovered[2][1] - recovered[0][1], recovered[2][2] - recovered[0][2]]);
        assert_relative_eq!(orig_dist, rec_dist, epsilon = 1e-6);
        assert_relative_eq!(recovered[2][1], 0.0, epsilon = 1e-6); // PySCF convention

        // Fourth atom: internal coordinates preserved but orientation differs
        println!("Zmat entries for this case:");
        for (i, entry) in zmat.iter().enumerate() {
            println!(
                "Atom {}: bond_ref={}, bond={}, angle_ref={}, angle={}, dihedral_ref={}, dihedral={}",
                i + 1,
                entry.bond_ref,
                entry.bond,
                entry.angle_ref,
                entry.angle,
                entry.dihedral_ref,
                entry.dihedral
            );
        }
    }

    #[test]
    fn test_cart2zmat_roundtrip() {
        // Test that cart2zmat -> zmat2cart recovers original coordinates
        // IMPORTANT: Roundtrip only works for coordinates following PySCF convention:
        // - Atom 1 at origin
        // - Atom 2 along x-axis
        // - Atom 3 in xz-plane (y = 0)
        let coords: Vec<[f64; 3]> = vec![[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [1.0, 0.0, 1.0], [0.5, 0.5, 1.0]];

        let zmat = cart2zmat(&coords);
        let recovered = zmat2cart(&zmat, Unit::Bohr).unwrap();

        for (orig, rec) in coords.iter().zip(recovered.iter()) {
            println!("Original: {:?}, Recovered: {:?}", orig, rec);
            assert_relative_eq!(orig[0], rec[0], epsilon = 1e-6);
            assert_relative_eq!(orig[1], rec[1], epsilon = 1e-6);
            assert_relative_eq!(orig[2], rec[2], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_cart2zmat_roundtrip_water() {
        // Test roundtrip for water molecule
        let coords: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],                                                                                             // O at origin
            [0.94 * ANG_TO_BOHR, 0.0, 0.0],                                                                              // H1 along x
            [0.94 * ANG_TO_BOHR * 104.5_f64.to_radians().cos(), 0.0, 0.94 * ANG_TO_BOHR * 104.5_f64.to_radians().sin()], // H2
        ];

        let zmat = cart2zmat(&coords);
        let recovered = zmat2cart(&zmat, Unit::Bohr).unwrap();

        for (orig, rec) in coords.iter().zip(recovered.iter()) {
            assert_relative_eq!(orig[0], rec[0], epsilon = 1e-6);
            assert_relative_eq!(orig[1], rec[1], epsilon = 1e-6);
            assert_relative_eq!(orig[2], rec[2], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_cart2zmat_roundtrip_complex() {
        // Test roundtrip for a more complex geometry
        // This is the PySCF example
        let coords: Vec<[f64; 3]> =
            vec![[0.0, 0.0, 0.0], [2.67247631, 0.0, 0.0], [2.67247631, 0.0, 3.27310166], [0.53449526, 0.30859098, 2.83668811]];

        let zmat = cart2zmat(&coords);
        println!("Zmat entries:");
        for (i, entry) in zmat.iter().enumerate() {
            println!(
                "Atom {}: bond_ref={}, bond={}, angle_ref={}, angle={}, dihedral_ref={}, dihedral={}",
                i + 1,
                entry.bond_ref,
                entry.bond,
                entry.angle_ref,
                entry.angle,
                entry.dihedral_ref,
                entry.dihedral
            );
        }

        let recovered = zmat2cart(&zmat, Unit::Bohr).unwrap();

        for (i, (orig, rec)) in coords.iter().zip(recovered.iter()).enumerate() {
            println!("Atom {}: Original {:?}, Recovered {:?}", i + 1, orig, rec);
            assert_relative_eq!(orig[0], rec[0], epsilon = 1e-5);
            assert_relative_eq!(orig[1], rec[1], epsilon = 1e-5);
            assert_relative_eq!(orig[2], rec[2], epsilon = 1e-5);
        }
    }

    #[test]
    fn test_cart2zmat_roundtrip_angstrom() {
        // Test roundtrip with Angstrom unit
        // IMPORTANT: Atom 3 must be in xz-plane (y = 0) for roundtrip to work
        let coords_bohr: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.89, 0.0, 0.0],  // ~1 Angstrom
            [1.89, 0.0, 1.89], // in xz-plane (y = 0)
        ];

        let zmat = cart2zmat(&coords_bohr);
        // Bond distances are in Bohr from cart2zmat
        // Convert to Angstrom for zmat2cart
        let zmat_ang: Vec<ZmatEntry> = zmat
            .iter()
            .map(|e| ZmatEntry {
                bond_ref: e.bond_ref,
                bond: e.bond / ANG_TO_BOHR, // Convert to Angstrom
                angle_ref: e.angle_ref,
                angle: e.angle,
                dihedral_ref: e.dihedral_ref,
                dihedral: e.dihedral,
            })
            .collect();

        let recovered = zmat2cart(&zmat_ang, Unit::Angstrom).unwrap();

        for (orig, rec) in coords_bohr.iter().zip(recovered.iter()) {
            assert_relative_eq!(orig[0], rec[0], epsilon = 1e-6);
            assert_relative_eq!(orig[1], rec[1], epsilon = 1e-6);
            assert_relative_eq!(orig[2], rec[2], epsilon = 1e-6);
        }
    }

    #[test]
    fn test_cart2zmat_collinear() {
        // Test cart2zmat for collinear atoms
        let coords: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0], // collinear with first two
            [0.0, 1.0, 0.0], // perpendicular
        ];

        let zmat = cart2zmat(&coords);
        assert_eq!(zmat.len(), 4);

        // Third atom should have angle 0 (collinear)
        assert_relative_eq!(zmat[2].angle, 0.0, epsilon = 1e-6);

        // Fourth atom: dihedral might be 0 due to collinear reference
        println!("Fourth atom zmat: {:?}", zmat[3]);
    }

    #[test]
    fn test_parse_zmatrix_with_zmat_entry() {
        // Test that we can use ZmatEntry to reconstruct coordinates
        // References are 0-indexed
        let zmat = vec![
            ZmatEntry::first(),
            ZmatEntry::second(2.67247631453057),
            ZmatEntry::third(4.22555607338457, 50.7684795164077),
            ZmatEntry::full(0, 2.90305235726773, 1, 79.3904651036893, 2, 6.20854462618583),
        ];

        let coords = zmat2cart(&zmat, Unit::Angstrom).unwrap();

        // Expected coordinates (in Bohr) from PySCF
        let expected: Vec<[f64; 3]> = vec![
            [0.0, 0.0, 0.0],
            [2.67247631 * ANG_TO_BOHR, 0.0, 0.0],
            [2.67247631 * ANG_TO_BOHR, 0.0, 3.27310166 * ANG_TO_BOHR],
            [0.53449526 * ANG_TO_BOHR, 0.30859098 * ANG_TO_BOHR, 2.83668811 * ANG_TO_BOHR],
        ];

        for (rec, exp) in coords.iter().zip(expected.iter()) {
            println!("Recovered {:?}, Expected {:?}", rec, exp);
            assert_relative_eq!(rec[0], exp[0], epsilon = 1e-5);
            assert_relative_eq!(rec[1], exp[1], epsilon = 1e-5);
            assert_relative_eq!(rec[2], exp[2], epsilon = 1e-5);
        }
    }
}
