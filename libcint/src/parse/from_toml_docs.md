# Building Molecules from TOML

Build a molecule object from a TOML string for integral calculations.

## Quick Start

```rust
use libcint::prelude::*;

let toml = r#"
atom = "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987"
basis = "STO-3G"
"#;

let mol = CIntMol::from_toml(toml);

// Use mol.cint for integral calculations
let overlap = mol.cint.integrate("int1e_ovlp", None, None).out.unwrap();

let nao = mol.cint.nao();
// print indices
println!("{:<6}{:>8?}", "OVLP", (0..nao).collect::<Vec<usize>>());
for i in 0..nao {
      println!("AO {}: {:>8.4?}", i, &overlap[i * nao..(i + 1) * nao]);
}
```

This will give output:

```text
OVLP  [       0,        1,        2,        3,        4,        5,        6]
AO 0: [  1.0000,   0.2367,   0.0000,   0.0000,   0.0000,   0.0540,   0.0540]
AO 1: [  0.2367,   1.0000,   0.0000,   0.0000,   0.0000,   0.4748,   0.4748]
AO 2: [  0.0000,   0.0000,   1.0000,   0.0000,   0.0000,   0.0000,   0.0000]
AO 3: [  0.0000,   0.0000,   0.0000,   1.0000,   0.0000,   0.0000,   0.3809]
AO 4: [  0.0000,   0.0000,   0.0000,   0.0000,   1.0000,   0.3935,  -0.0987]
AO 5: [  0.0540,   0.4748,   0.0000,   0.0000,   0.3935,   1.0000,   0.2517]
AO 6: [  0.0540,   0.4748,   0.0000,   0.3809,  -0.0987,   0.2517,   1.0000]
```

## BSE Data Source

This library uses the rust crate [`bse`](https://docs.rs/bse/latest/bse/) to handle basis set data.

It is recommended to also install the `basis-set-exchange` python library.
This stores a local copy of the BSE database.
Alternatively, you can set `BSE_DATA_DIR` environment variable to point to a local copy of the BSE data.
The data can be cloned from [github repo](https://github.com/MolSSI-BSE/basis_set_exchange).

If encountered any problem, please consult [bse data source configuration](https://docs.rs/bse/latest/bse/#data-source-configuration).
This provides more details on how to configure the data source.

## TOML Input Fields

| Field | Type | Default | Description |
|-------|------|---------|-------------|
| `atom` | string | required | Atom coordinates (Cartesian or Z-matrix) |
| `basis` | string/dict | required | Basis set specification |
| `ecp` | string/dict | none | ECP/pseudopotential specification |
| `unit` | string | "angstrom" | Coordinate unit (angstrom/bohr) |
| `cart` | bool | false | Use Cartesian GTOs (default: spherical) |
| `ghost_ecp` | bool | false | Assign ECP to ghost atoms |
| `allow_empty_basis` | bool | false | Allow atoms without basis functions |

## Atom Definition

### Cartesian Coordinates

Atom coordinates can be specified in multiple formats:

**Semicolon-separated (inline):**
```toml
atom = "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987"
```

**Newline-separated (multi-line):**
```toml
atom = """
O  0.0  0.0  1.0
H  1.2  0.5  1.3
H -0.7  1.2  0.1
"""
```

**Comma-separated coordinates:**
```toml
atom = "O,0,0,0; H,0,0,0.9572"
```

Format per atom: `<label> <x> <y> <z>` or `<label> <x> <y> <z> <unit>`

The unit can be overridden per atom:
```toml
atom = "O 0 0 0 angstrom; H 0 0 1.89 bohr"
```

### Z-Matrix Format

Z-matrix defines atoms using internal coordinates (bond distances, angles, dihedrals):

```toml
atom = """
O
H 1 0.94
H 1 0.94 2 104.5
"""
```

**Format specification:**
- **Atom 1:** Just the element symbol (placed at origin)
- **Atom 2:** `<symbol> <ref_atom> <bond_distance>`
- **Atom 3:** `<symbol> <bond_ref> <bond> <angle_ref> <angle_degrees>`
- **Atom 4+:** `<symbol> <bond_ref> <bond> <angle_ref> <angle> <dihedral_ref> <dihedral_degrees>`

Reference atoms can be specified by:
- **Index:** 1-based integer (1 = first atom, 2 = second, etc.)
- **Label:** The atom symbol/label used earlier

**Example with labeled references:**
```toml
atom = """
C
H1 1 1.1
H  1 1.1 2 109.47
H2 1 1.1 2 109.47 3 120
H1 1 1.1 2 109.47 4 120
"""
```

**Z-Matrix Convention**

The Z-matrix conversion follows PySCF's convention:
- Atom 1 is placed at origin `(0, 0, 0)`
- Atom 2 is placed along positive x-axis `(bond, 0, 0)`
- Atom 3 lies in the xz-plane (`y = 0`)
- Angles are in **degrees** (not radians)
- Bond distances use the specified `unit`.
- Per-atom unit override is not supported in Z-matrix.
  All distances must be in the same unit.

This convention differs from some other software.
For example, Gaussian uses different orientation rules (atom 2 placed along z-axis).
The roundtrip `cart2zmat` → `zmat2cart` only recovers exact coordinates when the input follows this PySCF convention.

**Example: Methane in Z-matrix**
```toml
# Tetrahedral methane: H-C-H angle = 109.47°
atom = """
C
H 1 1.09
H 1 1.09 2 109.47
H 1 1.09 2 109.47 3 120
H 1 1.09 2 109.47 4 120
"""
basis = "STO-3G"
```

The expected Cartesian coordinates (in Bohr) will be:
```text
 0.000,  0.000,  0.000
 1.090,  0.000,  0.000
-0.363,  0.000,  1.028
-0.363,  0.890, -0.514
-0.363, -0.890, -0.514
```

### Ghost Atoms

Ghost atoms (used for embedding, counterpoise corrections, etc.) are specified with prefix markers.

**Prefix formats:**
- `GHOST-<symbol>` (recommended)
- `X-<symbol>` (alternative)
- `ghost.<symbol>` (PySCF-compatible)

The separator after `GHOST`/`X` must be non-alphanumeric (e.g., `-`, `_`, `.`).

**Examples:**
```toml
atom = "O 0 0 0; GHOST-H 0 0 1.0; H 0 0 -1.0"
atom = "O 0 0 0; X-H 0 0 1.0"
atom = "O 0 0 0; ghost.H 0 0 1.0"
```

Ghost atoms have:
- Nuclear charge = 0
- Same basis functions as the corresponding real atom.
  For example, `GHOST-H` has same basis as `H`.
- No ECP assigned (unless `ghost_ecp = true`)

### Atom Labels

Atoms can be labeled to distinguish different atoms of the same element.

**Label syntax:**
- Numeric suffix: `H1`, `H2`, `H-1`, `H@2`
- Any non-alphabetic suffix after the element symbol

Labels are used for:
- Assigning different basis sets to different atoms
- Z-matrix reference identification

**Basis parsing hierarchy:**

For each atom, three identifiers are extracted:
1. **label**: Full input string (e.g., `GHOST-H1`, `H@2`)
2. **identifier**: Element + suffix, without ghost prefix (e.g., `H1`, `H@2`)
3. **symbol**: Pure element symbol (e.g., `H`, `O`)

Basis matching uses this hierarchy (see [Basis Set Definition](#basis-set-definition)).

## Basis Set Definition

### Uniform Basis (Same for All Atoms)

Specify a basis set name from the Basis Set Exchange (BSE) library:

```toml
basis = "def2-TZVP"
```

### Per-Atom Basis (Dictionary Format)

Assign different basis sets to specific atoms:

**Using `[basis-custom]` table:**
```toml
basis = "custom"

[basis-custom]
O = "def2-TZVP"
H = "6-31G"
default = "STO-3G"
```

**Using inline dict:**
```toml
basis = { O = "def2-TZVP", H = "6-31G", default = "STO-3G" }
```

### Basis Matching Rules

When using dictionary format, basis is matched in this order:

1. **Exact label match**: `H1`, `GHOST-O`, etc.
2. **Identifier match**: `H1`, `O@2` (without ghost prefix)
3. **Element symbol match**: `H`, `O`, `C`
4. **Default**: If no match found, use `default` key
5. **Error**: If no match and no `default`, raise error

**Note:** Keys in dictionary is case-sensitive.
Defining two keys with the same element symbol but different cases will result in undefined behavior.

**Example:**
```toml
atom = """
O  0 0 0
H1 1 0 0
H2 0 1 0
"""
basis = "custom"

[basis-custom]
H1 = "aug-cc-pVQZ"    # Matches atom "H1" exactly
H  = "6-31G"          # Matches "H2" by element symbol
O  = "def2-TZVP"      # Matches "O" by element symbol
default = "STO-3G"    # Fallback for unmatched atoms
```

### Explicit Basis Definition

Instead of BSE names, provide explicit basis in various formats:

**Supported formats:**
- `nwchem` (NWChem format)
- `gaussian94` / `g94` (Gaussian format)
- `cp2k` (CP2K format)
- `gamess_us` (GAMESS-US format)
- `turbomole` (Turbomole format)
- `molcas` (Molcas/OpenMolcas format)
- `cfour` (CFOUR format)
- `molpro` (Molpro format)
- `crystal` (CRYSTAL format)
- `json` (BSE JSON format)

**NWChem format example:**
```toml
basis = """
BASIS "ao basis" SPHERICAL PRINT
H    S
      3.425250914E+00    0.1543289673E+00
      0.6239137298E+00   0.5353281423E+00
      0.1688554040E+00   0.4446345422E+00
O    S
      130.7093214E+00    0.1543289673E+00
      23.80886605E+00    0.5353281423E+00
      6.443608313E+00    0.4446345422E+00
O    SP
      5.033151319E+00   -0.9996722919E-01   0.1559162750E+00
      1.169596125E+00    0.3995128261E+00   0.6076837186E+00
      0.3803889600E+00   0.7001154689E+00   0.3919573931E+00
END
"""
```

**Gaussian94 format example:**
```toml
basis = """
H     0
S    3   1.00
      3.425250914D+00    0.1543289673D+00
      0.6239137298D+00   0.5353281423D+00
      0.1688554040D+00   0.4446345422D+00
****
O     0
S    3   1.00
      130.7093214D+03    0.1543289673D+00
      23.80886605D+00    0.5353281423D+00
      6.443608313D+00    0.4446345422D+00
SP   3   1.00
      5.033151319D+00   -0.9996722919D-01   0.1559162750D+00
      1.169596125D+00    0.3995128261D+00   0.6076837186D+00
      0.3803889600D+00   0.7001154689D+00   0.3919573931D+00
****
"""
```

**Mixed format with custom table:**
```toml
basis = "custom"

[basis-custom]
H1 = """
H     0
S    3   1.00
      3.425250914D+00    0.1543289673D+00
      0.6239137298D+00   0.5353281423D+00
      0.1688554040D+00   0.4446345422D+00
****
"""
default = "STO-3G"
```

### Differences from PySCF

1. **BSE data source**: This library uses `basis-set-exchange` Rust crate.
   PySCF uses `basis-set-exchange` Python library.
   For exact numerical comparison, use PySCF while fetch basis in NWChem format:
   `bse.get_basis("name", fmt="nwchem")`.

2. **Dict key precedence**: In this library, matching order is:
   `label → identifier → symbol → default`.
   PySCF uses: `label → symbol → default`.
   The identifier is not used in PySCF, so different may happen for labeled atoms and ghost atoms.

## ECP (Effective Core Potential)

### Basic ECP Usage

Assign ECP/pseudopotential for heavy elements:

```toml
basis = "def2-TZVP"
ecp = "def2-TZVP"    # ECP from same family
```

### Per-Atom ECP

Like basis, ECP can be assigned per atom:

```toml
basis = "def2-TZVP"
ecp = "custom"

[ecp-custom]
Au = "def2-TZVP"
Sb = "lanl2dz"
default = ""         # No ECP for other atoms
```

### ECP Precedence Rules

ECP data is resolved separately from basis:

1. First, try to get ECP from the `basis` specification.
2. Then, check `ecp` specification.
   If ECP data exists there, it **overrides** the ECP from `basis`.

**Example:**
```toml
# ECP from basis (def2-TZVP has built-in ECP for Sb)
basis = "def2-TZVP"
# Override Sb's ECP with lanl2dz
ecp = "custom"

[ecp-custom]
Sb = "lanl2dz"
```

### Ghost Atom ECP Behavior

By default, ECP is **not** assigned to ghost atoms.
Use `ghost_ecp = true` to enable ECP for ghost atoms:

```toml
atom = "GHOST-Au 0 0 0"
basis = "def2-TZVP"
ecp = "def2-TZVP"
ghost_ecp = true     # Ghost Au will have ECP
```

## Coordinate Units

### Default Unit

Coordinates are in **Angstrom** by default.

```toml
atom = "O 0 0 0"     # Units: Angstrom
basis = "STO-3G"
```

### Specify Unit

Set unit globally or per-atom:

```toml
# Global unit
unit = "bohr"
atom = "O 0 0 0; H 0 0 1.89"

# Per-atom override
atom = "O 0 0 0 angstrom; H 0 0 1.89 bohr"
```

**Accepted unit strings:**
- Angstrom: `angstrom`, `ang`, `ANG`, `Angstrom`
- Bohr: `bohr`, `BOHR`, `au`, `AU`, `a.u.`, `A.U.`

### Conversion

Internally, all coordinates are stored in Bohr (libcint convention).
Conversion factor: 1 Å = 0.529177210544 Bohr (NIST 2022 CODATA).

## Cartesian vs Spherical GTOs

### Default: Spherical Harmonics

By default, use spherical (pure) angular functions:
- s: 1 function
- p: 3 functions (px, py, pz become px, py, pz spherical)
- d: 5 functions (dxx, dyy, dzz, dxy, dxz, dyz → d0, d+1, d-1, d+2, d-2)
- f: 7 functions

```toml
basis = "6-31G"      # Spherical by default
```

### Enable Cartesian

For Cartesian (6-component d, 10-component f, etc.):

```toml
basis = "6-31G"
cart = true
```

## Empty Basis Handling

Some atoms (especially ghost atoms) may intentionally have no basis functions.
By default, this raises an error. Allow with:

```toml
atom = """
C 0 0 0
H1 1 0 0
X-H2 0 1 0
"""
basis = "custom"
allow_empty_basis = true

[basis-custom]
C = "6-31G"
H1 = "STO-3G"
H = ""            # Empty basis for X-H2 (ghost)
```

## Complete Example

A comprehensive TOML input file:

```toml
# Water molecule with custom basis
atom = """
O  0.000000  0.000000  0.000000
H1 0.758602  0.000000  0.504284
H2 -0.758602  0.000000  0.504284
"""

basis = "custom"
ecp = ""               # No ECP needed for H2O
unit = "angstrom"
cart = false           # Spherical GTOs

[basis-custom]
O = "def2-TZVP"
H1 = "6-31G*"
H2 = "6-31G*"
default = "STO-3G"

# Ghost atoms for counterpoise correction
# GHOST-H1 1.0 0 0  # Uncomment to add ghost
```

## JSON Alternative

For JSON input, use `CIntMol::from_json`:

```text
let json = r#"{
    "atom": "O 0 0 0; H 0 0 0.9572",
    "basis": "STO-3G",
    "unit": "angstrom"
}"#;
let mol = CIntMol::from_json(json);
```

**Note:** JSON does not support the `[basis-custom]` table mechanism.
For per-atom basis in JSON, use inline dictionary:

```text
let json = r#"{
    "atom": "O 0 0 0; H1 0 0 1; H2 0 0 -1",
    "basis": {"O": "def2-TZVP", "H1": "6-31G", "default": "STO-3G"}
}"#;
```

## Error Handling

Both `from_toml` and `from_json` panic on error. For fallible versions:

```text
// Returns Result<CIntMol, CIntError>
let result = CIntMol::from_toml_f(toml_str);
if let Ok(mol) = result {
    // Success
} else {
    // Handle error
}
```

Common errors:
- Unknown element symbol
- Invalid Z-matrix syntax
- Basis set not found in BSE
- No basis assigned to an atom (without `allow_empty_basis`)
- Custom table defined but `basis` not set to `"custom"`