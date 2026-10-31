# libcint FFI Bindings and Wrapper

This project contains libcint (C language) FFI bindings, wrapper and build-from-source.

[libcint](https://github.com/sunqm/libcint) is a C library for GTO (gaussian-type orbital) electronic integral, can be applied in computational chemistry, and has already been applied extensively in [PySCF](https://github.com/pyscf/pyscf).

| Resources | Badges |
|--|--|
| Crate | [![Crate](https://img.shields.io/crates/v/libcint.svg)](https://crates.io/crates/libcint) |
| API Document | [![API Documentation](https://docs.rs/libcint/badge.svg)](https://docs.rs/libcint) |
| FFI Binding (libcint) | [v6.2.0](https://github.com/sunqm/libcint/tree/v6.2.0) |
| FFI Binding (qcint) | [v6.2.0](https://github.com/sunqm/qcint/tree/v6.2.0) |
| ECP Source Code (PySCF) | [v2.14.0](https://github.com/pyscf/pyscf/tree/v2.14.0) |

This crate is not official bindgen project, neither [libcint](https://github.com/sunqm/libcint), [PySCF](https://github.com/pyscf/pyscf), nor [REST](https://gitee.com/RESTGroup/rest). It is originally intended to be some pioneer work for possible future development of [rest_libcint](https://gitee.com/RESTGroup/rest_libcint) wrapper.

## Minimal Example

### Basic integral calculation

The most important function is [`CInt::integrate`](https://docs.rs/libcint/latest/libcint/cint/struct.CInt.html#method.integrate) and [`CInt::integrate_row_major`](https://docs.rs/libcint/latest/libcint/cint/struct.CInt.html#method.integrate_row_major). It is somehow similar to PySCF's `mol.intor(intor, aosym, shls_slice)`, but the user shall check output shape and strides. For more information on usage of crate `libcint`, we refer to [API Documentation](https://docs.rs/libcint).

```rust
use libcint::prelude::*;

// This is for testing and development purpose only.
// For actual usage, you should initialize `CInt` with your own molecule data.
let cint_data: CInt = init_h2o_def2_tzvp();

// int1e_ipkin: the integral of kinetic energy operator with derivative on the first orbital

// [mu, nu, comp], column-major, same data memory layout to PySCF's 2/3-center intor
let (out, shape): (Vec<f64>, Vec<usize>) = cint_data.integrate("int1e_ipkin", None, None).into();
assert_eq!(shape, [43, 43, 3]);

// [comp, mu, nu], row-major, same shape to PySCF intor
let (out, shape): (Vec<f64>, Vec<usize>) = cint_data.integrate_row_major("int1e_ipkin", None, None).into();
assert_eq!(shape, [3, 43, 43]);
```

### Building molecule from TOML/JSON (feature `bse`)

With the `bse` feature enabled, you can build a `CIntMol` object directly from TOML or JSON strings, similar to PySCF's `gto.M` format:

```rust
use libcint::prelude::*;

let toml = r#"
atom = "O 0 0 0; H 0 0 0.9572; H 0 0.9266 -0.239987"
basis = "STO-3G"
"#;

let mol = CIntMol::from_toml(toml).unwrap();
let ovlp = mol.cint.integrate("int1e_ovlp", None, None).out.unwrap();
```

For JSON input:

```rust
let json = r#"{"atom": "O 0 0 0; H 0 0 0.9572", "basis": "STO-3G"}"#;
let mol = CIntMol::from_json(json).unwrap();
```

Supported fields: `atom` (Cartesian/Z-matrix), `basis` (BSE name or custom), `ecp`, `unit`, `cart`, `ghost_ecp`. See [from_toml_docs.md](libcint/src/parse/from_toml_docs.md) for full specification.

### Important parts

- [`CInt`](https://docs.rs/libcint/latest/libcint/cint/struct.CInt.html): basic struct defining a molecule (atoms, shells, environment) for electron integrals.
- [`CInt::integrate`](https://docs.rs/libcint/latest/libcint/cint/struct.CInt.html#method.integrate): main integral driver, partly PySCF's `mol.intor` counterpart (column-major output).
- [`CInt::eval_gto`](https://docs.rs/libcint/latest/libcint/cint/struct.CInt.html#method.eval_gto): evaluate GTO values on grids (PySCF's `mol.eval_ao` counterpart).
- [`CInt::integrate_row_major`](https://docs.rs/libcint/latest/libcint/cint/struct.CInt.html#method.integrate_row_major): integral driver with row-major output (PySCF shape convention).
- [`CInt::integrate_cross`](https://docs.rs/libcint/latest/libcint/cint/struct.CInt.html#method.integrate_cross): integrals between several [`CInt`](https://docs.rs/libcint/latest/libcint/cint/struct.CInt.html) instances, useful with auxiliary basis sets.
- [`CInt::integrate_with_args`](https://docs.rs/libcint/latest/libcint/cint/struct.CInt.html#method.integrate_with_args), [`CInt::integrate_args_builder`](https://docs.rs/libcint/latest/libcint/cint/struct.CInt.html#method.integrate_args_builder): integral driver with an arguments builder (advanced options such as `out`, `row_major` and `cintopt`).
- [`CInt::get_cintopt`](https://docs.rs/libcint/latest/libcint/cint/struct.CInt.html#method.get_cintopt), [`CInt::get_cintopt_cross`](https://docs.rs/libcint/latest/libcint/cint/struct.CInt.html#method.get_cintopt_cross), [`CIntOpt`](https://docs.rs/libcint/latest/libcint/cint/struct.CIntOpt.html): build a reusable libcint optimizer (for a single molecule, or from the `mols` list of a cross integral) and pass it through `IntegrateArgs::cintopt` / `IntorCrossArgs::cintopt`, instead of rebuilding it for every call.
- `CIntMol` (feature `bse`): molecule builder from TOML/JSON, a weaker version of PySCF's `Mole` class.
- `ffi::cint_wrapper` and `ffi::cecp_wrapper`: supported integrators.

## Installation and Cargo Features

### Install with pre-compiled `libcint.so` (recommended)

If you have already compiled `libcint.so`, then put path of this shared library in `CINT_DIR` or `LD_LIBRARY_PATH` (or `REST_EXT_DIR`). Then you just use this library in your `Cargo.toml` file by

```toml
[dependencies]
libcint = { version = "0.5" }
```

### Install and also build-from-source

If you have not compiled `libcint.so` or `libcint.a`, then you are suggested to use this library by specifying some cargo features:

```toml
[dependencies]
libcint = { version = "0.5", features = ["build_from_source", "static"] }
```

The source code will be automatically downloaded from github, and cargo will handle the building process.

Build-from-source installs only libcint's library directory into Cargo's output prefix; downloaded sources and intermediate build files remain in the build directory.

If github is not accessible (offline builds, vendored sources), `CINT_SRC` also accepts a local source directory, a source archive file, or the URL of a source archive of [sunqm/libcint](https://github.com/sunqm/libcint) or [sunqm/qcint](https://github.com/sunqm/qcint); pin an archive with `CINT_SRC_SHA256`. A git URL or local repository still works as a source mirror (`CINT_VER` selects the tag). Relative paths are resolved against the `libcint-src` package directory, and a local git repository is cloned at `CINT_VER`, so use a non-git directory to build modified sources.

### Cargo features

- Default features: `bse` (molecule builder from TOML/JSON with Basis Set Exchange support).
- `bse`: Enable [`CIntMol::from_toml`] and [`CIntMol::from_json`] for building molecules from text input. Requires `bse` crate for basis set data.
- `api-v6_2`: Enable integrals added in libcint v6.2.0 (`int1e_ipr_origj`, `int1e_r_origj_ip`, `int1e_iprr_origj`, `int1e_rr_origj_ip`, `int3c1e_ipip1`, `int3c1e_ip1ip2`, `int3c1e_ipvip1`, `int3c1e_ipip2`). This feature is **not** enabled by default, so the default API stays compatible with libcint/qcint v6.1.2; enable it (e.g. `features = ["api-v6_2"]`) when using v6.2.0 and needing these integrals. With `qcint`, it additionally enables the integrals that qcint gained in v6.2.0. For a v6.1.2 source build, set `CINT_VER=v6.1.2`.
- `build_from_source`: Trigger of C language library libcint building. This performs by CMake; source code will be automatically downloaded from github (if environment variable `CINT_SRC` not specified).
- `static`: Use static library for linking. This will require static link `libcint.a`, and dynamic link `libquadmath.so` on targets other than macOS (unless `no_quadmath`).
- `qcint`: Use [sunqm/qcint](https://github.com/sunqm/qcint) instead of [sunqm/libcint](https://github.com/sunqm/libcint). Some integrals will not be available if `qcint` does not supports that. This will also change URL source if cargo feature `build_from_source` specified.
- `with_f12`: Whether F12 integrals (`int2e_stg`, `int2e_yp`, etc.) are supported.
- `with_4c1e`: Whether 4c1e integrals (`int4c1e`, etc.) are supported.
- `no_linkage`: Disable automatic linkage of `libcint` in `build.rs`. By default, this crate will not activate this feature and handle linkage properly. This feature is useful when you want to link `libcint` manually or use `libcint` in a different way (probably beneficial to conda, @jeanwsr). However, notice that ECP part will still be built anyway, which will require `libcint` to be linked properly by your self.
- `no_quadmath`: Do not link `libquadmath` with `static`, and with `build_from_source` build libcint without it (libcint then uses `long double` where it would use `__float128`). Needed for targets without libquadmath, such as aarch64 Linux.

### Shell environment variables

- `CINT_DIR`, `LD_LIBRARY_PATH`, `REST_EXT_DIR`: Your compiled library path of `libcint.so` and `libcint.a`. This crate will try to find if this library is in directory root, or `directory_root/lib`. Not searched with cargo feature `build_from_source`, which links the library it built.
- `CINT_SRC`: Source of libcint or qcint: a git repository (URL or local path), a source archive (local path, or URL ending in `.tar.gz`, `.tgz`, `.tar.bz2`, `.tar.xz` or `.zip`), or an unpacked source directory. Only works with cargo feature `build_from_source`.
- `CINT_VER`: Version of libcint or qcint (default `v6.2.0`; must start with `v` prefix), as a git tag. Only works with cargo feature `build_from_source` and a git `CINT_SRC`.
- `CINT_SRC_SHA256`: Expected SHA256 of the `CINT_SRC` archive; the build fails on a mismatch. Only works with cargo feature `build_from_source` and an archive `CINT_SRC`; a git source or an unpacked directory with this variable set is rejected.

## 50 lines RHF with Rust

This is a full example code to compute the RHF/def2-TZVP energy of H2O molecule, using
- This crate, `libcint`, as electronic integral library (corresponding to some parts of `pyscf.gto`);
- [RSTSR](https://github.com/RESTGroup/rstsr) as tensor library (corresponding to NumPy and SciPy).

You will see that except for import and preparation, the **code with core algorithms is 27 lines** (see also directory [examples/h2o_rhf](examples/h2o_rhf)). For comparison, the **same code in Python is 23 lines** (see also [pyscf_h2o_rhf.py](libcint/scripts/pyscf_h2o_rhf.py)).

```rust
use libcint::prelude::*;
use rstsr::prelude::*;

pub type Tsr = Tensor<f64, DeviceBLAS>;
pub type TsrView<'a> = TensorView<'a, f64, DeviceBLAS>;

/// Obtain integrals (in row-major, same to PySCF but reverse of libcint)
pub fn intor_row_major(cint_data: &CInt, intor: &str) -> Tsr {
    // use up all rayon available threads for tensor operations
    let device = DeviceBLAS::default();
    // intor, "s1", full_shls_slice
    let (out, shape) = cint_data.integrate(intor, None, None).into();
    // row-major by transposition of col-major shape
    rt::asarray((out, shape.f(), &device)).into_reverse_axes()
}

fn main() {
    let device = DeviceBLAS::default();
    // Assuming H2O/def2-TZVP data for `CInt` has been prepared
    let cint_data = init_h2o_def2_tzvp();

    /* #region rhf total energy algorithms */
    let atom_coords = {
        let coords = cint_data.atom_coords();
        let coords = coords.into_iter().flatten().collect::<Vec<f64>>();
        rt::asarray((coords, &device)).into_shape((-1, 3))
    };
    let atom_charges = rt::asarray((cint_data.atom_charges(), &device));
    let mut dist = rt::sci::cdist((atom_coords.view(), atom_coords.view()));
    dist.diagonal_mut(None).fill(f64::INFINITY);
    let eng_nuc = 0.5 * (&atom_charges * atom_charges.i((.., None)) / dist).sum();
    println!("Nuclear repulsion energy: {eng_nuc}");

    let hcore = intor_row_major(&cint_data, "int1e_kin") + intor_row_major(&cint_data, "int1e_nuc");
    let ovlp = intor_row_major(&cint_data, "int1e_ovlp");
    let int2e = intor_row_major(&cint_data, "int2e");
    let nocc = 5; // hardcoded for H2O, 5 occupied orbitals

    let mut dm = ovlp.zeros_like();
    for _ in 0..40 {
        // hardcoded SCF iterations
        let fock = &hcore + ((1.0_f64 * &int2e - 0.5_f64 * int2e.swapaxes(1, 2)) * &dm).sum_axes([-1, -2]);
        let (_, c) = rt::linalg::eigh((&fock, &ovlp)).into();
        dm = 2.0_f64 * c.i((.., ..nocc)) % c.i((.., ..nocc)).t();
    }
    let eng_scratch = &hcore + ((0.5_f64 * &int2e - 0.25_f64 * int2e.swapaxes(1, 2)) * &dm).sum_axes([-1, -2]);
    let eng_elec = (dm * &eng_scratch).sum();
    println!("Total elec energy: {eng_elec}");
    println!("Total RHF energy: {}", eng_nuc + eng_elec);
    /* #endregion */
}
```

## License

This repository is licensed under Apache-2.0, the same to PySCF.

- ECP part of code is directly copied and slightly modified from PySCF.
- Code for filling shell blocks of integrals is modified from an early version of rest_libcint.
