# libcint FFI Bindings and Wrapper

This project contains libcint (C language) FFI bindings, wrapper and build-from-source.

[libcint](https://github.com/sunqm/libcint) is a C library for GTO (gaussian-type orbital) electronic integral, can be applied in computational chemistry, and has already been applied extensively in [PySCF](https://github.com/pyscf/pyscf).

| Resources | Badges |
|--|--|
| Crate | [![Crate](https://img.shields.io/crates/v/libcint.svg)](https://crates.io/crates/libcint) |
| API Document | [![API Documentation](https://docs.rs/libcint/badge.svg)](https://docs.rs/libcint) |
| FFI Binding (libcint) | [v6.1.2](https://github.com/sunqm/libcint/tree/v6.1.2) |
| FFI Binding (qcint) | [v6.1.2](https://github.com/sunqm/qcint/tree/v6.1.2) |
| ECP Source Code (PySCF) | [v2.9.0](https://github.com/pyscf/pyscf/tree/v2.9.0) |

This crate is not official bindgen project, nither [libcint](https://github.com/sunqm/libcint), [PySCF](https://github.com/pyscf/pyscf), nor [REST](https://gitee.com/RESTGroup/rest). It is originally intended to be some pioneer work for possible future development of [rest_libcint](https://gitee.com/RESTGroup/rest_libcint) wrapper.

## Minimal Example

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

## Installation and Cargo Features

### Install with pre-compiled `libcint.so` (recommended)

If you have already compiled `libcint.so`, then put path of this shared library in `CINT_DIR` or `LD_LIBRARY_PATH` (or `REST_EXT_DIR`). Then you just use this library in your `Cargo.toml` file by

```toml
[dependencies]
libcint = { version = "0.1" }
```

### Install and also build-from-source

If you have not compiled `libcint.so` or `libcint.a`, then you are suggested to use this library by specifying some cargo features:

```toml
[dependencies]
libcint = { version = "0.1", features = ["build_from_source", "static"] }
```

The source code will be automatically downloaded from github, and cargo will handle the building process.

If access to github is not available, you can use environment variable `CINT_SRC` to specify source mirror of [sunqm/libcint](https://github.com/sunqm/libcint) or [sunqm/qcint](https://github.com/sunqm/qcint).

### Cargo features

- Default features: None of any listed below (use library provided by system or user, using [sunqm/libcint](https://github.com/sunqm/libcint), dynamic linking, without F12 and 4c1e support).
- `build_from_source`: Trigger of C language library libcint building. This performs by CMake; source code will be automatically downloaded from github (if environment variable `CINT_SRC` not specified).
- `static`: Use static library for linking. This will require static link `libcint.a`, and dynamic link `libquadmath.so`.
- `qcint`: Use [sunqm/qcint](https://github.com/sunqm/qcint) instead of [sunqm/libcint](https://github.com/sunqm/libcint). Some integrals will not be available if `qcint` does not supports that. This will also change URL source if cargo feature `build_from_source` specified.
- `with_f12`: Whether F12 integrals (`int2e_stg`, `int2e_yp`, etc.) are supported.
- `with_4c1e`: Whether 4c1e integrals (`int4c1e`, etc.) are supported. 

### Shell environment variables

- `CINT_DIR`, `LD_LIBRARY_PATH`, `REST_EXT_DIR`: Your compiled library path of `libcint.so` and `libcint.a`. This crate will try to find if this library is in directory root, or `directory_root/lib`. May override the library built by cargo feature `build_from_source`.
- `CINT_SRC`: Source of libcint or qcint (must be a git repository). Only works with cargo feature `build_from_source`.
- `CINT_VIR`: Version of libcint or qcint (e.g. `v6.1.2`, must starts with `v` prefix). Only works with cargo feature `build_from_source`.

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

This repository is licensed under Apache-2.0, the same to PySCF and REST.

- ECP part of code is directly copied and slightly modified from PySCF.
- Code for filling shell blocks of integrals is modified from an early version of rest_libcint.
