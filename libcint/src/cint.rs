//! Main CInt struct definition, and some essential utilities.

use crate::prelude::*;

use crate::ffi::cecp_wrapper::get_ecp_integrator;
use crate::ffi::cint_wrapper::get_cint_integrator;

/* #region structs */

/// CInt data structure, which contains all the necessary information for
/// integral evaluation.
#[derive(Debug, Clone, PartialEq)]
pub struct CInt {
    /// Slots of atoms.
    ///
    /// Names of this field can be retrived in [`cint_ffi`].
    ///
    /// | Index | Name          | Description |
    /// |-------|---------------|-------------|
    /// | 0 | `CHARGE_OF`       | atomic charge (ECP core electrons excluded)       |
    /// | 1 | `PTR_COORD`       | coordinates location (pointer to `env` field)     |
    /// | 2 | `NUC_MOD_OF`      | nuclear mode of the atom                          |
    /// | 3 | `PTR_ZETA`        | pointer to zeta values (pointer to `env` field)   |
    /// | 4 | `PTR_FRAC_CHARGE` | fractional charge (pointer to `env` field)        |
    /// | 5 | `RESERVE_ATMSLOT` |
    /// | 6 | `ATM_SLOTS`       |
    ///
    /// - 2 - `NUC_MOD_OF` can be three kinds:
    ///     - 1 - `POINT_NUC` point charge, which is the most common case.
    ///     - 2 - `GAUSSIAN_NUC` nuclear charge is represented by a gaussian
    ///       distribution, with a zeta value (related to `PTR_ZETA`).
    ///     - 3 - `MM_NUC` is a point charge with fractional charge, which is
    ///       used in molecular mechanics (MM) calculations, with a fractional
    ///       value (related to `PTR_FRAC_CHARGE`).
    /// - 3 - Zeta value is atomic charge that is gaussian-distributed.
    ///   (`NUC_MOD_OF` refers to `GAUSSIAN_NUC`). In other cases, this value
    ///   should be zero.
    pub atm: Vec<[c_int; 6]>,

    /// Slot of shells (as minimal building block of contracted GTO).
    ///
    /// | Index | Name         | Description |
    /// |-------|--------------|-------------|
    /// | 0 | `ATOM_OF`        | 0-based index of corresponding atom                            |
    /// | 1 | `ANG_OF`         | angular momentum                                               |
    /// | 2 | `NPRIM_OF`       | number of primitive GTO in basis                               |
    /// | 3 | `NCTR_OF`        | number of contracted GTO in basis                              |
    /// | 4 | `KAPPA_OF`       | kappa for spinor GTO                                           |
    /// | 5 | `PTR_EXP`        | exponents of primitive GTOs (pointer to `env` field)           |
    /// | 6 | `PTR_COEFF`      | column-major contraction coefficients (pointer to `env` field) |
    /// | 7 | `RESERVE_BASLOT` |
    /// | 8 | `BAS_SLOTS`      |
    pub bas: Vec<[c_int; 8]>,

    /// Slot of shells for ECP (effective core potential).
    ///
    /// Other entries are the same as `bas` field, but with some additional:
    ///
    /// | Index | Name         | Description |
    /// |-------|--------------|-------------|
    /// | 1 | `ANG_OF`         | angular momentum (refers to angular of U_L) |
    /// | 3 | `RADI_POWER`     | number of maximum radial power to be summed |
    /// | 4 | `SO_TYPE_OF`     | spin-orb type                               |
    pub ecpbas: Vec<[c_int; 8]>,

    /// Slot of floating-point variables.
    ///
    /// | Index | Name         | Description |
    /// |-------|--------------|-------------|
    /// |  0 | `PTR_EXPCUTOFF`     | Overall cutoff for integral prescreening, value needs to be ln (threshold) |
    /// |  1 | `PTR_COMMON_ORIG`   | R_C of (r-R_C) in dipole, GIAO operators                                   |
    /// |  4 | `PTR_RINV_ORIG`     | R_O in 1/(r-R_O)                                                           |
    /// |  7 | `PTR_RINV_ZETA`     | ZETA parameter for Gaussian charge distribution (Gaussian nuclear model)   |
    /// |  8 | `PTR_RANGE_OMEGA`   | omega parameter in range-separated coulomb operator (>0 for LR, <0 for SR) |
    /// |  9 | `PTR_F12_ZETA`      | Yukawa potential and Slater-type geminal e^(-zeta r)                       |
    /// | 10 | `PTR_GTG_ZETA`      | Gaussian type geminal e^(-zeta r^2)                                        |
    /// | 11 | `NGRIDS`            | For hybrid integrals with grids                                            |
    /// | 12 | `PTR_GRIDS`         | Location of grids for `int1e_grids` or variants                            |
    /// | 17 | `AS_RINV_ORIG_ATOM` | Position of atom to be used as origin in 1/(r-R_O) for ECP derivatives     |
    /// | 18 | `AS_ECPBAS_OFFSET`  | Offset of ECP basis in `bas` field, used for ECP integrals                 |
    /// | 19 | `AS_NECPBAS`        | Number of ECP shells in `ecpbas` field                                     |
    /// | 20 | `PTR_ENV_START`     | Start of data                                                              |
    pub env: Vec<f64>,

    /// Type of integral.
    ///
    /// - `Spheric`: spherical harmonics, which is the most common case.
    /// - `Cartesian`: cartesian coordinates.
    /// - `Spinor`: spinor integrals, which are used in relativistic
    ///   calculations.
    pub cint_type: CIntType,
}

/// Distinguish between general integral or ECP (effective core potential)
/// integral.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CIntKind {
    /// General integral type, which is used for most of the integrals.
    Int,
    /// ECP integral type, which is used for effective core potential integrals.
    Ecp,
}

/// Integral types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CIntType {
    #[default]
    /// Spherical integrals, output double float `f64`.
    Spheric,
    /// Cartesian integrals, output double float `f64`.
    Cartesian,
    /// Spinor integrals, output `Complex<f64>`, and should be called by
    /// functions with `_spinor` suffix.
    Spinor,
}

/// Symmetry of the integral.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CIntSymm {
    #[default]
    /// No symmetry.
    S1,
    /// Symmetric for the first two indices (in F-contiguous).
    S2ij,
    /// Symmetric for the last two indices in 4-center integrals (in
    /// F-contiguous).
    S2kl,
    /// Symmetric for the first and last indices in 4-center integrals (in
    /// F-contiguous).
    S4,
    /// Symmetric for all indices in 4-center integrals (in F-contiguous).
    S8,
}

/// CInt optimizer, which is used to optimize the integral evaluation.
///
/// Note that currently, this is not intended to be used directly by users.
#[derive(Debug, PartialEq)]
pub enum CIntOptimizer {
    Int(*mut CINTOpt),
    Ecp(*mut ECPOpt),
}

/// Output of the integral evaluation.
///
/// You can use `.into()` to convert it to a tuple of `(Vec<F>, Vec<usize>)`,
/// where `F` is `f64` for sph and cart integrals, and `Complex<f64>` for spinor
/// integrals.
pub struct CIntOutput<F> {
    pub out: Option<Vec<F>>,
    pub shape: Vec<usize>,
}

#[derive(Builder, Debug)]
#[builder(pattern = "owned", build_fn(error = "CIntError"))]
pub struct IntorArgs<'a, 'b, 'c, F> {
    /// Integrator name, such as "int1e_ovlp", "int2e", etc.
    #[builder(setter(into))]
    pub intor: &'a str,

    /// Shell slices for evaluating sub-tensor of the total integral.
    #[builder(default = "&[]")]
    pub shls_slice: &'b [[usize; 2]],

    /// Symmetry of the integral. This will affect shape of output tensor.
    #[builder(default, setter(into))]
    pub aosym: CIntSymm,

    /// Output buffer for the integral evaluation.
    ///
    /// If this is `None`, the integral engine will generate a buffer for
    /// output.
    #[builder(default, setter(strip_option))]
    pub out: Option<&'c mut [F]>,
}

/* #endregion */

/* #region get_integrator */

/// Obtain integrator by name.
///
/// # Panics
///
/// - integrator not found
pub fn get_integrator(intor: &str) -> Box<dyn Integrator> {
    get_integrator_f(intor).unwrap()
}

pub fn get_integrator_f(intor: &str) -> Result<Box<dyn Integrator>, CIntError> {
    // explicitly check cint and ecp differently
    let intor = if let Some(intor) = get_cint_integrator(&intor.to_lowercase()) {
        Ok(intor)
    } else if let Some(intor) = get_ecp_integrator(&intor.to_lowercase()) {
        Ok(intor)
    } else {
        Err(CIntError::IntegratorNotFound(intor.to_string()))
    }?;

    // should have checked anything, gracefully returns
    Ok(intor)
}

/* #endregion */

/* #region CIntSymm impl */

impl From<&str> for CIntSymm {
    #[inline]
    fn from(aosym: &str) -> Self {
        match aosym {
            "S1" | "s1" => CIntSymm::S1,
            "S2ij" | "s2ij" => CIntSymm::S2ij,
            "S2kl" | "s2kl" => CIntSymm::S2kl,
            "S4" | "s4" => CIntSymm::S4,
            "S8" | "s8" => CIntSymm::S8,
            _ => panic!("Unknown symmetry: {aosym}"),
        }
    }
}

impl From<Option<&str>> for CIntSymm {
    #[inline]
    fn from(aosym: Option<&str>) -> Self {
        aosym.unwrap_or("s1").into()
    }
}

/* #endregion */

/* #region CIntOptimizer impl */

unsafe impl Send for CIntOptimizer {}
unsafe impl Sync for CIntOptimizer {}

impl Drop for CIntOptimizer {
    fn drop(&mut self) {
        match self {
            Self::Int(opt) => unsafe { cint_ffi::CINTdel_optimizer(opt) },
            Self::Ecp(opt) => unsafe { cecp_ffi::ECPdel_optimizer(opt) },
        }
    }
}

impl CIntOptimizer {
    #[inline]
    pub fn as_ptr(&self) -> *const c_void {
        match self {
            Self::Int(opt) => *opt as *const c_void,
            Self::Ecp(opt) => *opt as *const c_void,
        }
    }
}

/* #endregion */

/* #region CIntOutput impl */

impl<F> From<CIntOutput<F>> for (Vec<F>, Vec<usize>) {
    fn from(output: CIntOutput<F>) -> Self {
        if output.out.is_none() {
            panic!(
                "You have called integral with output buffer, and should not call `into` to get the integral buffer, but use your mutable reference instead."
            );
        }
        (output.out.unwrap_or_default(), output.shape)
    }
}

/* #endregion */

/// Merge ECP data for integral evaluation.
impl CInt {
    /// Creates a new `CIntData` instance, with ECP integral information
    /// written to `bas` field, and properly initializes values in `env` field.
    ///
    /// Should be called before integral calculations. Be careful this function
    /// should not called twice in consecutive.
    pub fn merge_ecpbas(&self) -> CInt {
        println!("merging ECP data into CIntData...");
        if self.is_ecp_merged() {
            self.clone()
        } else {
            let mut merged = self.clone();
            merged.bas.extend_from_slice(&self.ecpbas);
            merged.env[cecp_ffi::AS_ECPBAS_OFFSET as usize] = self.bas.len() as f64;
            merged.env[cecp_ffi::AS_NECPBAS as usize] = self.ecpbas.len() as f64;
            merged
        }
    }

    /// Check whether the `CInt` instance has been merged with ECP data.
    pub fn is_ecp_merged(&self) -> bool {
        self.env[cecp_ffi::AS_ECPBAS_OFFSET as usize] != 0.0
    }
}

/// Obtaining integrator and argument builder.
impl CInt {
    /// Obtain integrator by name.
    ///
    /// This function does not require `CInt` instance. Make it here only for
    /// convenience.
    ///
    /// # Panics
    ///
    /// - integrator not found
    pub fn get_integrator(intor: &str) -> Box<dyn Integrator> {
        get_integrator(intor)
    }

    pub fn get_integrator_f(intor: &str) -> Result<Box<dyn Integrator>, CIntError> {
        get_integrator_f(intor)
    }

    pub fn intor_args_builder(&self) -> IntorArgsBuilder<'static, 'static, 'static, f64> {
        IntorArgsBuilder::default()
    }

    pub fn intor_args_builder_spinor(
        &self,
    ) -> IntorArgsBuilder<'static, 'static, 'static, Complex<f64>> {
        IntorArgsBuilder::default()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cint_data() {
        let cint_data = init_h2o_def2_tzvp();
        let opt = cint_data.optimizer("int2e");
        println!("int2e opt: {opt:?}");
        if let CIntOptimizer::Int(opt) = opt {
            println!("int2e opt inner {:?}", unsafe { opt.read() });
        }
    }
}
