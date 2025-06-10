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
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CIntKind {
    Int,
    Ecp,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CIntType {
    #[default]
    Sph,
    Cart,
    Spinor,
}

#[derive(Debug, PartialEq)]
pub enum CIntOptimizer {
    Int(NonNull<*mut CINTOpt>),
    Ecp(NonNull<*mut ECPOpt>),
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
    let intor = if let Some(intor) = get_cint_integrator(intor) {
        Ok(intor)
    } else if let Some(intor) = get_ecp_integrator(intor) {
        Ok(intor)
    } else {
        Err(CIntError::IntegratorNotFound(intor.to_string()))
    }?;

    // should have checked anything, gracefully returns
    Ok(intor)
}

/* #endregion */

impl Drop for CIntOptimizer {
    fn drop(&mut self) {
        match self {
            Self::Int(opt) => unsafe { cint_ffi::CINTdel_optimizer(opt.as_mut()) },
            Self::Ecp(opt) => unsafe { cecp_ffi::ECPdel_optimizer(opt.as_mut()) },
        }
    }
}

/// Merge ECP data for integral evaluation.
impl CInt {
    /// Creates a new `CIntData` instance, with ECP integral information
    /// written to `bas` field, and properly initializes values in `env` field.
    ///
    /// Should be called before integral calculations. Be careful this function
    /// should not called twice in consecutive.
    pub fn merge_ecpbas(&self) -> CInt {
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

/// Obtaining integrator and optimizer.
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

    /// Obtain integrator optimizer.
    ///
    /// # Panics
    ///
    /// - integrator not found
    /// - `env` field not properly initialized for ECP integrator (should call
    ///   `merge_ecpbas` before this function)
    ///
    /// # See also
    ///
    /// PySCF `make_cintopt`
    pub fn make_optimizer(&self, intor: &str) -> CIntOptimizer {
        self.make_optimizer_f(intor).unwrap()
    }

    pub fn make_optimizer_f(&self, intor: &str) -> Result<CIntOptimizer, CIntError> {
        let intor = get_integrator_f(intor)?;
        match intor.kind() {
            CIntKind::Int => {
                let atm_ptr = self.atm.as_ptr() as *const c_int;
                let bas_ptr = self.bas.as_ptr() as *const c_int;
                let env_ptr = self.env.as_ptr();
                let n_atm = self.atm.len() as c_int;
                let n_bas = self.bas.len() as c_int;
                let mut c_opt_ptr: *mut CINTOpt = null_mut();
                unsafe {
                    intor.optimizer(
                        &mut c_opt_ptr as *mut *mut CINTOpt as *mut *mut c_void,
                        atm_ptr,
                        n_atm,
                        bas_ptr,
                        n_bas,
                        env_ptr,
                    );
                };
                Ok(CIntOptimizer::Int(NonNull::new(&mut c_opt_ptr).unwrap()))
            },
            CIntKind::Ecp => {
                // Check if this cint_data has been merged with ECP data
                let merged = if self.is_ecp_merged() { self } else { &self.merge_ecpbas() };
                let atm_ptr = merged.atm.as_ptr() as *const c_int;
                let bas_ptr = merged.bas.as_ptr() as *const c_int;
                let env_ptr = merged.env.as_ptr();
                let n_atm = merged.atm.len() as c_int;
                let n_bas = merged.bas.len() as c_int;
                let mut c_opt_ptr: *mut ECPOpt = null_mut();
                unsafe {
                    intor.optimizer(
                        &mut c_opt_ptr as *mut *mut ECPOpt as *mut *mut c_void,
                        atm_ptr,
                        n_atm,
                        bas_ptr,
                        n_bas,
                        env_ptr,
                    );
                };
                Ok(CIntOptimizer::Ecp(NonNull::new(&mut c_opt_ptr).unwrap()))
            },
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cint_data() {
        let cint_data = init_h2o_def2_tzvp();
        let opt = cint_data.make_optimizer("int2e");
        println!("int2e opt: {opt:?}");
        if let CIntOptimizer::Int(opt) = opt {
            println!("int2e opt inner {:?}", unsafe { opt.read() });
        }
    }
}
