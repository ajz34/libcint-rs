use crate::prelude::*;

use crate::ffi::cecp_wrapper::get_ecp_integrator;
use crate::ffi::cint_wrapper::get_cint_integrator;

/* #region structs */

#[derive(Debug, Clone, PartialEq)]
pub struct CInt {
    pub atm: Vec<[c_int; 6]>,    // ATM_SLOTS
    pub bas: Vec<[c_int; 8]>,    // BAS_SLOTS
    pub ecpbas: Vec<[c_int; 8]>, // BAS_SLOTS
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
/// - rare cases that integrator not available for the specified `CIntType`
///   (i.e., for `with_f12` integrators such as `int2e_stg`, the cart and spinor
///   versions are not available)
pub fn get_integrator(intor: &str) -> (Box<dyn Integrator>, CIntType) {
    get_integrator_f(intor).unwrap()
}

pub fn get_integrator_f(intor: &str) -> Result<(Box<dyn Integrator>, CIntType), CIntError> {
    // check if suffix is present
    let (intor, cint_type) = if let Some(intor) = intor.strip_suffix("_sph") {
        (intor, CIntType::Sph)
    } else if let Some(intor) = intor.strip_suffix("_cart") {
        (intor, CIntType::Cart)
    } else if let Some(intor) = intor.strip_suffix("_spinor") {
        (intor, CIntType::Spinor)
    } else {
        (intor, CIntType::default()) // default to sph
    };

    // explicitly check cint and ecp differently
    let intor = if let Some(intor) = get_cint_integrator(intor) {
        Ok(intor)
    } else if let Some(intor) = get_ecp_integrator(intor) {
        Ok(intor)
    } else {
        Err(CIntError::IntegratorNotFound(intor.to_string()))
    }?;

    // check if the integrator and cint_type is available
    match cint_type {
        CIntType::Sph => {
            if !intor.is_sph_available() {
                return Err(CIntError::IntegratorNotAvailable(format!(
                    "{} is not available for spherical integrals",
                    intor.name()
                )));
            }
        },
        CIntType::Cart => {
            if !intor.is_cart_available() {
                return Err(CIntError::IntegratorNotAvailable(format!(
                    "{} is not available for Cartesian integrals",
                    intor.name()
                )));
            }
        },
        CIntType::Spinor => {
            if !intor.is_spinor_available() {
                return Err(CIntError::IntegratorNotAvailable(format!(
                    "{} is not available for spinor integrals",
                    intor.name()
                )));
            }
        },
    }

    // should have checked anything, gracefully returns
    Ok((intor, cint_type))
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
    pub fn merge_ecp_data(&self) -> CInt {
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
    /// convenience. You can also call
    ///
    /// # Panics
    ///
    /// - integrator not found
    /// - rare cases that integrator not available for the specified `CIntType`
    ///   (i.e., for `with_f12` integrators such as `int2e_stg`, the cart and
    ///   spinor versions are not available)
    pub fn get_integrator(intor: &str) -> (Box<dyn Integrator>, CIntType) {
        get_integrator(intor)
    }

    pub fn get_integrator_f(intor: &str) -> Result<(Box<dyn Integrator>, CIntType), CIntError> {
        get_integrator_f(intor)
    }

    /// Obtain integrator optimizer.
    ///
    /// # Panics
    ///
    /// - integrator not found
    /// - `env` field not properly initialized for ECP integrator (should call
    ///   `merge_ecp_data` before this function)
    ///
    /// # See also
    ///
    /// PySCF `make_cintopt`
    pub fn make_optimizer(&self, intor: &str) -> CIntOptimizer {
        self.make_optimizer_f(intor).unwrap()
    }

    pub fn make_optimizer_f(&self, intor: &str) -> Result<CIntOptimizer, CIntError> {
        let (intor, _) = get_integrator_f(intor)?;
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
                let merged = if self.is_ecp_merged() { self } else { &self.merge_ecp_data() };
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
