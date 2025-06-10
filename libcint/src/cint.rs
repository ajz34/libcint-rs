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
                let merged = self.merge_ecp_data();
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
        let cint_data = initialize();
        let opt = cint_data.make_optimizer("int2e");
        println!("int2e opt: {opt:?}");
        if let CIntOptimizer::Int(opt) = opt {
            println!("int2e opt inner {:?}", unsafe { opt.read() });
        }
    }

    fn initialize() -> CInt {
        // mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5",
        // basis="def2-TZVP").build()
        let atm = vec![[8, 20, 1, 23, 0, 0], [1, 24, 1, 27, 0, 0], [1, 28, 1, 31, 0, 0]];
        let bas = vec![
            [0, 0, 6, 1, 0, 44, 50, 0],
            [0, 0, 2, 1, 0, 56, 58, 0],
            [0, 0, 1, 1, 0, 60, 61, 0],
            [0, 0, 1, 1, 0, 62, 63, 0],
            [0, 0, 1, 1, 0, 64, 65, 0],
            [0, 1, 4, 1, 0, 66, 70, 0],
            [0, 1, 1, 1, 0, 74, 75, 0],
            [0, 1, 1, 1, 0, 76, 77, 0],
            [0, 2, 1, 1, 0, 78, 79, 0],
            [0, 2, 1, 1, 0, 80, 81, 0],
            [0, 3, 1, 1, 0, 82, 83, 0],
            [1, 0, 3, 1, 0, 32, 35, 0],
            [1, 0, 1, 1, 0, 38, 39, 0],
            [1, 0, 1, 1, 0, 40, 41, 0],
            [1, 1, 1, 1, 0, 42, 43, 0],
            [2, 0, 3, 1, 0, 32, 35, 0],
            [2, 0, 1, 1, 0, 38, 39, 0],
            [2, 0, 1, 1, 0, 40, 41, 0],
            [2, 1, 1, 1, 0, 42, 43, 0],
        ];
        let c_env = vec![
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            1.7763425570911580e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            0.0000000000000000e+00,
            -4.4476065664656128e-01,
            0.0000000000000000e+00,
            1.7197618551510188e+00,
            0.0000000000000000e+00,
            3.4061340999999999e+01,
            5.1235746000000004e+00,
            1.1646626000000000e+00,
            9.0618446120248586e-01,
            1.6354784928239057e+00,
            2.4145128304249659e+00,
            3.2723041000000003e-01,
            1.0930883523645869e+00,
            1.0307241000000000e-01,
            4.5959135109675275e-01,
            8.0000000000000004e-01,
            2.2072263710762661e+00,
            2.7032382631000000e+04,
            4.0523871392000001e+03,
            9.2232722709999996e+02,
            2.6124070989000001e+02,
            8.5354641350999998e+01,
            3.1035035245000000e+01,
            3.0481181169845928e+00,
            5.6914576328642115e+00,
            9.7338835744432526e+00,
            1.5238733819733028e+01,
            2.0843228934131737e+01,
            2.2391049059992991e+01,
            1.2260860728000001e+01,
            4.9987076005000004e+00,
            1.0568131135849375e+01,
            3.3391469496791393e+00,
            1.1703108158000000e+00,
            2.8427648592056753e+00,
            4.6474740994000002e-01,
            1.4220922112658689e+00,
            1.8504536357000001e-01,
            7.1280983010446131e-01,
            6.3274954801000000e+01,
            1.4627049379000001e+01,
            4.4501223455999996e+00,
            1.5275799646999999e+00,
            6.2570323747894276e+00,
            6.9268656235998423e+00,
            6.0323599265415284e+00,
            3.5035168827833356e+00,
            5.2935117942999999e-01,
            1.3172379939563448e+00,
            1.7478421270000000e-01,
            3.2969483673949351e-01,
            2.3140000000000001e+00,
            1.1328313432935008e+01,
            6.4500000000000002e-01,
            1.2113199965714336e+00,
            1.4279999999999999e+00,
            4.3969226782656516e+00,
        ];

        CInt { atm, bas, ecpbas: vec![], env: c_env }
    }
}
