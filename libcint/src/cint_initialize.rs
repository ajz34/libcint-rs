//! Initializers for the CINT library.

use crate::prelude::*;
use serde::{Deserialize, Serialize};

/// Serde intermediate structure for conversion between `CInt` to other
/// serializable formats.
#[derive(Serialize, Deserialize, Debug, Clone)]
pub struct CIntSerdeIntermediate {
    pub _atm: Vec<[c_int; 6]>,
    pub _bas: Vec<[c_int; 8]>,
    pub _env: Vec<f64>,
    pub _ecpbas: Vec<[c_int; 8]>,

    #[serde(default = "default_cart")]
    pub cart: bool,
}

fn default_cart() -> bool {
    false
}

impl CIntSerdeIntermediate {
    /// Parse a JSON string (or path of JSON file) to `CIntSerdeIntermediate`.
    pub fn parse_from_json(token: &str) -> Result<CIntSerdeIntermediate, CIntError> {
        let result = serde_json::from_str::<CIntSerdeIntermediate>(token);
        match result {
            Ok(cint_data) => Ok(cint_data),
            Err(err) => {
                // let token to be the path of the file
                let err_msg = err.to_string();
                let file_content = std::fs::read_to_string(token).map_err(|_| CIntError::ParseError(err_msg))?;
                serde_json::from_str::<CIntSerdeIntermediate>(&file_content).map_err(|err| CIntError::ParseError(err.to_string()))
            },
        }
    }
}

impl From<CIntSerdeIntermediate> for CInt {
    fn from(cint_data: CIntSerdeIntermediate) -> Self {
        CInt {
            atm: cint_data._atm,
            bas: cint_data._bas,
            env: cint_data._env,
            ecpbas: cint_data._ecpbas,
            cint_type: match cint_data.cart {
                true => CIntType::Cartesian,
                false => CIntType::Spheric,
            },
        }
    }
}

impl From<CInt> for CIntSerdeIntermediate {
    fn from(cint_data: CInt) -> Self {
        CIntSerdeIntermediate {
            _atm: cint_data.atm,
            _bas: cint_data.bas,
            _env: cint_data.env,
            _ecpbas: cint_data.ecpbas,
            cart: matches!(cint_data.cint_type, CIntType::Cartesian),
        }
    }
}

/// Serde conversion between `CInt` and other serializable formats (JSON
/// currently).
impl CInt {
    /// Convert `CInt` to JSON string.
    pub fn to_json(&self) -> String {
        self.to_json_f().unwrap()
    }

    /// Convert JSON string (or the path of JSON file) to `CInt`.
    ///
    /// This function accepts the JSON string from PySCF's `mol.dumps()`.
    pub fn from_json(token: &str) -> CInt {
        CInt::from_json_f(token).unwrap()
    }

    pub fn to_json_f(&self) -> Result<String, CIntError> {
        let cint_data: CIntSerdeIntermediate = self.clone().into();
        serde_json::to_string(&cint_data).map_err(|err| CIntError::ParseError(err.to_string()))
    }

    pub fn from_json_f(token: &str) -> Result<CInt, CIntError> {
        let cint_data = CIntSerdeIntermediate::parse_from_json(token)?;
        Ok(cint_data.into())
    }
}

/// Implementation of legacy initializers for `CInt`.
impl CInt {
    /// Create a new empty `CInt` instance (not recommended).
    ///
    /// In most cases, you are encouraged to directly construct the `CInt`
    /// instance.
    ///
    /// If you are considering transforming PySCF's `mol` to `CInt`, you can use
    /// [`CInt::from_json`] to parse the JSON string from PySCF's `mol.
    /// dumps()`.
    #[allow(clippy::new_without_default)]
    pub fn new() -> Self {
        CInt { atm: vec![], bas: vec![], env: vec![], ecpbas: vec![], cint_type: CIntType::Spheric }
    }

    /// Legacy initializer for `CInt` (without ECP).
    #[allow(clippy::ptr_arg)]
    pub fn initial_r2c(&mut self, atm: &Vec<Vec<i32>>, natm: i32, bas: &Vec<Vec<i32>>, nbas: i32, env: &Vec<f64>) {
        assert!(atm.len() as i32 == natm);
        assert!(bas.len() as i32 == nbas);
        self.atm = atm.iter().map(|x| x.clone().try_into().unwrap()).collect();
        self.bas = bas.iter().map(|x| x.clone().try_into().unwrap()).collect();
        self.env = env.clone();
    }

    /// Legacy initializer for `CInt` (with ECP).
    #[allow(clippy::ptr_arg)]
    pub fn initial_r2c_with_ecp(
        &mut self,
        atm: &Vec<Vec<i32>>,
        natm: i32,
        bas: &Vec<Vec<i32>>,
        nbas: i32,
        ecp: &Vec<Vec<i32>>,
        necp: i32,
        env: &Vec<f64>,
    ) {
        self.initial_r2c(atm, natm, bas, nbas, env);
        assert!(ecp.len() as i32 == necp);
        self.ecpbas = ecp.iter().map(|x| x.clone().try_into().unwrap()).collect();
    }
}
