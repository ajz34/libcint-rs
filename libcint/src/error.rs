use derive_builder::UninitializedFieldError;
use std::error::Error;
use std::fmt::{Debug, Display};

#[derive(Debug, Clone)]
pub enum CIntError {
    IntegratorNotFound(String),
    IntegratorNotAvailable(String),
    RuntimeError(String),
    InvalidValue(String),
    Miscellaneous(String),
    UninitializedFieldError(UninitializedFieldError),
}

impl Display for CIntError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(self, f)
    }
}

impl Error for CIntError {}

impl From<UninitializedFieldError> for CIntError {
    fn from(err: UninitializedFieldError) -> Self {
        CIntError::UninitializedFieldError(err)
    }
}
