use std::error::Error;
use std::fmt::{Debug, Display};

#[derive(Debug, Clone, PartialEq)]
pub enum CintError {
    IntegratorNotFound(String),
    IntegratorNotAvailable(String),
    RuntimeError(String),
    InvalidValue(String),
    Miscellaneous(String),
}

impl Display for CintError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(self, f)
    }
}

impl Error for CintError {}
