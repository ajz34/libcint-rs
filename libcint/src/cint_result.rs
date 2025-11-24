use crate::prelude::*;

/// Error for [`CInt`] operations.
#[derive(Debug, Clone)]
pub enum CIntErrorBase {
    IntegratorNotFound(String),
    IntegratorNotAvailable(String),
    RuntimeError(String),
    InvalidValue(String),
    UninitializedFieldError(UninitializedFieldError),
    ParseError(String),
    Miscellaneous(String),
}

#[derive(Debug)]
pub struct CIntError {
    pub kind: CIntErrorBase,
    pub backtrace: Option<Backtrace>,
}

/* #region CIntError impl */

impl Display for CIntErrorBase {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(self, f)
    }
}

impl Display for CIntError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(self, f)
    }
}

impl Error for CIntErrorBase {}
impl Error for CIntError {}

impl From<UninitializedFieldError> for CIntError {
    fn from(err: UninitializedFieldError) -> Self {
        CIntError { kind: CIntErrorBase::UninitializedFieldError(err), backtrace: Some(Backtrace::capture()) }
    }
}

pub trait CIntResultAPI<T> {
    fn cint_unwrap(self) -> T;
}

impl<T> CIntResultAPI<T> for Result<T, CIntError> {
    fn cint_unwrap(self) -> T {
        match self {
            Ok(val) => val,
            Err(err) => {
                let CIntError { kind, backtrace } = err;

                if let Some(backtrace) = backtrace {
                    std::eprintln!("\n====== CInt Backtrace ======\n{:}", backtrace);
                }
                panic!("CInt Error: {:?}", kind)
            },
        }
    }
}

#[macro_export]
macro_rules! cint_trace {
    () => {
        concat!(file!(), ":", line!(), ":", column!(), ": ")
    };
}

#[macro_export]
macro_rules! cint_error {
    ($errtype: ident, $($arg:tt)*) => {{
        use $crate::prelude::*;
        let mut s = String::new();
        write!(s, cint_trace!()).unwrap();
        write!(s, concat!("CIntError ", stringify!($errtype), ": ")).unwrap();
        write!(s, $($arg)*).unwrap();
        CIntError { kind: CIntErrorBase::$errtype(s), backtrace: Some(Backtrace::capture()) }
    }};
}

#[macro_export]
macro_rules! cint_raise {
    ($errtype: ident, $($arg:tt)*) => {{
        use $crate::prelude::*;
        let mut s = String::new();
        write!(s, cint_trace!()).unwrap();
        write!(s, concat!("CIntError ", stringify!($errtype), ": ")).unwrap();
        write!(s, $($arg)*).unwrap();
        Err(CIntError { kind: CIntErrorBase::$errtype(s), backtrace: Some(Backtrace::capture()) })
    }};
}

/* #endregion */
