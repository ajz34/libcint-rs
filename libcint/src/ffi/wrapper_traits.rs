use crate::cint_data::CintKind;
use crate::error::CintError;
use std::any::Any;
use std::ffi::{c_int, c_void};

/* #region integrator from cint */

pub trait Integrator {
    /// # Safety
    unsafe fn optimizer(
        &self,
        opt: *mut *mut c_void,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const f64,
    );
    /// # Safety
    unsafe fn integral_sph(
        &self,
        out: *mut f64,
        dims: *const c_int,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const f64,
        opt: *const c_void,
        cache: *mut f64,
    ) -> Result<c_int, CintError>;
    /// # Safety
    unsafe fn integral_cart(
        &self,
        out: *mut f64,
        dims: *const c_int,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const f64,
        opt: *const c_void,
        cache: *mut f64,
    ) -> Result<c_int, CintError>;
    /// # Safety
    unsafe fn integral_spinor(
        &self,
        out: *mut c_void,
        dims: *const c_int,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const f64,
        opt: *const c_void,
        cache: *mut f64,
    ) -> Result<c_int, CintError>;
    fn n_comp(&self) -> usize;
    fn n_spinor_comp(&self) -> usize;
    fn n_center(&self) -> usize;
    fn ng(&self) -> Vec<i32>;
    fn integrator_type(&self) -> &'static str;
    fn name(&self) -> &'static str;
    fn kind(&self) -> CintKind;
    fn as_any(&self) -> &dyn Any;
}

#[macro_export]
macro_rules! impl_integrator {
    (
        $intor: ident,
        $optimizer: ident,
        $integral_sph: ident,
        $integral_cart: ident,
        $integral_spinor: ident,
        $n_comp: expr,
        $n_spinor_comp: expr,
        $n_center: expr,
        $ng: expr,
        $integrator_type: literal,
        $name: literal,
        $kind: expr
    ) => {
        #[allow(non_camel_case_types)]
        pub struct $intor;
        impl Integrator for $intor {
            unsafe fn optimizer(
                &self,
                opt: *mut *mut c_void,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
            ) {
                unsafe { $optimizer(opt as _, atm, natm, bas, nbas, env) }
            }

            unsafe fn integral_sph(
                &self,
                out: *mut f64,
                dims: *const c_int,
                shls: *const c_int,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
                opt: *const c_void,
                cache: *mut f64,
            ) -> Result<c_int, CintError> {
                unsafe {
                    Ok($integral_sph(out, dims, shls, atm, natm, bas, nbas, env, opt as _, cache))
                }
            }

            unsafe fn integral_cart(
                &self,
                out: *mut f64,
                dims: *const c_int,
                shls: *const c_int,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
                opt: *const c_void,
                cache: *mut f64,
            ) -> Result<c_int, CintError> {
                unsafe {
                    Ok($integral_cart(out, dims, shls, atm, natm, bas, nbas, env, opt as _, cache))
                }
            }

            unsafe fn integral_spinor(
                &self,
                out: *mut c_void,
                dims: *const c_int,
                shls: *const c_int,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
                opt: *const c_void,
                cache: *mut f64,
            ) -> Result<c_int, CintError> {
                unsafe {
                    Ok($integral_spinor(
                        out as _, dims, shls, atm, natm, bas, nbas, env, opt as _, cache,
                    ))
                }
            }

            fn n_comp(&self) -> usize {
                $n_comp as usize
            }

            fn n_spinor_comp(&self) -> usize {
                $n_spinor_comp as usize
            }

            fn n_center(&self) -> usize {
                $n_center as usize
            }

            fn ng(&self) -> Vec<i32> {
                $ng
            }

            fn integrator_type(&self) -> &'static str {
                $integrator_type
            }

            fn name(&self) -> &'static str {
                $name
            }

            fn kind(&self) -> CintKind {
                $kind
            }

            fn as_any(&self) -> &dyn Any {
                self
            }
        }
    };
}

#[macro_export]
macro_rules! impl_integrator_sph_cart_only {
    (
        $intor: ident,
        $optimizer: ident,
        $integral_sph: ident,
        $integral_cart: ident,
        $n_comp: expr,
        $n_spinor_comp: expr,
        $n_center: expr,
        $ng: expr,
        $integrator_type: literal,
        $name: literal,
        $kind: expr
    ) => {
        #[allow(non_camel_case_types)]
        pub struct $intor;
        impl Integrator for $intor {
            unsafe fn optimizer(
                &self,
                opt: *mut *mut c_void,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
            ) {
                unsafe { $optimizer(opt as _, atm, natm, bas, nbas, env) }
            }

            unsafe fn integral_sph(
                &self,
                out: *mut f64,
                dims: *const c_int,
                shls: *const c_int,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
                opt: *const c_void,
                cache: *mut f64,
            ) -> Result<c_int, CintError> {
                unsafe {
                    Ok($integral_sph(out, dims, shls, atm, natm, bas, nbas, env, opt as _, cache))
                }
            }

            unsafe fn integral_cart(
                &self,
                out: *mut f64,
                dims: *const c_int,
                shls: *const c_int,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
                opt: *const c_void,
                cache: *mut f64,
            ) -> Result<c_int, CintError> {
                unsafe {
                    Ok($integral_cart(out, dims, shls, atm, natm, bas, nbas, env, opt as _, cache))
                }
            }

            unsafe fn integral_spinor(
                &self,
                _out: *mut c_void,
                _dims: *const c_int,
                _shls: *const c_int,
                _atm: *const c_int,
                _natm: c_int,
                _bas: *const c_int,
                _nbas: c_int,
                _env: *const f64,
                _opt: *const c_void,
                _cache: *mut f64,
            ) -> Result<c_int, CintError> {
                Err(CintError::IntegratorNotAvailable(format!(
                    "Integral spinor for {} is not available",
                    stringify!($name)
                )))
            }

            fn n_comp(&self) -> usize {
                $n_comp as usize
            }

            fn n_spinor_comp(&self) -> usize {
                $n_spinor_comp as usize
            }

            fn n_center(&self) -> usize {
                $n_center as usize
            }

            fn ng(&self) -> Vec<i32> {
                $ng
            }

            fn integrator_type(&self) -> &'static str {
                $integrator_type
            }

            fn name(&self) -> &'static str {
                $name
            }

            fn kind(&self) -> CintKind {
                $kind
            }

            fn as_any(&self) -> &dyn Any {
                self
            }
        }
    };
}

#[macro_export]
macro_rules! impl_integrator_sph_only {
    (
        $intor: ident,
        $optimizer: ident,
        $integral_sph: ident,
        $n_comp: expr,
        $n_spinor_comp: expr,
        $n_center: expr,
        $ng: expr,
        $integrator_type: literal,
        $name: literal,
        $kind: expr
    ) => {
        #[allow(non_camel_case_types)]
        pub struct $intor;
        impl Integrator for $intor {
            unsafe fn optimizer(
                &self,
                opt: *mut *mut c_void,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
            ) {
                unsafe { $optimizer(opt as _, atm, natm, bas, nbas, env) }
            }

            unsafe fn integral_sph(
                &self,
                out: *mut f64,
                dims: *const c_int,
                shls: *const c_int,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
                opt: *const c_void,
                cache: *mut f64,
            ) -> Result<c_int, CintError> {
                unsafe {
                    Ok($integral_sph(out, dims, shls, atm, natm, bas, nbas, env, opt as _, cache))
                }
            }

            unsafe fn integral_cart(
                &self,
                _out: *mut f64,
                _dims: *const c_int,
                _shls: *const c_int,
                _atm: *const c_int,
                _natm: c_int,
                _bas: *const c_int,
                _nbas: c_int,
                _env: *const f64,
                _opt: *const c_void,
                _cache: *mut f64,
            ) -> Result<c_int, CintError> {
                Err(CintError::IntegratorNotAvailable(format!(
                    "Integral cart for {} is not available",
                    stringify!($name)
                )))
            }

            unsafe fn integral_spinor(
                &self,
                _out: *mut c_void,
                _dims: *const c_int,
                _shls: *const c_int,
                _atm: *const c_int,
                _natm: c_int,
                _bas: *const c_int,
                _nbas: c_int,
                _env: *const f64,
                _opt: *const c_void,
                _cache: *mut f64,
            ) -> Result<c_int, CintError> {
                Err(CintError::IntegratorNotAvailable(format!(
                    "Integral spinor for {} is not available",
                    stringify!($name)
                )))
            }

            fn n_comp(&self) -> usize {
                $n_comp as usize
            }

            fn n_spinor_comp(&self) -> usize {
                $n_spinor_comp as usize
            }

            fn n_center(&self) -> usize {
                $n_center as usize
            }

            fn ng(&self) -> Vec<i32> {
                $ng
            }

            fn integrator_type(&self) -> &'static str {
                $integrator_type
            }

            fn name(&self) -> &'static str {
                $name
            }

            fn kind(&self) -> CintKind {
                $kind
            }

            fn as_any(&self) -> &dyn Any {
                self
            }
        }
    };
}

/* #endregion */
