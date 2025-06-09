use crate::ffi::cecp::ECPOpt;
use crate::ffi::cint;
use crate::ffi::cint::CINTOpt;
use std::any::Any;
use std::ffi::c_int;

/* #region integrator from cint */

pub trait Integrator {
    /// # Safety
    unsafe fn optimizer(
        &self,
        opt: *mut *mut CINTOpt,
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
        opt: *const CINTOpt,
        cache: *mut f64,
    ) -> c_int;
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
        opt: *const CINTOpt,
        cache: *mut f64,
    ) -> c_int;
    /// # Safety
    unsafe fn integral_spinor(
        &self,
        out: *mut cint::__BindgenComplex<f64>,
        dims: *const c_int,
        shls: *const c_int,
        atm: *const c_int,
        natm: c_int,
        bas: *const c_int,
        nbas: c_int,
        env: *const f64,
        opt: *const CINTOpt,
        cache: *mut f64,
    ) -> c_int;
    fn n_comp(&self) -> usize;
    fn n_spinor_comp(&self) -> usize;
    fn n_center(&self) -> usize;
    fn ng(&self) -> Vec<i32>;
    fn integrator_type(&self) -> &'static str;
    fn name(&self) -> &'static str;
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
        $name: literal
    ) => {
        #[allow(non_camel_case_types)]
        pub struct $intor;
        impl Integrator for $intor {
            unsafe fn optimizer(
                &self,
                opt: *mut *mut CINTOpt,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
            ) {
                unsafe { cint::$optimizer(opt, atm, natm, bas, nbas, env) }
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
                opt: *const CINTOpt,
                cache: *mut f64,
            ) -> c_int {
                unsafe {
                    cint::$integral_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache)
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
                opt: *const CINTOpt,
                cache: *mut f64,
            ) -> c_int {
                unsafe {
                    cint::$integral_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache)
                }
            }
            unsafe fn integral_spinor(
                &self,
                out: *mut cint::__BindgenComplex<f64>,
                dims: *const c_int,
                shls: *const c_int,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
                opt: *const CINTOpt,
                cache: *mut f64,
            ) -> c_int {
                unsafe {
                    cint::$integral_spinor(out, dims, shls, atm, natm, bas, nbas, env, opt, cache)
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
            fn as_any(&self) -> &dyn Any {
                self
            }
        }
    };
}

/* #endregion */

/* #region integrator from cecp */

pub trait ECPIntegrator {
    /// # Safety
    unsafe fn optimizer(
        &self,
        opt: *mut *mut ECPOpt,
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
        opt: *const ECPOpt,
        cache: *mut f64,
    ) -> c_int;
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
        opt: *const ECPOpt,
        cache: *mut f64,
    ) -> c_int;
    fn n_comp(&self) -> usize;
    fn integrator_type(&self) -> &'static str;
    fn name(&self) -> &'static str;
    fn as_any(&self) -> &dyn Any;
}

#[macro_export]
macro_rules! impl_ecpintegrator {
    (
        $integrator: ident,
        $optimizer: ident,
        $integral_sph: ident,
        $integral_cart: ident,
        $n_comp: expr,
        $integrator_type: literal,
        $name: literal
    ) => {
        #[allow(non_camel_case_types)]
        pub struct $integrator;
        impl ECPIntegrator for $integrator {
            unsafe fn optimizer(
                &self,
                opt: *mut *mut ECPOpt,
                atm: *const c_int,
                natm: c_int,
                bas: *const c_int,
                nbas: c_int,
                env: *const f64,
            ) {
                unsafe { cecp::$optimizer(opt, atm, natm, bas, nbas, env) }
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
                opt: *const ECPOpt,
                cache: *mut f64,
            ) -> c_int {
                unsafe {
                    cecp::$integral_sph(out, dims, shls, atm, natm, bas, nbas, env, opt, cache)
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
                opt: *const ECPOpt,
                cache: *mut f64,
            ) -> c_int {
                unsafe {
                    cecp::$integral_cart(out, dims, shls, atm, natm, bas, nbas, env, opt, cache)
                }
            }
            fn n_comp(&self) -> usize {
                $n_comp as usize
            }
            fn integrator_type(&self) -> &'static str {
                $integrator_type
            }
            fn name(&self) -> &'static str {
                $name
            }
            fn as_any(&self) -> &dyn Any {
                self
            }
        }
    };
}

/* #endregion */
