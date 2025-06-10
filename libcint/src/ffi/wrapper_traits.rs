use crate::prelude::*;

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
        opt: *const c_void,
        cache: *mut f64,
    ) -> c_int;
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
    ) -> c_int;
    fn is_sph_available(&self) -> bool;
    fn is_cart_available(&self) -> bool;
    fn is_spinor_available(&self) -> bool;
    fn n_comp(&self) -> usize;
    fn n_spinor_comp(&self) -> usize;
    fn n_center(&self) -> usize;
    fn ng(&self) -> Vec<i32>;
    fn integrator_category(&self) -> &'static str;
    fn name(&self) -> &'static str;
    fn kind(&self) -> CIntKind;
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
        $is_sph_available: expr,
        $is_cart_available: expr,
        $is_spinor_available: expr,
        $n_comp: expr,
        $n_spinor_comp: expr,
        $n_center: expr,
        $ng: expr,
        $integrator_category: literal,
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
            ) -> c_int {
                unsafe {
                    $integral_sph(out, dims, shls, atm, natm, bas, nbas, env, opt as _, cache)
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
            ) -> c_int {
                unsafe {
                    $integral_cart(out, dims, shls, atm, natm, bas, nbas, env, opt as _, cache)
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
            ) -> c_int {
                unsafe {
                    $integral_spinor(
                        out as _, dims, shls, atm, natm, bas, nbas, env, opt as _, cache,
                    )
                }
            }

            fn is_sph_available(&self) -> bool {
                $is_sph_available
            }

            fn is_cart_available(&self) -> bool {
                $is_cart_available
            }

            fn is_spinor_available(&self) -> bool {
                $is_spinor_available
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

            fn integrator_category(&self) -> &'static str {
                $integrator_category
            }

            fn name(&self) -> &'static str {
                $name
            }

            fn kind(&self) -> CIntKind {
                $kind
            }

            fn as_any(&self) -> &dyn Any {
                self
            }
        }
    };
}

/* #endregion */

#[allow(dead_code)]
pub(crate) unsafe fn panic_cart(
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
) -> c_int {
    panic!("Integral for cart is not implemented for this integrator.")
}

#[allow(dead_code)]
pub(crate) unsafe fn panic_spinor(
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
) -> c_int {
    panic!("Integral for spinor is not implemented for this integrator.")
}
