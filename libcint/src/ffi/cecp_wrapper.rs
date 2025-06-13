use crate::cint::CIntKind;
use crate::ffi::cecp_ffi::*;
use crate::ffi::wrapper_traits::Integrator;
use crate::ffi::wrapper_traits::panic_spinor;
use crate::impl_integrator;
use core::any::Any;
use core::ffi::{c_int, c_void};

impl_integrator!(
    ECPscalar,
    ECPscalar_optimizer,
    ECPscalar_sph,
    ECPscalar_cart,
    panic_spinor,
    true,
    true,
    false,
    1,
    1,
    2,
    vec![0, 0, 0, 0, 0, 1, 0, 1],
    "ECP",
    "ECPscalar",
    CIntKind::Ecp
);
impl_integrator!(
    ECPso,
    ECPso_optimizer,
    ECPso_sph,
    ECPso_cart,
    ECPso_spinor,
    true,
    true,
    true,
    3,
    1,
    2,
    vec![0, 0, 0, 0, 0, 3, 0, 1],
    "ECP",
    "ECPso",
    CIntKind::Ecp
);
impl_integrator!(
    ECPscalar_ignuc,
    ECPscalar_ignuc_optimizer,
    ECPscalar_ignuc_sph,
    ECPscalar_ignuc_cart,
    panic_spinor,
    true,
    true,
    false,
    3,
    3,
    2,
    vec![1, 0, 0, 0, 1, 1, 0, 3],
    "ECP",
    "ECPscalar_ignuc",
    CIntKind::Ecp
);
impl_integrator!(
    ECPscalar_ipnuc,
    ECPscalar_ipnuc_optimizer,
    ECPscalar_ipnuc_sph,
    ECPscalar_ipnuc_cart,
    panic_spinor,
    true,
    true,
    false,
    3,
    3,
    2,
    vec![1, 0, 0, 0, 1, 1, 0, 3],
    "ECP",
    "ECPscalar_ipnuc",
    CIntKind::Ecp
);
impl_integrator!(
    ECPscalar_ipipnuc,
    ECPscalar_ipipnuc_optimizer,
    ECPscalar_ipipnuc_sph,
    ECPscalar_ipipnuc_cart,
    panic_spinor,
    true,
    true,
    false,
    9,
    9,
    2,
    vec![2, 0, 0, 0, 2, 1, 0, 9],
    "ECP",
    "ECPscalar_ipipnuc",
    CIntKind::Ecp
);
impl_integrator!(
    ECPscalar_ipnucip,
    ECPscalar_ipnucip_optimizer,
    ECPscalar_ipnucip_sph,
    ECPscalar_ipnucip_cart,
    panic_spinor,
    true,
    true,
    false,
    9,
    9,
    2,
    vec![1, 1, 0, 0, 2, 1, 0, 9],
    "ECP",
    "ECPscalar_ipnucip",
    CIntKind::Ecp
);
impl_integrator!(
    ECPscalar_iprinv,
    ECPscalar_iprinv_optimizer,
    ECPscalar_iprinv_sph,
    ECPscalar_iprinv_cart,
    panic_spinor,
    true,
    true,
    false,
    3,
    3,
    2,
    vec![1, 1, 0, 0, 2, 1, 0, 3],
    "ECP",
    "ECPscalar_iprinv",
    CIntKind::Ecp
);
impl_integrator!(
    ECPscalar_ipiprinv,
    ECPscalar_ipiprinv_optimizer,
    ECPscalar_ipiprinv_sph,
    ECPscalar_ipiprinv_cart,
    panic_spinor,
    true,
    true,
    false,
    9,
    9,
    2,
    vec![2, 0, 0, 0, 2, 1, 0, 9],
    "ECP",
    "ECPscalar_ipiprinv",
    CIntKind::Ecp
);
impl_integrator!(
    ECPscalar_iprinvip,
    ECPscalar_iprinvip_optimizer,
    ECPscalar_iprinvip_sph,
    ECPscalar_iprinvip_cart,
    panic_spinor,
    true,
    true,
    false,
    9,
    9,
    2,
    vec![2, 0, 0, 0, 2, 1, 0, 9],
    "ECP",
    "ECPscalar_iprinvip",
    CIntKind::Ecp
);

pub fn get_ecp_integrator(name: &str) -> Option<Box<dyn Integrator>> {
    match name {
        "ecpscalar" => Some(Box::new(ECPscalar)),
        "ecpso" => Some(Box::new(ECPso)),
        "ecpscalar_ignuc" => Some(Box::new(ECPscalar_ignuc)),
        "ecpscalar_ipnuc" => Some(Box::new(ECPscalar_ipnuc)),
        "ecpscalar_ipipnuc" => Some(Box::new(ECPscalar_ipipnuc)),
        "ecpscalar_ipnucip" => Some(Box::new(ECPscalar_ipnucip)),
        "ecpscalar_iprinv" => Some(Box::new(ECPscalar_iprinv)),
        "ecpscalar_ipiprinv" => Some(Box::new(ECPscalar_ipiprinv)),
        "ecpscalar_iprinvip" => Some(Box::new(ECPscalar_iprinvip)),
        _ => None,
    }
}
