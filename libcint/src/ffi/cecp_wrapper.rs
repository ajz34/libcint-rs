use crate::cint_data::CintKind;
use crate::error::CintError;
use crate::ffi::cecp::*;
use crate::ffi::wrapper_traits::Integrator;
use crate::impl_integrator_sph_cart_only;
use core::any::Any;
use core::ffi::{c_int, c_void};

impl_integrator_sph_cart_only!(
    ECPscalar,
    ECPscalar_optimizer,
    ECPscalar_sph,
    ECPscalar_cart,
    1,
    1,
    2,
    vec![0, 0, 0, 0, 0, 1, 1, 1],
    "ECP",
    "ECPIntegratorBase",
    CintKind::Ecp
);
impl_integrator_sph_cart_only!(
    ECPscalar_ignuc,
    ECPscalar_ignuc_optimizer,
    ECPscalar_ignuc_sph,
    ECPscalar_ignuc_cart,
    3,
    3,
    2,
    vec![1, 0, 0, 0, 1, 1, 0, 3],
    "ECP",
    "ECPscalar_ignuc",
    CintKind::Ecp
);
impl_integrator_sph_cart_only!(
    ECPscalar_ipnuc,
    ECPscalar_ipnuc_optimizer,
    ECPscalar_ipnuc_sph,
    ECPscalar_ipnuc_cart,
    3,
    3,
    2,
    vec![1, 0, 0, 0, 1, 1, 0, 3],
    "ECP",
    "ECPscalar_ipnuc",
    CintKind::Ecp
);
impl_integrator_sph_cart_only!(
    ECPscalar_ipipnuc,
    ECPscalar_ipipnuc_optimizer,
    ECPscalar_ipipnuc_sph,
    ECPscalar_ipipnuc_cart,
    9,
    9,
    2,
    vec![2, 0, 0, 0, 2, 1, 0, 9],
    "ECP",
    "ECPscalar_ipipnuc",
    CintKind::Ecp
);
impl_integrator_sph_cart_only!(
    ECPscalar_ipnucip,
    ECPscalar_ipnucip_optimizer,
    ECPscalar_ipnucip_sph,
    ECPscalar_ipnucip_cart,
    9,
    9,
    2,
    vec![1, 1, 0, 0, 2, 1, 0, 9],
    "ECP",
    "ECPscalar_ipnucip",
    CintKind::Ecp
);
impl_integrator_sph_cart_only!(
    ECPscalar_iprinv,
    ECPscalar_iprinv_optimizer,
    ECPscalar_iprinv_sph,
    ECPscalar_iprinv_cart,
    9,
    9,
    2,
    vec![1, 1, 0, 0, 2, 1, 0, 9],
    "ECP",
    "ECPscalar_iprinv",
    CintKind::Ecp
);
impl_integrator_sph_cart_only!(
    ECPscalar_ipiprinv,
    ECPscalar_ipiprinv_optimizer,
    ECPscalar_ipiprinv_sph,
    ECPscalar_ipiprinv_cart,
    9,
    9,
    2,
    vec![2, 0, 0, 0, 2, 1, 0, 9],
    "ECP",
    "ECPscalar_ipiprinv",
    CintKind::Ecp
);
impl_integrator_sph_cart_only!(
    ECPscalar_iprinvip,
    ECPscalar_iprinvip_optimizer,
    ECPscalar_iprinvip_sph,
    ECPscalar_iprinvip_cart,
    9,
    9,
    2,
    vec![2, 0, 0, 0, 2, 1, 0, 9],
    "ECP",
    "ECPscalar_iprinvip",
    CintKind::Ecp
);

pub fn get_ecp_integrator(name: &str) -> Option<Box<dyn Integrator>> {
    match name {
        "ecpscalar" => Some(Box::new(ECPscalar)),
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
