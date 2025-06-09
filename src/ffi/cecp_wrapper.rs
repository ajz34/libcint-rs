use crate::ffi::cecp;
use crate::ffi::cecp::ECPOpt;
use crate::ffi::wrapper_traits::ECPIntegrator;
use crate::impl_ecpintegrator;
use core::any::Any;
use core::ffi::c_int;

impl_ecpintegrator!(
    ECPscalar,
    ECPscalar_optimizer,
    ECPscalar_sph,
    ECPscalar_cart,
    1,
    "ECP",
    "ECPIntegratorBase"
);
impl_ecpintegrator!(
    ECPscalar_ignuc,
    ECPscalar_ignuc_optimizer,
    ECPscalar_ignuc_sph,
    ECPscalar_ignuc_cart,
    3,
    "ECP",
    "ECPscalar_ignuc"
);
impl_ecpintegrator!(
    ECPscalar_ipnuc,
    ECPscalar_ipnuc_optimizer,
    ECPscalar_ipnuc_sph,
    ECPscalar_ipnuc_cart,
    3,
    "ECP",
    "ECPscalar_ipnuc"
);
impl_ecpintegrator!(
    ECPscalar_ipipnuc,
    ECPscalar_ipipnuc_optimizer,
    ECPscalar_ipipnuc_sph,
    ECPscalar_ipipnuc_cart,
    9,
    "ECP",
    "ECPscalar_ipipnuc"
);
impl_ecpintegrator!(
    ECPscalar_ipnucip,
    ECPscalar_ipnucip_optimizer,
    ECPscalar_ipnucip_sph,
    ECPscalar_ipnucip_cart,
    9,
    "ECP",
    "ECPscalar_ipnucip"
);
impl_ecpintegrator!(
    ECPscalar_iprinv,
    ECPscalar_iprinv_optimizer,
    ECPscalar_iprinv_sph,
    ECPscalar_iprinv_cart,
    3,
    "ECP",
    "ECPscalar_iprinv"
);
impl_ecpintegrator!(
    ECPscalar_ipiprinv,
    ECPscalar_ipiprinv_optimizer,
    ECPscalar_ipiprinv_sph,
    ECPscalar_ipiprinv_cart,
    9,
    "ECP",
    "ECPscalar_ipiprinv"
);
impl_ecpintegrator!(
    ECPscalar_iprinvip,
    ECPscalar_iprinvip_optimizer,
    ECPscalar_iprinvip_sph,
    ECPscalar_iprinvip_cart,
    9,
    "ECP",
    "ECPscalar_iprinvip"
);

pub fn get_ecp_integrator(name: &str) -> Option<Box<dyn ECPIntegrator>> {
    match name {
        "ECPscalar" => Some(Box::new(ECPscalar)),
        "ECPscalar_ignuc" => Some(Box::new(ECPscalar_ignuc)),
        "ECPscalar_ipnuc" => Some(Box::new(ECPscalar_ipnuc)),
        "ECPscalar_ipipnuc" => Some(Box::new(ECPscalar_ipipnuc)),
        "ECPscalar_ipnucip" => Some(Box::new(ECPscalar_ipnucip)),
        "ECPscalar_iprinv" => Some(Box::new(ECPscalar_iprinv)),
        "ECPscalar_ipiprinv" => Some(Box::new(ECPscalar_ipiprinv)),
        "ECPscalar_iprinvip" => Some(Box::new(ECPscalar_iprinvip)),
        _ => None,
    }
}
