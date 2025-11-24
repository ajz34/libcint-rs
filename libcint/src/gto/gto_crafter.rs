use crate::prelude::*;

#[derive(Builder, Debug)]
#[builder(pattern = "owned", build_fn(error = "CIntError"))]
pub struct GTOArgs<'l, F> {
    #[builder(setter(into))]
    pub eval_name: &'l str,

    pub coord: &'l [[f64; 3]],

    #[builder(default = None)]
    pub shls_slice: Option<[usize; 2]>,

    #[builder(default = None)]
    pub non0tab: Option<&'l [u8]>,

    #[builder(default = false)]
    pub fill_zero: bool,

    #[builder(default = None)]
    pub out: Option<&'l mut [F]>,
}

pub fn get_gto_eval_name_f(eval_name: &str) -> Option<Box<dyn GTOEvalAPI>> {
    let name = eval_name.to_lowercase().trim_start_matches("gtoval_").trim_start_matches("gto_").to_string();
    match name.as_str() {
        "" | "deriv0" => Some(Box::new(GTOEvalDeriv0)),
        "deriv1" => Some(Box::new(GTOEvalDeriv1)),
        _ => None,
    }
}

impl CInt {
    pub fn gto_args_builder(&self) -> GTOArgsBuilder<'_, f64> {
        GTOArgsBuilder::default()
    }

    pub fn get_gto_eval_name_f(&self, eval_name: &str) -> Result<Box<dyn GTOEvalAPI>, CIntError> {
        get_gto_eval_name_f(eval_name).ok_or(cint_error!(IntegratorNotFound, "{eval_name}"))
    }

    pub fn get_gto_eval_name(&self, eval_name: &str) -> Box<dyn GTOEvalAPI> {
        self.get_gto_eval_name_f(eval_name).unwrap()
    }
}

#[test]
fn playground() {
    let mol = init_h2o_def2_tzvp();
    let args = mol.gto_args_builder().build();
    println!("{:?}", args);
}
