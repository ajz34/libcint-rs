use crate::prelude::*;

#[derive(Builder, Debug)]
#[builder(pattern = "owned", build_fn(error = "CIntError"))]
pub struct GtoArgs<'l, F: num::Num> {
    #[builder(setter(into))]
    pub eval_name: &'l str,

    pub coord: &'l [[f64; 3]],

    #[builder(default, setter(into))]
    pub shls_slice: Option<[usize; 2]>,

    #[builder(default, setter(strip_option))]
    pub non0tab: Option<&'l [u8]>,

    #[builder(default = Some(CUTOFF), setter(into))]
    pub cutoff: Option<f64>,

    #[builder(default, setter(into))]
    pub nbins: Option<u8>,

    #[builder(default = true)]
    pub fill_zero: bool,

    #[builder(default = F::one())]
    pub fac: F,

    #[builder(default, setter(strip_option))]
    pub out: Option<&'l mut [F]>,
}

pub fn get_gto_eval_name_f(eval_name: &str) -> Option<(Box<dyn GtoEvalAPI>, Option<CIntType>)> {
    let name = eval_name.to_lowercase().replace(['_', '-'], " ");
    let mut name_list = name.split_whitespace().collect_vec();

    // drop "gto" prefix
    name_list.retain(|&x| x != "gto" && x != "gtoval");

    // determine cint_type (sph/cart/spinor)
    let mut cint_type: Option<CIntType> = None;
    if let Some(pos) = name_list.iter().position(|&x| x == "sph") {
        name_list.remove(pos);
        cint_type = Some(Spheric);
    } else if let Some(pos) = name_list.iter().position(|&x| x == "cart") {
        name_list.remove(pos);
        cint_type = Some(Cartesian);
    } else if let Some(pos) = name_list.iter().position(|&x| x == "spinor") {
        name_list.remove(pos);
        cint_type = Some(Spinor);
    }

    match name_list.len() {
        0 => Some((Box::new(GtoEvalDeriv0), cint_type)),
        1 => match name_list[0] {
            "deriv0" => Some((Box::new(GtoEvalDeriv0), cint_type)),
            "deriv1" => Some((Box::new(GtoEvalDeriv1), cint_type)),
            "deriv2" => Some((Box::new(GtoEvalDeriv2), cint_type)),
            "deriv3" => Some((Box::new(GtoEvalDeriv3), cint_type)),
            "deriv4" => Some((Box::new(GtoEvalDeriv4), cint_type)),
            "ip" => Some((Box::new(GtoEvalDerivIp), cint_type)),
            "ig" => Some((Box::new(GtoEvalDerivIg), cint_type)),
            "ipig" => Some((Box::new(GtoEvalDerivIpIg), cint_type)),
            "ipr" => Some((Box::new(GtoEvalDerivIpR), cint_type)),
            "iprc" => Some((Box::new(GtoEvalDerivIpRc::default()), cint_type)),
            "sp" => Some((Box::new(GtoEvalDerivSp), cint_type)),
            "ipsp" => Some((Box::new(GtoEvalDerivIpSp), cint_type)),
            "ipipsp" => Some((Box::new(GtoEvalDerivIpIpSp), cint_type)),
            _ => None,
        },
        _ => None,
    }
}

impl CInt {
    pub fn gto_args_builder(&self) -> GtoArgsBuilder<'_, f64> {
        GtoArgsBuilder::default()
    }

    pub fn get_gto_eval_name_f(&self, eval_name: &str) -> Result<(Box<dyn GtoEvalAPI>, Option<CIntType>), CIntError> {
        get_gto_eval_name_f(eval_name).ok_or(cint_error!(IntegratorNotFound, "{eval_name}"))
    }

    pub fn get_gto_eval_name(&self, eval_name: &str) -> (Box<dyn GtoEvalAPI>, Option<CIntType>) {
        self.get_gto_eval_name_f(eval_name).cint_unwrap()
    }

    pub fn eval_gto_with_args_f(&self, args: GtoArgs<f64>) -> Result<CIntOutput<f64>, CIntError> {
        let mut data = self.clone();
        let GtoArgs { eval_name, coord, shls_slice, non0tab, cutoff, nbins, mut fill_zero, fac, out } = args;

        // check name, sph/cart, set evaluator
        let (mut evaluator, cint_type) = data.get_gto_eval_name_f(eval_name)?;
        if let Some(cint_type) = cint_type {
            data.set_cint_type(cint_type);
        }
        evaluator.init(&data);

        // handle output and check its sanity

        let shls_slice = match shls_slice {
            Some(slice) => slice,
            None => [0, data.nbas()],
        };
        let [sh0, sh1] = shls_slice;
        if sh0 > sh1 || sh1 > data.nbas() {
            cint_raise!(InvalidValue, "shls_slice [{sh0}, {sh1}] is invalid for nbas = {}", data.nbas())?
        }

        let nbas = sh1 - sh0;
        let ao_loc = &data.ao_loc();
        let ngrids = coord.len();
        let ncomp = evaluator.ncomp();
        let nao = ao_loc[sh1] - ao_loc[sh0];
        let out_size = ngrids * nao * ncomp;
        let out_shape = vec![ngrids, nao, ncomp];

        let mut out_vec = match out {
            Some(_) => None,
            None => {
                // rust is fast on zero-value allocation, so we always allocate zeroed vec here
                // also see `alloc::alloc::__rust_alloc_zeroed` or C `calloc`
                fill_zero = false;
                Some(vec![0.0; out_size])
            },
        };
        let out = match out {
            Some(out) => out,
            None => out_vec.as_mut().unwrap(),
        };
        if out.len() < out_size {
            cint_raise!(InvalidValue, "Output vector size {} is smaller than required size {out_size} of shape {out_shape:?}", out.len())?;
        }

        // handle non0tab and check its sanity
        if let Some(tab) = non0tab {
            let nblk = ngrids.div_ceil(BLKSIZE);
            let expected_len = nblk * nbas;
            let tab_len = tab.len();
            if tab_len != expected_len {
                cint_raise!(InvalidValue, "non0tab length {tab_len} is smaller than required size {expected_len}")?
            }
        }
        let non0tab_with_cutoff = if non0tab.is_none() && cutoff.is_some() {
            Some(gto_screen_index(coord, shls_slice, nbins, cutoff, &data.atm, &data.bas, &data.env).0)
        } else {
            None
        };
        let non0tab = non0tab.or(non0tab_with_cutoff.as_deref());

        // go to low-level evaluation function
        gto_eval_loop(&data, evaluator.as_ref(), out, coord, fac, shls_slice, non0tab, fill_zero)?;
        Ok(CIntOutput { out: out_vec, shape: out_shape })
    }

    pub fn eval_gto_with_args(&self, args: GtoArgs<f64>) -> CIntOutput<f64> {
        self.eval_gto_with_args_f(args).cint_unwrap()
    }

    pub fn eval_gto_f(&self, eval_name: &str, coord: &[[f64; 3]]) -> Result<CIntOutput<f64>, CIntError> {
        self.eval_gto_with_args_f(GtoArgsBuilder::default().eval_name(eval_name).coord(coord).build()?)
    }

    pub fn eval_gto(&self, eval_name: &str, coord: &[[f64; 3]]) -> CIntOutput<f64> {
        self.eval_gto_f(eval_name, coord).cint_unwrap()
    }
}

#[test]
fn playground() {
    let mol = init_h2o_def2_tzvp();
    let ngrid = 2048;
    let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    let args = mol.gto_args_builder().eval_name("deriv0_sph").coord(&coord).build().unwrap();
    let output = mol.eval_gto_with_args(args);
    println!("output: {:?}", output.shape);
    println!("fp = {}", cint_fp(output.out.as_ref().unwrap()));
}
