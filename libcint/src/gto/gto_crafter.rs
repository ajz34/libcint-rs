use crate::prelude::*;

/// GTO evaluation arguments for [`CInt::eval_gto_with_args`] and
/// [`CInt::eval_gto_with_args_spinor`].
///
/// Builder of this struct can be retrieved by [`CInt::gto_args_builder`] and
/// [`CInt::gto_args_builder_spinor`].
#[derive(Builder, Debug)]
#[builder(pattern = "owned", build_fn(error = "CIntError"))]
pub struct GtoArgs<'l, F> {
    /// Evaluator name, such as `"GTOval_sph_deriv1"`, `"GTOval_cart_ipig"`,
    /// etc.
    #[builder(setter(into))]
    pub eval_name: &'l str,

    /// Coordinates of grids, shape `(ngrid, 3)`.
    ///
    /// If you have `&[f64]` array instead of `&[[f64; 3]]`, you can use
    /// `core::slice::as_chunks` to convert it (rustc 1.88.0+).
    pub coord: &'l [[f64; 3]],

    /// Slice of shells to evaluate, such as `[0, 5]` to evaluate shells 0 to 4.
    ///
    /// Default is `None`, which means all shells will be evaluated.
    #[builder(default, setter(into))]
    pub shls_slice: Option<[usize; 2]>,

    /// Non-zero table for screening, shape `(nbas, nblk)`, where
    /// $n_\mathrm{blk} = \lceil \frac{n_\mathrm{grid}}{\texttt{BLKSIZE}}
    /// \rceil$.
    ///
    /// Default is `None`, and will be determined by value of struct field
    /// `cutoff`.
    #[builder(default, setter(strip_option))]
    pub non0tab: Option<&'l [u8]>,

    /// Cutoff for screening and generating `non0tab`.
    ///
    /// Default to [`CUTOFF`] or 1e-22.
    #[builder(default = Some(CUTOFF), setter(into))]
    pub cutoff: Option<f64>,

    /// Number of bins for screening.
    ///
    /// Default to [`NBINS`] or 100.
    #[builder(default, setter(into))]
    pub nbins: Option<u8>,

    /// Whether to fill the output array with zero during GTO value evaluation.
    ///
    /// Default to `true`.
    ///
    /// Please note that if the user does not provide output array, then we will
    /// always use something similar to [`alloc_zeroed`] (rust's internal trait
    /// function) to fast allocate a zeroed buffer, and will not fill zero
    /// during GTO value evaluation.
    ///
    /// [`alloc_zeroed`]: core::alloc::GlobalAlloc::alloc_zeroed
    #[builder(default = true)]
    pub fill_zero: bool,

    /// Scaling factor for GTO values.
    ///
    /// Default to 1.0.
    #[builder(default = 1.0)]
    pub fac: f64,

    /// Output array for GTO values, shape `(ngrid, nao, ncomp)` in column-major
    /// order.
    ///
    /// If not provided, the output array will be allocated internally. The
    /// allocation should call something similar to [`alloc_zeroed`] (rust's
    /// internal trait function).
    ///
    /// [`alloc_zeroed`]: core::alloc::GlobalAlloc::alloc_zeroed
    #[builder(default, setter(strip_option))]
    pub out: Option<&'l mut [F]>,
}

/// Get GTO evaluator by name.
///
/// # Returns
///
/// - [`Box<dyn GtoEvalAPI>`](GtoEvalAPI): The GTO evaluator instance.
/// - [`Option<CIntType>`](CIntType): The CIntType (Spheric/Cartesian/Spinor) if
///   specified in the name, None if not specified.
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
    /// Get GTO evaluator by name.
    ///
    /// This function is fallible.
    ///
    /// # See also
    ///
    /// - [`get_gto_eval_name`](CInt::get_gto_eval_name) for non-fallible
    ///   counterpart.
    pub fn get_gto_eval_name_f(&self, eval_name: &str) -> Result<(Box<dyn GtoEvalAPI>, Option<CIntType>), CIntError> {
        get_gto_eval_name_f(eval_name).ok_or(cint_error!(IntegratorNotFound, "{eval_name}"))
    }

    /// Get GTO evaluator by name.
    ///
    /// # Returns
    ///
    /// - [`Box<dyn GtoEvalAPI>`](GtoEvalAPI): The GTO evaluator instance.
    /// - [`Option<CIntType>`](CIntType): The CIntType
    ///   (Spheric/Cartesian/Spinor) if specified in the name, None if not
    ///   specified.
    ///
    /// # See also
    ///
    /// - [`get_gto_eval_name_f`](CInt::get_gto_eval_name_f) for fallible
    ///   counterpart.
    pub fn get_gto_eval_name(&self, eval_name: &str) -> (Box<dyn GtoEvalAPI>, Option<CIntType>) {
        self.get_gto_eval_name_f(eval_name).cint_unwrap()
    }

    /// Get GTO screening tabulation index (`non0tab`) on grids.
    ///
    /// This will generate a screening table with shape `(nbas, nblk)`, where
    /// $n_\mathrm{blk} = \lceil \frac{n_\mathrm{grid}}{\texttt{BLKSIZE}}
    /// \rceil$. [`BLKSIZE`] is a predefined constant 48.
    ///
    /// The screening table represents whether each (basis shell, grid block) is
    /// significant enough to be evaluated. If the tabulated value is zero, then
    /// the corresponding (basis shell, grid block) pair can be skipped during
    /// GTO evaluation, within the `cutoff` threadhold that the API user
    /// provided.
    ///
    /// # PySCF equivalent
    ///
    /// `pyscf.gto.eval_gto.make_screen_index`; please note that `blksize`
    /// argument is not configurable for rust's version, which is defined as a
    /// const [`BLKSIZE`].
    ///
    /// # Arguments
    ///
    /// - `coords`: Coordinates of grids, shape `(ngrid, 3)`.
    /// - `shls_slice`: Slice of shells to evaluate. If None, evaluate all
    ///   shells in the molecule.
    /// - `nbins`: Number of bins for screening. Default to [`NBINS`] or 100.
    /// - `cutoff`: Cutoff (GTO exponent without derivatives) for screening.
    ///   Default to [`CUTOFF`] or 1e-22.
    pub fn gto_screen_index(
        &self,
        coords: &[[f64; 3]],
        shls_slice: Option<[usize; 2]>,
        nbins: Option<u8>,
        cutoff: Option<f64>,
    ) -> CIntOutput<u8> {
        let shls_slice = shls_slice.unwrap_or([0, self.nbas()]);
        let (out, shape) = gto_screen_index(coords, shls_slice, nbins, cutoff, &self.atm, &self.bas, &self.env);
        CIntOutput { out: Some(out), shape: shape.into() }
    }
}

impl CInt {
    /// Create a builder for GTO evaluation arguments for
    /// [`CInt::eval_gto_with_args`].
    ///
    /// # See also
    ///
    /// - [`CInt::eval_gto_with_args`] for usage of the built arguments.
    /// - [`GtoArgsBuilder`] for the builder struct.
    /// - [`CInt::gto_args_builder_spinor`] for spinor version.
    pub fn gto_args_builder(&self) -> GtoArgsBuilder<'_, f64> {
        GtoArgsBuilder::default()
    }

    /// Evaluate GTO values on grids with given arguments.
    ///
    /// # See also
    ///
    /// - [`eval_gto_with_args`](CInt::eval_gto_with_args) for non-fallible
    ///   counterpart.
    pub fn eval_gto_with_args_f(&self, args: GtoArgs<f64>) -> Result<CIntOutput<f64>, CIntError> {
        let mut data = self.clone();
        let GtoArgs { eval_name, coord, shls_slice, non0tab, cutoff, nbins, mut fill_zero, fac, out } = args;

        // check name, sph/cart, set evaluator
        let (mut evaluator, cint_type) = data.get_gto_eval_name_f(eval_name)?;
        if let Some(cint_type) = cint_type {
            if cint_type == Spinor {
                cint_raise!(InvalidValue, "Spinor type is not supported for `eval_gto` or `eval_gto_with_args`. Use `eval_gto_spinor` or `eval_gto_with_args_spinor` instead.")?
            }
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

    /// Evaluate GTO values on grids with given arguments.
    ///
    /// For simple usage and more documentation, please refer to associated
    /// function [`eval_gto`](CInt::eval_gto).
    ///
    /// # PySCF equivalent
    ///
    /// `mol.eval_ao` with full arguments. Please note that `comp` and `ao_loc`
    /// arguments in PySCF is not covered in rust's version, since they are
    /// probably redundant.
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    ///
    /// // generate a grid of 2048 points
    /// let ngrid = 2048;
    /// let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    /// let args = cint_data.gto_args_builder()
    ///     .eval_name("deriv1_sph")
    ///     .coord(&coord)
    ///     .cutoff(1e-15)    
    ///     .build()
    ///     .unwrap();
    /// let (out, shape) = cint_data.eval_gto_with_args(args).into();
    /// assert_eq!(shape, vec![2048, 43, 4]); // (ngrid, nao, ncomp)
    /// ```
    ///
    /// # See also
    ///
    /// - [`eval_gto_with_args_f`](CInt::eval_gto_with_args_f) for fallible
    ///   counterpart.
    /// - [`GtoArgs`] for the arguments struct.
    /// - [`gto_args_builder`](CInt::gto_args_builder) for the builder of
    ///   arguments.
    /// - [`eval_gto`](CInt::eval_gto) for simpler usage.
    pub fn eval_gto_with_args(&self, args: GtoArgs<f64>) -> CIntOutput<f64> {
        self.eval_gto_with_args_f(args).cint_unwrap()
    }

    /// Evaluate GTO values on grids.
    ///
    /// This function is fallible.
    ///
    /// # See also
    ///
    /// - [`eval_gto`](CInt::eval_gto) for non-fallible counterpart.
    pub fn eval_gto_f(&self, eval_name: &str, coord: &[[f64; 3]]) -> Result<CIntOutput<f64>, CIntError> {
        self.eval_gto_with_args_f(GtoArgsBuilder::default().eval_name(eval_name).coord(coord).build()?)
    }

    /// Evaluate GTO values on grids.
    ///
    /// <div class="warning">
    ///
    /// **This function follows column-major convention.**
    ///
    /// For GTO value evaluation on grids, we do not provide another function
    /// for row-major convention, different to the case of [`integrate`] and
    /// [`integrate_row_major`] for electronic integrals.
    ///
    /// We always put the grids to be the most contiguous dimension, then the
    /// atomic orbitals (basis), then the components. For spinors, the spin
    /// dimension is the least contiguous dimension.
    ///
    /// We recommend the users who wish to have row-major output, only revert
    /// the list of `shape` (struct field of [`CIntOutput`]) after getting the
    /// output. We do not recommend API user of this function to perform
    /// explicit memory layout transposition; this will slow down your
    /// program: both in unnecessary memory bandwidth consumption, and cache
    /// efficiency reduction in DFT numerical integration.
    ///
    /// [`integrate`]: CInt::integrate
    /// [`integrate_row_major`]: CInt::integrate_row_major
    ///
    /// </div>
    ///
    /// # PySCF equivalent
    ///
    /// `mol.eval_ao(eval_name, coords)`
    ///
    /// # Arguments
    ///
    /// - `eval_name`: Evaluator name, such as `"GTOval_sph_deriv1"`.
    /// - `coord`: Coordinates of grids, shape `(ngrid, 3)`.
    ///
    /// # Notes on argument `eval_name`
    ///
    /// The `eval_name` can contain the following parts:
    /// - Component or operator: `"deriv0"`, `"deriv1"`, `"deriv2"`, `"deriv3"`,
    ///   `"ip"`, `"ig"`, etc.
    /// - Shell type: `"sph"` or `"cart"`. This affects the number of basis
    ///   functions per shell. If not specified, the type will be determined by
    ///   the current `cint_type` (struct field) of the [`CInt`] instance. This
    ///   function does not support `"spinor"` type; please use
    ///   [`eval_gto_spinor`] instead.
    /// - Prefix: `"GTO"` or `"GTOval"`. This is for PySCF naming convention
    ///   compatibility, not necessary.
    ///
    /// See also [`get_gto_eval_name_f`] for the actual implementation.
    ///
    /// [`eval_gto_spinor`]: CInt::eval_gto_spinor
    ///
    /// # Output
    ///
    /// Returns a [`CIntOutput<f64>`] for the evaluated GTO values, which
    /// contains
    /// - `out` (`Vec<f64>`): The output buffer for GTO values.
    /// - `shape` (`Vec<usize>`): The shape of output array in column-major
    ///   order `(ngrid, nao, ncomp)`.
    ///
    /// # Example
    ///
    /// The following example uses a pre-defined water molecule with def2-TZVP
    /// basis set.
    ///
    /// The construction of [`CInt`] is not performed in this crate,
    /// but should be provided by user. We refer to function
    /// [`init_h2o_def2_tzvp`] for how to define a [`CInt`] instance.
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    ///
    /// // generate a grid of 2048 points
    /// let ngrid = 2048;
    /// let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    /// let (out, shape) = cint_data.eval_gto("GTOval_sph_deriv1", &coord).into();
    /// assert_eq!(shape, vec![2048, 43, 4]); // (ngrid, nao, ncomp)
    /// ```
    ///
    /// # Supported evaluators
    ///
    /// | Component     | $n_\mathrm{comp}$ | Description |
    /// |---------------|-------------------|-------------|
    /// | `deriv0`      | 1                | GTO values |
    /// | `deriv1`      | 4                | GTO values and derivatives to 1st order (useful for GGA/meta-GGA Fock matrix formulation) |
    /// | `deriv2`      | 10               | GTO values and derivatives to 2nd order |
    /// | `deriv3`      | 20               | GTO values and derivatives to 3rd order |
    /// | `deriv4`      | 35               | GTO values and derivatives to 4th order |
    /// | `ip`          | 3                | derivative operator $\nabla$, or equivalently momemtum operator $\hat p$ scaled by imaginary unit $i$ |
    /// | `ig`          | 3                | GIAO operator at 1st order $\hat U_\mathrm{g}$ |
    /// | `ipig`        | 9                | Combined operator of `ip` and `ig` |
    /// | `ipr`         | 3                | derivative operator $\nabla$ multiplied by electronic position (taking atomic center as origin) |
    /// | `iprc`        | 3                | derivative operator $\nabla$ multiplied by electronic position (use certain common origin specified by [`set_common_origin`]) |
    /// | `sp`          | 1 × 4            | spinor sigma operator |
    /// | `ipsp`        | 3 × 4            | combined operator of `ip` and `sp` |
    /// | `ipipsp`      | 9 × 4            | combined operator of `ip`, `ip`, and `sp` |
    ///
    /// Note that for `sp` (sigma) related operators, for example `ipsp`, the
    /// component for non-spinor GTO is 3 × 4 or 12, but will be 3 for spinor
    /// GTO.
    ///
    /// [`set_common_origin`]: CInt::set_common_origin
    ///
    /// <div class="warning">
    ///
    /// **This function cannot handle spinor GTO evaluation**
    ///
    /// Due to rust's strict type system, it is very difficult to handle both
    /// spinor and spheric/cartesian GTO evaluation in the same function. So
    /// spinor GTO evaluation should be called by
    /// [`eval_gto_spinor`](Self::eval_gto_spinor), which output is
    /// [`Complex<f64>`].
    ///
    /// ```should_panic
    /// # use libcint::prelude::*;
    /// # let cint_data = init_h2o_def2_tzvp();
    /// # let coord: Vec<[f64; 3]> = (0..2048).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    /// let (out, shape) = cint_data.eval_gto("GTOval_spinor_deriv2", &coord).into();
    /// // InvalidValue: Spinor type is not supported
    /// ```
    ///
    /// </div>
    ///
    /// # See also
    ///
    /// - [`eval_gto_f`](CInt::eval_gto_f) for fallible counterpart.
    /// - [`eval_gto_spinor`](CInt::eval_gto_spinor) for spinor GTO evaluation.
    /// - [`eval_gto_with_args`](CInt::eval_gto_with_args) for more advanced
    ///   usage (full arguments that this crate supports).
    pub fn eval_gto(&self, eval_name: &str, coord: &[[f64; 3]]) -> CIntOutput<f64> {
        self.eval_gto_f(eval_name, coord).cint_unwrap()
    }
}

impl CInt {
    /// Create a builder for GTO evaluation arguments for spinor GTO evaluation.
    ///
    /// # See also
    ///
    /// - [`CInt::eval_gto_with_args_spinor`] for usage of the built arguments.
    /// - [`GtoArgsBuilder`] for the builder struct.
    /// - [`CInt::gto_args_builder`] for non-spinor version.
    pub fn gto_args_builder_spinor(&self) -> GtoArgsBuilder<'_, Complex<f64>> {
        GtoArgsBuilder::default()
    }

    /// Evaluate spinor GTO values on grids with given arguments.
    ///
    /// # See also
    ///
    /// - [`eval_gto_with_args_spinor`](CInt::eval_gto_with_args_spinor) for
    ///   non-fallible counterpart.
    pub fn eval_gto_with_args_spinor_f(&self, args: GtoArgs<Complex<f64>>) -> Result<CIntOutput<Complex<f64>>, CIntError> {
        let mut data = self.clone();
        let GtoArgs { eval_name, coord, shls_slice, non0tab, cutoff, nbins, mut fill_zero, fac, out } = args;

        // check name, sph/cart, set evaluator
        let (mut evaluator, cint_type) = data.get_gto_eval_name_f(eval_name)?;
        if cint_type.is_some_and(|cint_type| cint_type != Spinor) {
            cint_raise!(InvalidValue, "Only Spinor type is supported for `eval_gto_spinor` or `eval_gto_with_args_spinor`. For other types, use `eval_gto` or `eval_gto_with_args`.")?
        }
        data.set_cint_type(Spinor);
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
        let ntensor = evaluator.ntensor();
        let nao = ao_loc[sh1] - ao_loc[sh0];
        let out_size = ngrids * nao * ntensor * 2;
        let out_shape = vec![ngrids, nao, ntensor, 2];

        let mut out_vec = match out {
            Some(_) => None,
            None => {
                // rust is fast on zero-value allocation, so we always allocate zeroed vec here
                // also see `alloc::alloc::__rust_alloc_zeroed` or C `calloc`
                fill_zero = false;
                // for complex, note that `Complex::<f64>::zero()` may not be optimized by
                // rust's zeroed allocation, so need a transmute here
                let vec = vec![[0.0_f64, 0.0_f64]; out_size];
                Some(unsafe { transmute::<Vec<[f64; 2]>, Vec<Complex<f64>>>(vec) })
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
        gto_eval_spinor_loop(&data, evaluator.as_ref(), out, coord, fac, shls_slice, non0tab, fill_zero)?;
        Ok(CIntOutput { out: out_vec, shape: out_shape })
    }

    /// Evaluate spinor GTO values on grids with given arguments.
    ///
    /// # See also
    ///
    /// - [`eval_gto_with_args_spinor_f`](CInt::eval_gto_with_args_spinor_f) for
    ///   fallible counterpart.
    /// - [`GtoArgs`] for the arguments struct.
    /// - [`gto_args_builder_spinor`](CInt::gto_args_builder_spinor) for the
    ///   builder of arguments.
    /// - [`eval_gto_spinor`](CInt::eval_gto_spinor) for simpler usage.
    /// - [`eval_gto_with_args`](CInt::eval_gto_with_args) for non-spinor
    ///   version.
    pub fn eval_gto_with_args_spinor(&self, args: GtoArgs<Complex<f64>>) -> CIntOutput<Complex<f64>> {
        self.eval_gto_with_args_spinor_f(args).cint_unwrap()
    }

    /// Evaluate spinor GTO values on grids.
    ///
    /// This function is fallible.
    ///
    /// # See also
    ///
    /// - [`eval_gto_spinor`](CInt::eval_gto_spinor) for non-fallible
    ///   counterpart.
    pub fn eval_gto_spinor_f(&self, eval_name: &str, coord: &[[f64; 3]]) -> Result<CIntOutput<Complex<f64>>, CIntError> {
        self.eval_gto_with_args_spinor_f(GtoArgsBuilder::default().eval_name(eval_name).coord(coord).build()?)
    }

    /// Evaluate spinor GTO values on grids.
    ///
    /// For non-spinor GTO evaluation, please use [`eval_gto`](CInt::eval_gto).
    ///
    /// This function will return complex-valued GTO values, with shape `(ngrid,
    /// nao, ncomp, 2)` in column-major order, where the last dimension of
    /// size 2 represents the two spinor components (spin-up and spin-down).
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// let ngrid = 2048;
    /// let coord: Vec<[f64; 3]> = (0..ngrid).map(|i| [(i as f64).sin(), (i as f64).cos(), (i as f64 + 0.5).sin()]).collect();
    /// let (out, shape) = cint_data.eval_gto_spinor("GTOval_spinor_deriv1", &coord).into();
    /// assert_eq!(shape, vec![2048, 86, 4, 2]); // (ngrid, nao, ncomp, 2)
    /// let (out, shape) = cint_data.eval_gto_spinor("GTOval_spinor_ipsp", &coord).into();
    /// assert_eq!(shape, vec![2048, 86, 3, 2]); // (ngrid, nao, ncomp, 2)
    /// ```
    ///
    /// # See also
    ///
    /// - [`eval_gto_spinor_f`](CInt::eval_gto_spinor_f) for fallible
    ///   counterpart.
    /// - [`eval_gto`](CInt::eval_gto) for non-spinor GTO evaluation.
    /// - [`eval_gto_with_args_spinor`](CInt::eval_gto_with_args_spinor) for
    ///   more advanced usage (full arguments that this crate supports).
    pub fn eval_gto_spinor(&self, eval_name: &str, coord: &[[f64; 3]]) -> CIntOutput<Complex<f64>> {
        self.eval_gto_spinor_f(eval_name, coord).cint_unwrap()
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
