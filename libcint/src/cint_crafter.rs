//! Integral crafter for [`CInt`] instance (the main integral functions).
//!
//! Implementation in this module is mostly not intended for basic users.

#![allow(dead_code)]

use crate::cint::{IntegrateArgsBuilder, IntorCrossArgs};
use crate::prelude::*;
use rayon::prelude::*;

/// Implementation of integral at higher API level (for advanced user usage).
impl CInt {
    /// Obtain integrator optimizer.
    ///
    /// # PySCF equivalent
    ///
    /// `gto.moleintor.make_cintopt`
    ///
    /// # Panics
    ///
    /// - integrator not found
    /// - `env` field not properly initialized for ECP integrator (should call
    ///   `merge_ecpbas` before this function)
    pub fn optimizer(&self, intor: &str) -> CIntOptimizer {
        self.optimizer_f(intor).unwrap()
    }

    /// Obtain integrator optimizer.
    ///
    /// This function is failable.
    ///
    /// # See also
    ///
    /// [`CInt::optimizer`]
    pub fn optimizer_f(&self, intor: &str) -> Result<CIntOptimizer, CIntError> {
        let integrator = crate::cint::get_integrator_f(intor)?;
        // make dummy shell slices and symms, only for checking
        let n_center = integrator.n_center();
        let cint_symm = CIntSymm::S1;
        let shls_slice = vec![[0, 0]; n_center];
        self.check_shls_slice(&*integrator, &shls_slice, cint_symm)?;
        Ok(self.get_optimizer(&*integrator))
    }

    /// Perform integral by arguments from builder.
    ///
    /// This is advanced integral function. Additionally to
    /// [`integrate`](CInt::integrate), this function also allows you to specify
    /// - `out`: output buffer for integral result, if not specified, a new
    ///   buffer will be allocated.
    /// - `row_major`: whether the output is in row-major order (default is
    ///   false, which means column-major order).
    ///
    /// This function is not very convenient to use, since rust does not allow
    /// named arguments.
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let data = init_h2o_def2_tzvp();
    /// let mut out_buffer = vec![0.0; 100000];
    /// let shls_slice = [[0, data.nbas()], [0, data.nbas()], [5, 10]];
    ///
    /// // col-major, the same memory layout to PySCF's 2/3-center integrals
    /// let args = data.integrate_args_builder()
    ///     .intor("int3c2e_ip2")
    ///     .aosym("s2ij")
    ///     .shls_slice(&shls_slice)
    ///     .out(&mut out_buffer)
    ///     .row_major(false)
    ///     .build()
    ///     .unwrap();
    /// let output = data.integrate_with_args(args);
    /// assert!(output.out.is_none());
    /// assert_eq!(output.shape, [946, 19, 3]);
    /// assert!((cint_fp(&out_buffer) - -15.467133235509742).abs() < 1e-10);
    ///
    /// // row-major, the same shape to PySCF's integrals
    /// let args = data.integrate_args_builder()
    ///     .intor("int3c2e_ip2")
    ///     .aosym("s2ij")
    ///     .shls_slice(&shls_slice)
    ///     .out(&mut out_buffer)
    ///     .row_major(true)
    ///     .build()
    ///     .unwrap();
    /// let output = data.integrate_with_args(args);
    /// assert!(output.out.is_none());
    /// assert_eq!(output.shape, [3, 946, 19]);
    /// assert!((cint_fp(&out_buffer) - 16.69819847729145).abs() < 1e-10);
    /// ```
    ///
    /// PySCF equivalent (`out_buffer` column-major, `out` row-major for
    /// 3-center integrals in PySCF):
    ///
    /// ```python
    /// mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="def2-TZVP").build()
    /// out_buffer = np.zeros(100000)
    /// shls_slice = [0, mol.nbas, 0, mol.nbas, 5, 10]
    /// out = mol.intor("int3c2e_ip2", aosym="s2ij", shls_slice=shls_slice, out=out_buffer)
    /// assert abs(lib.fp(out_buffer) - -15.467133235509742) < 1e-10
    /// assert abs(lib.fp(out) - 16.69819847729145) < 1e-10
    /// ```
    ///
    /// # See also
    ///
    /// - [`CInt::integrate`], [`CInt::integrate_row_major`]: for basic integral
    ///   evaluation.
    /// - [`CInt::integrate_with_args_spinor`]: for spinor integrals.
    /// - [`CInt::integrate_cross_with_args`]: for cross integrals with multiple
    ///   `CInt` instances.
    /// - [`CInt::integrate_with_args_f`]: for failable version of this
    ///   function.
    pub fn integrate_with_args(&self, args: IntegrateArgs<f64>) -> CIntOutput<f64> {
        self.integrate_with_args_inner(args).unwrap()
    }

    /// Perform integral by arguments from builder, failable version.
    ///
    /// # See also
    ///
    /// [`CInt::integrate_with_args`]
    pub fn integrate_with_args_f(&self, args: IntegrateArgs<f64>) -> Result<CIntOutput<f64>, CIntError> {
        self.integrate_with_args_inner(args)
    }

    /// Perform integral by arguments from builder for spinor integrals.
    ///
    /// # See also
    ///
    /// [`CInt::integrate_with_args`]
    /// [`CInt::integrate_with_args_spinor_f`]
    pub fn integrate_with_args_spinor(&self, args: IntegrateArgs<Complex<f64>>) -> CIntOutput<Complex<f64>> {
        self.integrate_with_args_inner(args).unwrap()
    }

    /// Perform integral by arguments from builder for spinor integrals,
    /// failable version.
    ///
    /// # See also
    ///
    /// [`CInt::integrate_with_args_spinor`]
    pub fn integrate_with_args_spinor_f(&self, args: IntegrateArgs<Complex<f64>>) -> Result<CIntOutput<Complex<f64>>, CIntError> {
        self.integrate_with_args_inner(args)
    }

    pub fn integrate_with_args_inner<F>(&self, args: IntegrateArgs<F>) -> Result<CIntOutput<F>, CIntError>
    where
        F: ComplexFloat + Send + Sync,
    {
        let mut data = if self.is_ecp_merged() { self.clone() } else { self.merge_ecpbas() };
        if let Some(grids) = args.grids {
            data.set_grids(grids);
        }

        // parse intor sph/cart/spinor
        let mut intor = args.intor;
        if intor.ends_with("_sph") {
            intor = &intor[..intor.len() - 4];
            data.cint_type = Spheric;
        } else if intor.ends_with("_cart") {
            intor = &intor[..intor.len() - 5];
            data.cint_type = Cartesian;
        } else if intor.ends_with("_spinor") {
            intor = &intor[..intor.len() - 7];
            data.cint_type = Spinor;
        }

        // parse grids
        let has_grids = intor.starts_with("int1e_grids");
        if has_grids && (data.ngrids() == 0) {
            cint_raise!(IntegratorNotAvailable, "`int1e_grids` requires grids to integrate.")?
        }

        // unwrap and prepare (without output checking)
        let integrator = CInt::get_integrator(intor);
        let mut shls_slice = args.shls_slice.iter().map(|&[shl0, shl1]| [shl0 as c_int, shl1 as c_int]).collect_vec();
        if shls_slice.is_empty() {
            shls_slice = vec![[0, data.nbas() as c_int]; integrator.n_center()];
        }
        let aosym = args.aosym;
        let cint_opt = data.get_optimizer(&*integrator);

        // perform most checks
        data.check_float_type::<F>()?;
        data.check_shls_slice(&*integrator, &shls_slice, aosym)?;
        data.check_optimizer(&*integrator, &cint_opt)?;

        // prepare output and check size of output
        let mut out_shape = data.cgto_shape(&shls_slice, aosym);
        if has_grids {
            let ngrids = data.ngrids();
            out_shape.insert(0, ngrids);
        }
        let n_comp = match data.cint_type {
            Spheric | Cartesian => integrator.n_comp(),
            Spinor => integrator.n_spinor_comp(),
        };
        if n_comp > 1 {
            match args.row_major {
                false => out_shape.push(n_comp),
                true => out_shape.insert(0, n_comp),
            }
        }
        let out_size = out_shape.iter().product();

        let mut out_vec = match args.out {
            Some(_) => None,
            None => Some(unsafe { aligned_uninitialized_vec::<F>(out_size) }),
        };
        let out = match args.out {
            Some(out) => out,
            None => out_vec.as_mut().unwrap(),
        };
        if out.len() < out_size {
            cint_raise!(InvalidValue, "Output vector size {} is smaller than required size {}", out.len(), out_size)?;
        }

        // actual integral execution
        match (has_grids, args.row_major) {
            (false, false) => {
                // col-major, output is [shl0, shl1, ..., comp]
                data.integral_inplace(&*integrator, out, &shls_slice, Some(&cint_opt), aosym)?;
            },
            (false, true) => {
                // row-major, output is [comp, shl0, shl1, ...]
                data.integral_row_major_inplace(&*integrator, out, &shls_slice, Some(&cint_opt), aosym)?;
            },
            (true, false) => {
                // col-major with grids, output is [grids, shl0, shl1, comp]
                data.integral_grids_inplace(&*integrator, out, &shls_slice, Some(&cint_opt))?;
            },
            (true, true) => {
                // row-major with grids, output is [comp, grids, shl0, shl1]
                data.integral_grids_row_major_inplace(&*integrator, out, &shls_slice, Some(&cint_opt))?;
            },
        };

        Ok(CIntOutput { out: out_vec, shape: out_shape })
    }

    /// Perform cross integral with multiple `CInt` instances.
    ///
    /// This function is similar to [`CInt::integrate_with_args`], but it
    /// allows you to specify multiple `CInt` instances, which is useful for
    /// integrals with auxiliary basis sets.
    ///
    /// Corresponding builder is [`CInt::integrate_cross_args_builder`].
    ///
    /// # See also
    ///
    /// - [`CInt::integrate_cross`] and [`CInt::integrate_cross_row_major`]: for
    ///   basic integral evaluation.
    /// - [`CInt::integrate_cross_with_args_f`]: for failable version of this
    ///   function.
    /// - [`CInt::integrate_cross_with_args_spinor`]: for spinor integrals.
    pub fn integrate_cross_with_args(args: IntorCrossArgs<f64>) -> CIntOutput<f64> {
        CInt::integrate_cross_with_args_inner(args).unwrap()
    }

    /// Perform cross integral with multiple `CInt` instances, failable version.
    ///
    /// # See also
    ///
    /// [`CInt::integrate_cross_with_args`]
    pub fn integrate_cross_with_args_f<F>(args: IntorCrossArgs<f64>) -> Result<CIntOutput<f64>, CIntError> {
        CInt::integrate_cross_with_args_inner(args)
    }

    /// Perform cross integral with multiple `CInt` instances for spinor
    /// integrals.
    ///
    /// # See also
    ///
    /// [`CInt::integrate_cross_with_args`]
    pub fn integrate_cross_with_args_spinor(args: IntorCrossArgs<Complex<f64>>) -> CIntOutput<Complex<f64>> {
        CInt::integrate_cross_with_args_inner(args).unwrap()
    }

    /// Perform cross integral with multiple `CInt` instances for spinor
    /// integrals, failable version.
    ///
    /// # See also
    ///
    /// [`CInt::integrate_cross_with_args_spinor`]
    pub fn integrate_cross_with_args_spinor_f<F>(args: IntorCrossArgs<Complex<f64>>) -> Result<CIntOutput<Complex<f64>>, CIntError> {
        CInt::integrate_cross_with_args_inner(args)
    }

    pub fn integrate_cross_with_args_inner<F>(args: IntorCrossArgs<F>) -> Result<CIntOutput<F>, CIntError>
    where
        F: ComplexFloat + Send + Sync,
    {
        let IntorCrossArgs { intor, mols, shls_slice, aosym, out, row_major } = args;

        let integrator = CInt::get_integrator_f(intor)?;
        let n_center = integrator.n_center();
        if mols.len() != n_center {
            cint_raise!(InvalidValue, "Number of molecules ({}) does not match number of centers ({n_center})", mols.len())?;
        }
        let mut shls_slice = if shls_slice.is_empty() {
            mols.iter().map(|mol| [0, mol.nbas()]).collect_vec()
        } else {
            if shls_slice.len() != n_center {
                cint_raise!(InvalidValue, "Number of shell slices ({}) does not match number of centers ({n_center})", shls_slice.len())?;
            }
            // check each slice
            for i in 0..n_center {
                let mol = &mols[i];
                let slc = shls_slice[i];

                if !(slc[0] <= slc[1] && slc[1] <= mol.nbas()) {
                    cint_raise!(InvalidValue, "Shell slice {slc:?} at index {i} is not proper within [0, {}]", mol.nbas())?;
                }
            }
            shls_slice.to_vec()
        };

        // vec of [mol_to_concate, nbas_offset]
        let mut to_concate = vec![(&mols[0], 0)];
        let mut mol_concate = mols[0].clone();

        for i in 1..n_center {
            let mut has_equal = false;
            for &(&mol_to_concate, nbas_offset) in to_concate.iter() {
                if mols[i] == mol_to_concate {
                    has_equal = true;
                    // update nbas_offset
                    shls_slice[i][0] += nbas_offset;
                    shls_slice[i][1] += nbas_offset;
                    break;
                }
            }
            if !has_equal {
                // add new molecule to concatenate
                let nbas_offset = mol_concate.nbas();
                mol_concate = &mol_concate + mols[i];
                to_concate.push((&mols[i], nbas_offset));
                // update shls_slice
                shls_slice[i][0] += nbas_offset;
                shls_slice[i][1] += nbas_offset;
            }
        }

        let args = if let Some(out) = out {
            IntegrateArgsBuilder::default().intor(intor).aosym(aosym).shls_slice(&shls_slice).row_major(row_major).out(out).build()?
        } else {
            IntegrateArgsBuilder::default().intor(intor).aosym(aosym).shls_slice(&shls_slice).row_major(row_major).build()?
        };

        mol_concate.integrate_with_args_inner(args)
    }
}

/// Obtaining integrator and argument builder.
///
/// These functions are mostly not for basic users, but for internal use or
/// advanced users. For basic users, you can use [`CInt::integrate`] (and
/// [`CInt::integrate_spinor`] for spinor integrals) to evaluate integrals.
impl CInt {
    /// Obtain integrator by name.
    ///
    /// This function does not require `CInt` instance. Make it here only for
    /// convenience.
    ///
    /// # Panics
    ///
    /// - integrator not found
    pub fn get_integrator(intor: &str) -> Box<dyn Integrator> {
        get_integrator(intor)
    }

    /// Obtain integrator by name, failable version.
    ///
    /// # See also
    ///
    /// - [`CInt::get_integrator`]
    pub fn get_integrator_f(intor: &str) -> Result<Box<dyn Integrator>, CIntError> {
        get_integrator_f(intor)
    }

    /// Obtain the builder of arguments for advanced integral evaluation.
    ///
    /// # See also
    ///
    /// - [`CInt::integrate_with_args`]
    pub fn integrate_args_builder(&self) -> IntegrateArgsBuilder<'static, f64> {
        IntegrateArgsBuilder::default()
    }

    /// Obtain the builder of arguments for advanced integral evaluation
    /// for spinor integrals.
    ///
    /// # See also
    ///
    /// - [`CInt::integrate_with_args_spinor`]
    pub fn integrate_args_builder_spinor(&self) -> IntegrateArgsBuilder<'static, Complex<f64>> {
        IntegrateArgsBuilder::default()
    }

    /// Obtain the builder of arguments for cross integral evaluation.
    ///
    /// # See also
    ///
    /// - [`CInt::integrate_cross_with_args`]
    pub fn integrate_cross_args_builder(&self) -> IntorCrossArgsBuilder<'static, f64> {
        IntorCrossArgsBuilder::default()
    }

    /// Obtain the builder of arguments for cross integral evaluation
    ///
    /// # See also
    ///
    /// - [`CInt::integrate_cross_with_args_spinor`]
    pub fn integrate_cross_args_builder_spinor(&self) -> IntorCrossArgsBuilder<'static, Complex<f64>> {
        IntorCrossArgsBuilder::default()
    }
}

/// Implementation of integral at higher API level (legacy support).
impl CInt {
    /// Legacy support for integral evaluation for rest_libcint.
    ///
    /// This function can be replaced by [`CInt::integrate`].
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// use libcint::ffi::cint_wrapper::int1e_ipkin;
    /// let cint_data = init_h2o_def2_tzvp();
    ///
    /// let (out, shape) = cint_data.integral_s1::<int1e_ipkin>(None);
    /// assert_eq!(shape, [43, 43, 3]);
    /// ```
    pub fn integral_s1<T>(&self, shls_slice: Option<&[[c_int; 2]]>) -> (Vec<f64>, Vec<usize>)
    where
        T: Integrator + Default,
    {
        self.integral_inner::<T, f64>(shls_slice, CIntSymm::S1)
    }

    /// Legacy support for integral evaluation for rest_libcint with spinor
    /// integrals.
    ///
    /// This function can be replaced by [`CInt::integrate_spinor`].
    pub fn integral_s1_spinor<T>(&self, shls_slice: Option<&[[c_int; 2]]>) -> (Vec<Complex<f64>>, Vec<usize>)
    where
        T: Integrator + Default,
    {
        self.integral_inner::<T, Complex<f64>>(shls_slice, CIntSymm::S1)
    }

    /// Legacy support for integral evaluation for rest_libcint with s2ij
    /// symmetry.
    ///
    /// This function can be replaced by [`CInt::integrate`].
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// use libcint::ffi::cint_wrapper::int3c2e_ip2;
    /// let cint_data = init_h2o_def2_tzvp();
    ///
    /// let (out, shape) = cint_data.integral_s2ij::<int3c2e_ip2>(None);
    /// assert_eq!(shape, [946, 43, 3]);
    /// ```
    pub fn integral_s2ij<T>(&self, shls_slice: Option<&[[c_int; 2]]>) -> (Vec<f64>, Vec<usize>)
    where
        T: Integrator + Default,
    {
        self.integral_inner::<T, f64>(shls_slice, CIntSymm::S2ij)
    }

    /// Legacy support for integral evaluation for rest_libcint with s2ij
    /// symmetry and spinor integrals.
    ///
    /// This function can be replaced by [`CInt::integrate_spinor`].
    pub fn integral_s2ij_spinor<T>(&self, shls_slice: Option<&[[c_int; 2]]>) -> (Vec<Complex<f64>>, Vec<usize>)
    where
        T: Integrator + Default,
    {
        self.integral_inner::<T, Complex<f64>>(shls_slice, CIntSymm::S2ij)
    }

    pub fn integral_inner<T, F>(&self, shls_slice: Option<&[[c_int; 2]]>, aosym: CIntSymm) -> (Vec<F>, Vec<usize>)
    where
        T: Integrator + Default,
        F: ComplexFloat + Send + Sync,
    {
        // prepare and unwrap
        let data = if self.is_ecp_merged() { self } else { &self.merge_ecpbas() };
        let integrator = T::default();
        let shls_slice = match shls_slice {
            Some(shls_slice) => shls_slice,
            None => &vec![[0, data.nbas() as c_int]; integrator.n_center()],
        };

        // additional check
        data.check_float_type::<F>().unwrap();
        data.check_shls_slice(&integrator, shls_slice, aosym).unwrap();

        // integral preparation and execution
        let mut out_shape = data.cgto_shape(shls_slice, aosym);
        let n_comp = match self.cint_type {
            Spheric | Cartesian => integrator.n_comp(),
            Spinor => integrator.n_spinor_comp(),
        };
        if n_comp > 1 {
            out_shape.push(n_comp);
        }
        let out_size = out_shape.iter().product();
        let mut out = unsafe { aligned_uninitialized_vec::<F>(out_size) };
        let cint_opt = data.get_optimizer(&integrator);
        data.integral_inplace(&integrator, &mut out, shls_slice, Some(&cint_opt), aosym).unwrap();
        (out, out_shape)
    }
}

/// Merge ECP data for integral evaluation.
impl CInt {
    /// Creates a new `CIntData` instance, with ECP integral information
    /// written to `bas` field, and properly initializes values in `env` field.
    ///
    /// This function is not intended to be used directly by users. This
    /// function is often used internally for data preparation.
    pub fn merge_ecpbas(&self) -> CInt {
        const AS_ECPBAS_OFFSET: usize = cecp_ffi::AS_ECPBAS_OFFSET as usize;
        const AS_NECPBAS: usize = cecp_ffi::AS_NECPBAS as usize;

        if self.is_ecp_merged() {
            self.clone()
        } else {
            let mut merged = self.clone();
            merged.bas.extend_from_slice(&self.ecpbas);
            merged.env[AS_ECPBAS_OFFSET] = self.bas.len() as f64;
            merged.env[AS_NECPBAS] = self.ecpbas.len() as f64;
            merged
        }
    }

    /// Decouples ECP data from the `CInt` instance, returning a new `CInt`
    /// instance that have basis data and ecp data separated.
    ///
    /// This function is not intended to be used directly by users. This
    /// function is often used internally for data preparation.
    pub fn decopule_ecpbas(&self) -> CInt {
        const AS_ECPBAS_OFFSET: usize = cecp_ffi::AS_ECPBAS_OFFSET as usize;
        const AS_NECPBAS: usize = cecp_ffi::AS_NECPBAS as usize;

        if !self.is_ecp_merged() {
            self.clone()
        } else {
            let mut decoupled = self.clone();
            let nbas = self.env[AS_ECPBAS_OFFSET] as usize;
            let (bas, ecpbas) = decoupled.bas.split_at(nbas);
            let bas = bas.to_vec();
            let ecpbas = ecpbas.to_vec();
            decoupled.bas = bas;
            decoupled.ecpbas = ecpbas;
            decoupled.env[AS_ECPBAS_OFFSET] = 0.0;
            decoupled.env[AS_NECPBAS] = 0.0;
            decoupled
        }
    }

    /// Check whether the `CInt` instance has been merged with ECP data.
    pub fn is_ecp_merged(&self) -> bool {
        self.env[cecp_ffi::AS_ECPBAS_OFFSET as usize] != 0.0
    }
}

/// Implementation of integral at lower API level (integral preparation).
impl CInt {
    /* #region check_integratable */

    /// Check if the integral can be integrated with the given integrator.
    ///
    /// # Example
    ///
    /// Following examples always uses H2O/def2-tzvp basis set with spheric.
    ///
    /// ```no_run
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// ```
    ///
    /// Successful case:
    ///
    /// ```rust
    /// # use libcint::prelude::*;
    /// # let cint_data = init_h2o_def2_tzvp();
    /// let integrator = CInt::get_integrator("int2e");
    /// let shls_slice = &[[3, 13], [3, 13], [0, 19], [0, 19]];
    /// let check = cint_data.check_shls_slice(&*integrator, shls_slice, CIntSymm::S4);
    /// assert!(check.is_ok());
    /// ```
    ///
    /// Failed case (shls_slice not matching to cint_symm):
    ///
    /// ```rust
    /// # use libcint::prelude::*;
    /// # let cint_data = init_h2o_def2_tzvp();
    /// let integrator = CInt::get_integrator("int2e");
    /// let shls_slice = &[[0, 15], [3, 12], [0, 19], [0, 19]]; // first two slices different
    /// let check = cint_data.check_shls_slice(&*integrator, shls_slice, CIntSymm::S4);
    /// assert!(check.is_err());
    /// ```
    ///
    /// Failed case (number of shls_slice not matching to number of centers):
    ///
    /// ```rust
    /// # use libcint::prelude::*;
    /// # let cint_data = init_h2o_def2_tzvp();
    /// let integrator = CInt::get_integrator("int1e_ipkin");
    /// let shls_slice = &[[0, 15], [0, 15], [0, 19]]; // thr3e slices, but integrator requires 2
    /// let check = cint_data.check_shls_slice(&*integrator, shls_slice, CIntSymm::S2ij);
    /// assert!(check.is_err());
    /// ```
    pub fn check_shls_slice(&self, integrator: &dyn Integrator, shls_slice: &[[c_int; 2]], cint_symm: CIntSymm) -> Result<(), CIntError> {
        let n_center = integrator.n_center();

        // length of shls_slice must be the same value to number of centers of
        // integrator.
        if shls_slice.len() != n_center {
            cint_raise!(
                IntegratorNotAvailable,
                "Integrator {} requires {} centers, but got {}",
                integrator.name(),
                integrator.n_center(),
                shls_slice.len()
            )?;
        }

        // Some integrators may not support certain shell types.
        match self.cint_type {
            Spheric => {
                if !integrator.is_sph_available() {
                    cint_raise!(IntegratorNotAvailable, "Integrator {} does not support spherical integrals", integrator.name())?;
                }
            },
            Cartesian => {
                if !integrator.is_cart_available() {
                    cint_raise!(IntegratorNotAvailable, "Integrator {} does not support Cartesian integrals", integrator.name())?;
                }
            },
            Spinor => {
                if !integrator.is_spinor_available() {
                    cint_raise!(IntegratorNotAvailable, "Integrator {} does not support spinor integrals", integrator.name())?;
                }
            },
        }

        // Every shell in shls_slice must be within the range of number of shells.
        let nbas = self.nbas() as c_int;
        for (i, shl) in shls_slice.iter().enumerate() {
            if !(0 <= shl[0] && shl[0] <= shl[1] && shl[1] <= nbas) {
                cint_raise!(InvalidValue, "Shell slice {shl:?} at index {i} is not proper within [0, {nbas}]")?;
            }
        }

        // check symmetry
        match cint_symm {
            CIntSymm::S2ij => {
                if shls_slice[0] != shls_slice[1] {
                    cint_raise!(
                        InvalidValue,
                        "Integrator {} does not support S2ij symmetry with different shell slices ({:?}, {:?})",
                        integrator.name(),
                        shls_slice[0],
                        shls_slice[1]
                    )?;
                }
            },
            CIntSymm::S2kl => {
                if n_center != 4 {
                    cint_raise!(InvalidValue, "Integrator {} requires 4 centers for S2kl symmetry, but got {n_center}", integrator.name())?;
                }
                if shls_slice[2] != shls_slice[3] {
                    cint_raise!(
                        InvalidValue,
                        "Integrator {} does not support S2kl symmetry with different shell slices ({:?}, {:?})",
                        integrator.name(),
                        shls_slice[2],
                        shls_slice[3]
                    )?;
                }
            },
            CIntSymm::S4 => {
                if n_center != 4 {
                    cint_raise!(InvalidValue, "Integrator {} requires 4 centers for S4 symmetry, but got {n_center}", integrator.name())?;
                }
                if shls_slice[0] != shls_slice[1] || shls_slice[2] != shls_slice[3] {
                    cint_raise!(
                        InvalidValue,
                        "Integrator {} does not support S4 symmetry with different shell slices {shls_slice:?})",
                        integrator.name()
                    )?;
                }
            },
            CIntSymm::S8 => {
                if n_center != 4 {
                    cint_raise!(InvalidValue, "Integrator {} requires 4 centers for S8 symmetry, but got {n_center}", integrator.name())?;
                }
                if shls_slice[0] != shls_slice[1] || shls_slice[0] != shls_slice[2] || shls_slice[0] != shls_slice[3] {
                    cint_raise!(
                        InvalidValue,
                        "Integrator {} does not support S8 symmetry with different shell slices {shls_slice:?})",
                        integrator.name()
                    )?;
                }
            },
            CIntSymm::S1 => {
                // S1 symmetry is always available
            },
        }

        // grids only available for S1 symmetry
        if integrator.name().starts_with("int1e_grids") && cint_symm != CIntSymm::S1 {
            cint_raise!(InvalidValue, "Integrator {} does not support grids with symmetry other than s1", integrator.name())?;
        }

        Ok(())
    }

    /// Check if the float type is correct for the cint_type.
    ///
    /// - For `Spheric` and `Cartesian`, it should be `f64`.
    /// - For `Spinor`, it should be `Complex<f64>`.
    ///
    /// Note that we only check size of type, instead of its real type.
    pub fn check_float_type<F>(&self) -> Result<(), CIntError>
    where
        F: ComplexFloat,
    {
        let expected = match self.cint_type {
            Spheric | Cartesian => 8, // f64
            Spinor => 16,             // Complex<f64>
        };
        let actual = std::mem::size_of::<F>();
        if actual != expected {
            cint_raise!(
                InvalidValue,
                "Expected float type size {expected} bytes, but got {actual} bytes. This is probably because your `cint_type` does not match integral (`integral_spinor` for spinor, otherwise `integral` for spheric/cartesian)."
            )?;
        }
        Ok(())
    }

    /// Check if the optimizer is compatible with the integrator.
    ///
    /// This check only involves that ECP and general integrals have different
    /// optimizers.
    pub fn check_optimizer(&self, integrator: &dyn Integrator, optimizer: &CIntOptimizer) -> Result<(), CIntError> {
        match optimizer {
            CIntOptimizer::Int(_) => {
                if integrator.kind() != CIntKind::Int {
                    cint_raise!(InvalidValue, "Optimizer is for Int, but integrator is {:?}", integrator.kind())?;
                }
            },
            CIntOptimizer::Ecp(_) => {
                if integrator.kind() != CIntKind::Ecp {
                    cint_raise!(InvalidValue, "Optimizer is for Ecp, but integrator is {:?}", integrator.kind())?;
                }
            },
        }
        Ok(())
    }

    /* #endregion */

    /* #region cgto_shape, cgto_loc */

    pub fn cgto_shape(&self, shls_slice: &[[c_int; 2]], aosym: CIntSymm) -> Vec<usize> {
        match aosym {
            CIntSymm::S1 => self.cgto_shape_s1(shls_slice),
            CIntSymm::S2ij => self.cgto_shape_s2ij(shls_slice),
            CIntSymm::S2kl => self.cgto_shape_s2kl(shls_slice),
            CIntSymm::S4 => self.cgto_shape_s4(shls_slice),
            CIntSymm::S8 => self.cgto_shape_s8(shls_slice),
        }
    }

    /// Shape of integral (in atomic orbital basis, symmetry s1), for specified
    /// slices of shell, without integrator components.
    ///
    /// This function does not check the validity of the `shls_slice`.
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// let shls_slice = [[0, 15], [5, 19], [3, 12]];
    /// let cgto_shape = cint_data.cgto_shape_s1(&shls_slice);
    /// assert_eq!(cgto_shape, vec![37, 38, 29]);
    /// ```
    pub fn cgto_shape_s1(&self, shls_slice: &[[c_int; 2]]) -> Vec<usize> {
        let ao_loc = self.ao_loc();

        shls_slice
            .iter()
            .map(|&[shl0, shl1]| {
                let p0 = ao_loc[shl0 as usize];
                let p1 = ao_loc[shl1 as usize];
                p1 - p0
            })
            .collect()
    }

    /// Shape of integral (in atomic orbital basis, symmetry s2ij), for
    /// specified slices of shell, without integrator components.
    ///
    /// This function does not check the validity of the `shls_slice`.
    ///
    /// # Panics
    ///
    /// This function only checks if the first two elements of `shls_slice` are
    /// the same. If not, it will panic.
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// let shls_slice = [[0, 15], [0, 15], [3, 12]];
    /// let cgto_shape = cint_data.cgto_shape_s2ij(&shls_slice);
    /// assert_eq!(cgto_shape, vec![703, 29]);
    /// ```
    pub fn cgto_shape_s2ij(&self, shls_slice: &[[c_int; 2]]) -> Vec<usize> {
        let ao_loc = self.ao_loc();
        let n_center = shls_slice.len();

        if shls_slice[0] != shls_slice[1] {
            panic!("shls_slice must have the same first two elements for S2ij symmetry");
        }

        let mut shape = vec![];
        {
            let shl_slice = shls_slice[0];
            let p0 = ao_loc[shl_slice[0] as usize];
            let p1 = ao_loc[shl_slice[1] as usize];
            let np = p1 - p0;
            shape.push(np * (np + 1) / 2);
        }
        for n in 2..n_center {
            let shl_slice = shls_slice[n];
            let p0 = ao_loc[shl_slice[0] as usize];
            let p1 = ao_loc[shl_slice[1] as usize];
            shape.push(p1 - p0);
        }
        shape
    }

    /// Shape of integral (in atomic orbital basis, symmetry s2kl), for
    /// specified slices of shell, without integrator components.
    ///
    /// This function does not check the validity of the `shls_slice`.
    pub fn cgto_shape_s2kl(&self, shls_slice: &[[c_int; 2]]) -> Vec<usize> {
        let ao_loc = self.ao_loc();

        if shls_slice[2] != shls_slice[3] {
            panic!("shls_slice must have the same last two elements for S2kl symmetry");
        }

        let mut shape = vec![];
        {
            // push slice 1
            let shl_slice = shls_slice[0];
            let p0 = ao_loc[shl_slice[0] as usize];
            let p1 = ao_loc[shl_slice[1] as usize];
            let np = p1 - p0;
            shape.push(np);
        }
        {
            // push slice 2
            let shl_slice = shls_slice[1];
            let p0 = ao_loc[shl_slice[0] as usize];
            let p1 = ao_loc[shl_slice[1] as usize];
            let np = p1 - p0;
            shape.push(np);
        }
        {
            // push slice 3/4
            let shl_slice = shls_slice[2];
            let p0 = ao_loc[shl_slice[0] as usize];
            let p1 = ao_loc[shl_slice[1] as usize];
            let np = p1 - p0;
            shape.push(np * (np + 1) / 2);
        }
        shape
    }

    /// Shape of integral (in atomic orbital basis, symmetry s4), for
    /// specified slices of shell, without integrator components.
    ///
    /// This function does not check the validity of the `shls_slice`.
    pub fn cgto_shape_s4(&self, shls_slice: &[[c_int; 2]]) -> Vec<usize> {
        let ao_loc = self.ao_loc();

        if shls_slice[0] != shls_slice[1] || shls_slice[2] != shls_slice[3] {
            panic!("shls_slice must have the same first two elements and last two elements for S4 symmetry");
        }

        let mut shape = vec![];
        {
            // push slice 1
            let shl_slice = shls_slice[0];
            let p0 = ao_loc[shl_slice[0] as usize];
            let p1 = ao_loc[shl_slice[1] as usize];
            let np = p1 - p0;
            shape.push(np * (np + 1) / 2);
        }
        {
            // push slice 2
            let shl_slice = shls_slice[2];
            let p0 = ao_loc[shl_slice[0] as usize];
            let p1 = ao_loc[shl_slice[1] as usize];
            let np = p1 - p0;
            shape.push(np * (np + 1) / 2);
        }
        shape
    }

    /// Shape of integral (in atomic orbital basis, symmetry s8), for
    /// specified slices of shell, without integrator components.
    ///
    /// This function does not check the validity of the `shls_slice`.
    pub fn cgto_shape_s8(&self, shls_slice: &[[c_int; 2]]) -> Vec<usize> {
        let ao_loc = self.ao_loc();

        if shls_slice[0] != shls_slice[1] || shls_slice[0] != shls_slice[2] || shls_slice[0] != shls_slice[3] {
            panic!("shls_slice must have the same first four elements for S8 symmetry");
        }

        let shl_slice = shls_slice[0];
        let p0 = ao_loc[shl_slice[0] as usize];
        let p1 = ao_loc[shl_slice[1] as usize];
        let np = p1 - p0;

        let np_s2 = np * (np + 1) / 2;
        let np_s8 = np_s2 * (np_s2 + 1) / 2;
        vec![np_s8]
    }

    /// Obtain atomic orbital locations for integral, starting from 0 (instead
    /// of the first basis index mapped from first shell, or to say that is
    /// relative).
    ///
    /// This function will be used when coping small blocks of integral tensors
    /// into large output tensor.
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// let shls_slice = [[0, 7], [2, 6], [15, 19]];
    /// let cgto_locs = cint_data.cgto_locs(&shls_slice);
    /// assert_eq!(cgto_locs, vec![
    ///     vec![0, 1, 2, 3, 4, 5, 8, 11],
    ///     vec![0, 1, 2, 3, 6],
    ///     vec![0, 1, 2, 3, 6]
    /// ]);
    /// ```
    pub fn cgto_locs(&self, shls_slice: &[[c_int; 2]]) -> Vec<Vec<usize>> {
        let ao_loc = self.ao_loc();
        shls_slice
            .iter()
            .map(|&[shl0, shl1]| {
                let [shl0, shl1] = [shl0 as usize, shl1 as usize];
                ao_loc[shl0..shl1 + 1].iter().map(|&x| x - ao_loc[shl0]).collect_vec()
            })
            .collect_vec()
    }

    /* #endregion */

    /* #region other preparation */

    /// Obtain cache size for integral.
    ///
    /// If the shell slice is not known to you currently, just pass empty
    /// `shls_slice = vec![]`, then it should give the maximum cache size
    /// for this molecule/intor.
    ///
    /// Please also note that cache should be allocated per thread.
    ///
    /// # PySCF equivalent
    ///
    /// `GTOmax_cache_size` in `fill_int2e.c`
    ///
    /// # Panics
    ///
    /// - Does not check the validity of the `integrator`.
    /// - Panic if `shls_slice` exceeds the number of shells in the molecule.
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// let integrator = CInt::get_integrator("int2e_ip1");
    /// let cache_size = cint_data.max_cache_size(&*integrator, &[]);
    /// ```
    pub fn max_cache_size(&self, integrator: &dyn Integrator, shls_slice: &[[c_int; 2]]) -> usize {
        let natm = self.natm() as c_int;
        let nbas = self.nbas() as c_int;
        let shls_min = shls_slice.iter().map(|x| x[0]).min().unwrap_or(0);
        let shls_max = shls_slice.iter().map(|x| x[1]).max().unwrap_or(nbas);

        (shls_min..shls_max)
            .into_par_iter()
            .map(|shl| unsafe {
                let shls = [shl; 4];
                match self.cint_type {
                    Spheric => integrator.integral_sph(
                        null_mut(),
                        null(),
                        shls.as_ptr(),
                        self.atm_ptr(),
                        natm,
                        self.bas_ptr(),
                        nbas,
                        self.env_ptr(),
                        null(),
                        null_mut(),
                    ),
                    Cartesian => integrator.integral_cart(
                        null_mut(),
                        null(),
                        shls.as_ptr(),
                        self.atm_ptr(),
                        natm,
                        self.bas_ptr(),
                        nbas,
                        self.env_ptr(),
                        null(),
                        null_mut(),
                    ),
                    Spinor => integrator.integral_spinor(
                        null_mut(),
                        null(),
                        shls.as_ptr(),
                        self.atm_ptr(),
                        natm,
                        self.bas_ptr(),
                        nbas,
                        self.env_ptr(),
                        null(),
                        null_mut(),
                    ),
                }
            })
            .max()
            .unwrap_or(0) as usize
    }

    /// Obtain maximum buffer size for integral.
    ///
    /// # PySCF equivalent
    ///
    /// `GTOmax_shell_dim` in `fill_int2e.c` for similar functionality (but
    /// different in result)
    ///
    /// # Panics
    ///
    /// - Does not check the validity of the `integrator`.
    /// - Panic if `shls_slice` exceeds the number of shells in the molecule.
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// let integrator = CInt::get_integrator("int3c2e_ip1");
    /// let shls_slice = [[0, 7], [2, 6], [15, 19]];
    /// let max_buffer_size = cint_data.max_buffer_size(&*integrator, &shls_slice);
    /// assert_eq!(max_buffer_size, 81);
    /// ```
    pub fn max_buffer_size(&self, integrator: &dyn Integrator, shls_slice: &[[c_int; 2]]) -> usize {
        let ao_loc = self.ao_loc();
        let n_comp = match self.cint_type {
            Spheric | Cartesian => integrator.n_comp(),
            Spinor => integrator.n_spinor_comp(),
        };

        let mut result = n_comp;
        for &[shl0, shl1] in shls_slice {
            let mut cgto_max = 0;
            for shl in shl0..shl1 {
                let shl = shl as usize;
                let p0 = ao_loc[shl];
                let p1 = ao_loc[shl + 1];
                cgto_max = cgto_max.max(p1 - p0);
            }
            result *= cgto_max;
        }
        result
    }

    /// Obtain optimizer for integral.
    ///
    /// # Panics
    ///
    /// - Does not check the validity of the `integrator`.
    pub fn get_optimizer(&self, integrator: &dyn Integrator) -> CIntOptimizer {
        match integrator.kind() {
            CIntKind::Int => {
                let atm_ptr = self.atm.as_ptr() as *const c_int;
                let bas_ptr = self.bas.as_ptr() as *const c_int;
                let env_ptr = self.env.as_ptr();
                let n_atm = self.atm.len() as c_int;
                let n_bas = self.bas.len() as c_int;
                let mut c_opt_ptr: *mut CINTOpt = null_mut();
                unsafe {
                    integrator.optimizer(&mut c_opt_ptr as *mut *mut CINTOpt as *mut *mut c_void, atm_ptr, n_atm, bas_ptr, n_bas, env_ptr);
                };
                CIntOptimizer::Int(c_opt_ptr)
            },
            CIntKind::Ecp => {
                // Check if this cint_data has been merged with ECP data
                let merged = if self.is_ecp_merged() { self } else { &self.merge_ecpbas() };
                let atm_ptr = merged.atm.as_ptr() as *const c_int;
                let bas_ptr = merged.bas.as_ptr() as *const c_int;
                let env_ptr = merged.env.as_ptr();
                let n_atm = merged.atm.len() as c_int;
                let n_bas = merged.bas.len() as c_int;
                let mut c_opt_ptr: *mut ECPOpt = null_mut();
                unsafe {
                    integrator.optimizer(&mut c_opt_ptr as *mut *mut ECPOpt as *mut *mut c_void, atm_ptr, n_atm, bas_ptr, n_bas, env_ptr);
                };
                CIntOptimizer::Ecp(c_opt_ptr)
            },
        }
    }

    /* #endregion */

    /* #region integral of shell block */

    /// Smallest unit of electron-integral function from libcint.
    ///
    /// This is not a safe wrapper function, though it is pure rust. Use with
    /// caution.
    ///
    /// - `out` - Output integral buffer, need to be allocated enough space
    ///   before calling this function and properly offsetted.
    /// - `shls` - shell indices, which size should be n_center length.
    /// - `cint_opt` - optimizer that should be obtained before calling this
    ///   function. If `None`, then no optimization will be performed and will
    ///   pass null pointer.
    /// - `shape` - **In general cases, it is not recommended to use.** If
    ///   output larger than 4GB, then libcint internal realization may
    ///   overflow.
    /// - `cache` - cache buffer, need to be allocated enough space before
    ///   calling this function; simply using `vec![]` should also works, which
    ///   lets libcint manages cache and efficiency decreases. See Also
    ///   [`CInt::max_cache_size`] for guide of properly allocate cache.
    ///
    /// # Safety
    ///
    /// This function does not check anything. Caller should handle checks:
    ///
    /// - Does not check type (float or complex) of `out`. Type check will be
    ///   checked in caller.
    /// - Does not check `shls` length and values.
    /// - Does not check `shape` length and values.
    /// - Does not check whether `cint_opt` is ECP or general integrator, which
    ///   may or may not corresponds to `integrator`'s type.
    pub unsafe fn integral_block<F>(
        &self,
        integrator: &dyn Integrator,
        out: &mut [F],
        shls: &[c_int],
        shape: &[c_int],
        cint_opt: Option<&CIntOptimizer>,
        cache: &mut [f64],
    ) where
        F: ComplexFloat,
    {
        let cache_ptr = match cache.is_empty() {
            true => null_mut(),
            false => cache.as_mut_ptr(),
        };
        let shape_ptr = match shape.is_empty() {
            true => null_mut(),
            false => shape.as_ptr() as *const c_int,
        };
        let opt_ptr = cint_opt.map(|opt| opt.as_ptr()).unwrap_or(null());
        match self.cint_type {
            Spheric => unsafe {
                integrator.integral_sph(
                    out.as_mut_ptr() as *mut f64,
                    shape_ptr,
                    shls.as_ptr(),
                    self.atm_ptr(),
                    self.natm() as c_int,
                    self.bas_ptr(),
                    self.nbas() as c_int,
                    self.env_ptr(),
                    opt_ptr,
                    cache_ptr,
                );
            },
            Cartesian => unsafe {
                integrator.integral_cart(
                    out.as_mut_ptr() as *mut f64,
                    shape_ptr,
                    shls.as_ptr(),
                    self.atm_ptr(),
                    self.natm() as c_int,
                    self.bas_ptr(),
                    self.nbas() as c_int,
                    self.env_ptr(),
                    opt_ptr,
                    cache_ptr,
                );
            },
            Spinor => unsafe {
                integrator.integral_spinor(
                    out.as_mut_ptr() as *mut c_void,
                    shape_ptr,
                    shls.as_ptr(),
                    self.atm_ptr(),
                    self.natm() as c_int,
                    self.bas_ptr(),
                    self.nbas() as c_int,
                    self.env_ptr(),
                    opt_ptr,
                    cache_ptr,
                );
            },
        }
    }

    /* #endregion */
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn playground() {
        let cint_data = init_h2o_def2_tzvp();
        let integrator = CInt::get_integrator("int2e_ip1");
        let cache_size = cint_data.max_cache_size(&*integrator, &[]);
        println!("Cache size: {cache_size}");
    }
}
