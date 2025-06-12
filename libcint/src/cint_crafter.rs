//! Integral crafter for `CInt` instance (the main integral functions).
//!
//! Code naming conventions:
//! - `cgto_`: shell slices related variables, indicating shell -> basis
//!   mapping.

#![allow(dead_code)]

use crate::prelude::*;
use rayon::prelude::*;
use rstsr_common::layout::IndexedIterLayout;
use rstsr_common::prelude::*;

/// Implementation of integral at higher API level (for user usage).
impl CInt {
    /// Obtain integrator optimizer.
    ///
    /// # Panics
    ///
    /// - integrator not found
    /// - `env` field not properly initialized for ECP integrator (should call
    ///   `merge_ecpbas` before this function)
    ///
    /// # See also
    ///
    /// PySCF `make_cintopt`
    pub fn optimizer(&self, intor: &str) -> CIntOptimizer {
        self.optimizer_f(intor).unwrap()
    }

    pub fn optimizer_f(&self, intor: &str) -> Result<CIntOptimizer, CIntError> {
        let integrator = crate::cint::get_integrator_f(intor)?;
        // make dummy shell slices and symms, only for checking
        let n_center = integrator.n_center();
        let cint_symm = CIntSymm::S1;
        let shls_slice = vec![[0, 0]; n_center];
        self.check_shls_slice(&*integrator, &shls_slice, cint_symm)?;
        Ok(self.get_optimizer(&*integrator))
    }

    pub fn intor(
        &self,
        intor: &str,
        aosym: impl Into<CIntSymm>,
        shls_slice: Option<&[[usize; 2]]>,
    ) -> CIntOutput<f64> {
        let intor_args = self
            .intor_args_builder()
            .intor(intor)
            .aosym(aosym)
            .shls_slice(shls_slice.unwrap_or(&[]))
            .build()
            .unwrap();
        self.intor_with_args_inner(intor_args).unwrap()
    }

    pub fn intor_spinor(
        &self,
        intor: &str,
        aosym: impl Into<CIntSymm>,
        shls_slice: Option<&[[usize; 2]]>,
    ) -> CIntOutput<Complex<f64>> {
        let intor_args = self
            .intor_args_builder_spinor()
            .intor(intor)
            .aosym(aosym)
            .shls_slice(shls_slice.unwrap_or(&[]))
            .build()
            .unwrap();
        self.intor_with_args_inner(intor_args).unwrap()
    }

    pub fn intor_with_args_inner<F>(&self, args: IntorArgs<F>) -> Result<CIntOutput<F>, CIntError>
    where
        F: ComplexFloat + Send + Sync,
    {
        let data = if self.is_ecp_merged() { &self.merge_ecpbas() } else { self };

        // unwrap and prepare (without output checking)
        let integrator = CInt::get_integrator(args.intor);
        let mut shls_slice = args
            .shls_slice
            .iter()
            .map(|&[shl0, shl1]| [shl0 as c_int, shl1 as c_int])
            .collect_vec();
        if shls_slice.is_empty() {
            shls_slice = vec![[0, data.nbas() as c_int]];
        }
        let aosym = args.aosym;
        let cint_opt = data.get_optimizer(&*integrator);

        // perfrom most checks
        data.check_float_type::<F>()?;
        data.check_shls_slice(&*integrator, &shls_slice, aosym)?;
        data.check_optimizer(&*integrator, &cint_opt)?;

        // prepare output and check size of output
        let mut out_shape = data.cgto_shape(&shls_slice, aosym);
        if integrator.n_comp() > 1 {
            out_shape.push(integrator.n_comp());
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
            return Err(CIntError::InvalidValue(format!(
                "Output vector size {} is smaller than required size {}",
                out.len(),
                out_size
            )));
        }

        // actual integral execution
        data.integral_inplace(&*integrator, out, &shls_slice, Some(&cint_opt), aosym)?;

        Ok(CIntOutput { out: out_vec, shape: out_shape })
    }
}

/// Implementation of integral at higher API level (legacy support).
impl CInt {
    pub fn integral_s1<T>(&self, shls_slice: Option<&[[c_int; 2]]>) -> Vec<f64>
    where
        T: Integrator + Default + ComplexFloat + Send + Sync,
    {
        if self.cint_type == Spinor {
            panic!("Spinor should be called by `integral_s1_spinor<Integrator>`");
        }

        let (out, _shape) = self.integral_inner::<T, f64>(shls_slice, CIntSymm::S1);
        out
    }

    pub fn integral_s1_spinor<T>(&self, shls_slice: Option<&[[c_int; 2]]>) -> Vec<Complex<f64>>
    where
        T: Integrator + Default + ComplexFloat + Send + Sync,
    {
        if self.cint_type != Spinor {
            panic!("Non-spinor should be called by `integral_s1<Integrator>`");
        }

        let (out, _shape) = self.integral_inner::<T, Complex<f64>>(shls_slice, CIntSymm::S1);
        out
    }

    pub fn integral_s2ij<T>(&self, shls_slice: Option<&[[c_int; 2]]>) -> Vec<f64>
    where
        T: Integrator + Default + ComplexFloat + Send + Sync,
    {
        if self.cint_type == Spinor {
            panic!("Spinor should be called by `integral_s2ij_spinor<Integrator>`");
        }

        let (out, _shape) = self.integral_inner::<T, f64>(shls_slice, CIntSymm::S2ij);
        out
    }

    pub fn integral_s2ij_spinor<T>(&self, shls_slice: Option<&[[c_int; 2]]>) -> Vec<Complex<f64>>
    where
        T: Integrator + Default + ComplexFloat + Send + Sync,
    {
        if self.cint_type != Spinor {
            panic!("Non-spinor should be called by `integral_s2ij<Integrator>`");
        }
        eprintln!(
            "`integral_spinor_s2ij` should generally not be called, since spinor may not show s2ij symmetry."
        );

        let (out, _shape) = self.integral_inner::<T, Complex<f64>>(shls_slice, CIntSymm::S2ij);
        out
    }

    pub fn integral_inner<T, F>(
        &self,
        shls_slice: Option<&[[c_int; 2]]>,
        aosym: CIntSymm,
    ) -> (Vec<F>, Vec<usize>)
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
        data.check_shls_slice(&integrator, shls_slice, aosym).unwrap();

        // integral preparation and execution
        let mut out_shape = data.cgto_shape(shls_slice, aosym);
        if integrator.n_comp() > 1 {
            out_shape.push(integrator.n_comp())
        };
        let out_size = out_shape.iter().product();
        let mut out = unsafe { aligned_uninitialized_vec::<F>(out_size) };
        let cint_opt = data.get_optimizer(&integrator);
        data.integral_inplace(&integrator, &mut out, shls_slice, Some(&cint_opt), aosym).unwrap();
        (out, out_shape)
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
    pub fn check_shls_slice(
        &self,
        integrator: &dyn Integrator,
        shls_slice: &[[c_int; 2]],
        cint_symm: CIntSymm,
    ) -> Result<(), CIntError> {
        let n_center = integrator.n_center();

        // length of shls_slice must be the same value to number of centers of
        // integrator.
        if shls_slice.len() != n_center {
            return Err(CIntError::IntegratorNotAvailable(format!(
                "Integrator {} requires {} centers, but got {}",
                integrator.name(),
                integrator.n_center(),
                shls_slice.len()
            )));
        }

        // Some integrators may not support certain shell types.
        match self.cint_type {
            Spheric => {
                if !integrator.is_sph_available() {
                    return Err(CIntError::IntegratorNotAvailable(format!(
                        "Integrator {} does not support spherical integrals",
                        integrator.name()
                    )));
                }
            },
            Cartesian => {
                if !integrator.is_cart_available() {
                    return Err(CIntError::IntegratorNotAvailable(format!(
                        "Integrator {} does not support Cartesian integrals",
                        integrator.name()
                    )));
                }
            },
            Spinor => {
                if !integrator.is_spinor_available() {
                    return Err(CIntError::IntegratorNotAvailable(format!(
                        "Integrator {} does not support spinor integrals",
                        integrator.name()
                    )));
                }
            },
        }

        // Every shell in shls_slice must be within the range of number of shells.
        let nbas = self.nbas() as c_int;
        for (i, shl) in shls_slice.iter().enumerate() {
            if !(0 <= shl[0] && shl[0] <= shl[1] && shl[1] <= nbas) {
                return Err(CIntError::InvalidValue(format!(
                    "Shell slice {:?} at index {} is not proper within [{}, {}]",
                    shl, i, 0, nbas
                )));
            }
        }

        // check symmetry
        match cint_symm {
            CIntSymm::S2ij => {
                if shls_slice[0] != shls_slice[1] {
                    return Err(CIntError::InvalidValue(format!(
                        "Integrator {} does not support S2ij symmetry with different shell slices ({:?}, {:?})",
                        integrator.name(),
                        shls_slice[0],
                        shls_slice[1]
                    )));
                }
            },
            CIntSymm::S2kl => {
                if n_center != 4 {
                    return Err(CIntError::InvalidValue(format!(
                        "Integrator {} requires 4 centers for S2kl symmetry, but got {}",
                        integrator.name(),
                        n_center
                    )));
                }
                if shls_slice[2] != shls_slice[3] {
                    return Err(CIntError::InvalidValue(format!(
                        "Integrator {} does not support S2kl symmetry with different shell slices ({:?}, {:?})",
                        integrator.name(),
                        shls_slice[2],
                        shls_slice[3]
                    )));
                }
            },
            CIntSymm::S4 => {
                if n_center != 4 {
                    return Err(CIntError::InvalidValue(format!(
                        "Integrator {} requires 4 centers for S4 symmetry, but got {}",
                        integrator.name(),
                        n_center
                    )));
                }
                if shls_slice[0] != shls_slice[1] || shls_slice[2] != shls_slice[3] {
                    return Err(CIntError::InvalidValue(format!(
                        "Integrator {} does not support S4 symmetry with different shell slices {:?}",
                        integrator.name(),
                        shls_slice
                    )));
                }
            },
            CIntSymm::S8 => {
                if n_center != 4 {
                    return Err(CIntError::InvalidValue(format!(
                        "Integrator {} requires 4 centers for S8 symmetry, but got {}",
                        integrator.name(),
                        n_center
                    )));
                }
                if shls_slice[0] != shls_slice[1]
                    || shls_slice[0] != shls_slice[2]
                    || shls_slice[0] != shls_slice[3]
                {
                    return Err(CIntError::InvalidValue(format!(
                        "Integrator {} does not support S8 symmetry with different shell slices {:?}",
                        integrator.name(),
                        shls_slice
                    )));
                }
            },
            CIntSymm::S1 => {
                // S1 symmetry is always available
            },
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
            Err(CIntError::InvalidValue(format!(
                "Expected float type size {expected} bytes, but got {actual} bytes"
            )))
        } else {
            Ok(())
        }
    }

    /// Check if the optimizer is compatible with the integrator.
    ///
    /// This check only involves that ECP and general integrals have different
    /// optimizers.
    pub fn check_optimizer(
        &self,
        integrator: &dyn Integrator,
        optimizer: &CIntOptimizer,
    ) -> Result<(), CIntError> {
        match optimizer {
            CIntOptimizer::Int(_) => {
                if integrator.kind() != CIntKind::Int {
                    return Err(CIntError::InvalidValue(format!(
                        "Optimizer is for Int, but integrator is {:?}",
                        integrator.kind()
                    )));
                }
            },
            CIntOptimizer::Ecp(_) => {
                if integrator.kind() != CIntKind::Ecp {
                    return Err(CIntError::InvalidValue(format!(
                        "Optimizer is for Ecp, but integrator is {:?}",
                        integrator.kind()
                    )));
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
            _ => {
                panic!("Currently only S1 and S2ij symmetries are supported.");
            },
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
            .unwrap() as usize
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
    /// assert_eq!(max_buffer_size, 243);
    /// ```
    pub fn max_buffer_size(&self, integrator: &dyn Integrator, shls_slice: &[[c_int; 2]]) -> usize {
        let ao_loc = self.ao_loc();
        let n_comp = integrator.n_comp();

        let mut result = n_comp;
        for &[shl0, shl1] in shls_slice {
            for shl in shl0..shl1 {
                let shl = shl as usize;
                let p0 = ao_loc[shl];
                let p1 = ao_loc[shl + 1];
                result *= p1 - p0;
            }
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
                    integrator.optimizer(
                        &mut c_opt_ptr as *mut *mut CINTOpt as *mut *mut c_void,
                        atm_ptr,
                        n_atm,
                        bas_ptr,
                        n_bas,
                        env_ptr,
                    );
                };
                CIntOptimizer::Int(NonNull::new(&mut c_opt_ptr).unwrap())
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
                    integrator.optimizer(
                        &mut c_opt_ptr as *mut *mut ECPOpt as *mut *mut c_void,
                        atm_ptr,
                        n_atm,
                        bas_ptr,
                        n_bas,
                        env_ptr,
                    );
                };
                CIntOptimizer::Ecp(NonNull::new(&mut c_opt_ptr).unwrap())
            },
        }
    }

    /* #endregion */
}

/// Implementation of integral at lower API level (crafting).
impl CInt {
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
    /// - `shape` - **In general cases, it is not recommanded to use.** If
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
        let cache_ptr = match cache.len() {
            0 => null_mut(),
            _ => cache.as_mut_ptr(),
        };
        let shape_ptr = match shape.len() {
            0 => null_mut(),
            _ => shape.as_ptr() as *const c_int,
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

    pub fn integral_inplace<F>(
        &self,
        integrator: &dyn Integrator,
        out: &mut [F],
        shls_slice: &[[c_int; 2]],
        cint_opt: Option<&CIntOptimizer>,
        aosym: CIntSymm,
    ) -> Result<(), CIntError>
    where
        F: ComplexFloat + Send + Sync,
    {
        match aosym {
            CIntSymm::S1 => self.integral_s1_inplace(integrator, out, shls_slice, cint_opt),
            CIntSymm::S2ij => self.integral_s2ij_inplace(integrator, out, shls_slice, cint_opt),
            _ => Err(CIntError::InvalidValue(format!("CIntSymm {aosym:?} is not supported"))),
        }
    }

    #[allow(non_snake_case)]
    pub fn integral_s1_inplace<F>(
        &self,
        integrator: &dyn Integrator,
        out: &mut [F],
        shls_slice: &[[c_int; 2]],
        cint_opt: Option<&CIntOptimizer>,
    ) -> Result<(), CIntError>
    where
        F: ComplexFloat + Send + Sync,
    {
        /* #region sanity check and preparation */

        self.check_float_type::<F>()?;
        self.check_shls_slice(integrator, shls_slice, CIntSymm::S1)?;
        if let Some(cint_opt) = cint_opt {
            self.check_optimizer(integrator, cint_opt)?;
        }

        // dimensions

        let n_comp = integrator.n_comp(); // number of components for intor
        let n_center = integrator.n_center(); // atom center number for intor
        let cgto_shape = self.cgto_shape_s1(shls_slice); // AO shape, without intor component
        let cgto_locs = self.cgto_locs(shls_slice); // AO relative locations mapped to shells, 0-indexed

        // cache (thread local)
        let cache_size = self.max_cache_size(integrator, shls_slice);
        let buffer_size = self.max_buffer_size(integrator, shls_slice);
        let thread_cache = (0..rayon::current_num_threads())
            .map(|_| Mutex::new(unsafe { aligned_uninitialized_vec::<f64>(cache_size) }))
            .collect_vec();
        let thread_buffer = (0..rayon::current_num_threads())
            .map(|_| Mutex::new(unsafe { aligned_uninitialized_vec::<F>(buffer_size) }))
            .collect_vec();

        /* #endregion */

        /* #region parallel integration generation */

        const I: usize = 0; // index of first shell
        const J: usize = 1; // index of second shell
        const K: usize = 2; // index of third shell
        const L: usize = 3; // index of fourth shell

        match n_center {
            2 => {
                let out_shape = [cgto_shape[I], cgto_shape[J], n_comp];

                let iter_layout = [
                    (shls_slice[I][1] - shls_slice[I][0]) as usize,
                    (shls_slice[J][1] - shls_slice[J][0]) as usize,
                ]
                .f();
                let iter_indices = IndexedIterLayout::new(&iter_layout, ColMajor).unwrap();

                iter_indices.into_par_iter().for_each(|([idx_i, idx_j], _)| {
                    // idx refers to the index of shell for iteration
                    // shl refers to the index of shell in the basis set (real shell index)
                    let shl_i = idx_i as c_int + shls_slice[I][0];
                    let shl_j = idx_j as c_int + shls_slice[J][0];
                    // cgto (ao basis) location for each shell
                    let cgto_loc_i = cgto_locs[I][idx_i];
                    let cgto_loc_j = cgto_locs[J][idx_j];
                    // number of cgto (ao basis) for each shell
                    let cgto_i = cgto_locs[I][idx_i + 1] - cgto_loc_i;
                    let cgto_j = cgto_locs[J][idx_j + 1] - cgto_loc_j;

                    let shls = [shl_i, shl_j];

                    // prepare cache and buffer
                    let thread_idx = rayon::current_thread_index().unwrap_or(0);
                    let mut cache = thread_cache[thread_idx].lock().unwrap();
                    let mut buf = thread_buffer[thread_idx].lock().unwrap();

                    // call integral function
                    unsafe {
                        self.integral_block(integrator, &mut buf, &shls, &[], cint_opt, &mut cache);
                    }

                    // copy buffer to output slice
                    let out = unsafe { cast_mut_slice(&*out) };
                    let out_offsets = [cgto_loc_i, cgto_loc_j, 0];
                    let buf_shape = [cgto_i, cgto_j, n_comp];
                    copy_3d_s1(out, &out_offsets, &out_shape, &buf, &buf_shape);
                });
            },
            3 => {
                let out_shape = [cgto_shape[I], cgto_shape[J], cgto_shape[K], n_comp];

                let iter_layout = [
                    (shls_slice[I][1] - shls_slice[I][0]) as usize,
                    (shls_slice[J][1] - shls_slice[J][0]) as usize,
                    (shls_slice[K][1] - shls_slice[K][0]) as usize,
                ]
                .f();
                let iter_indices = IndexedIterLayout::new(&iter_layout, ColMajor).unwrap();

                iter_indices.into_par_iter().for_each(|([idx_i, idx_j, idx_k], _)| {
                    let shl_i = idx_i as c_int + shls_slice[I][0];
                    let shl_j = idx_j as c_int + shls_slice[J][0];
                    let shl_k = idx_k as c_int + shls_slice[K][0];
                    let cgto_loc_i = cgto_locs[I][idx_i];
                    let cgto_loc_j = cgto_locs[J][idx_j];
                    let cgto_loc_k = cgto_locs[K][idx_k];
                    let cgto_i = cgto_locs[I][idx_i + 1] - cgto_loc_i;
                    let cgto_j = cgto_locs[J][idx_j + 1] - cgto_loc_j;
                    let cgto_k = cgto_locs[K][idx_k + 1] - cgto_loc_k;

                    let shls = [shl_i, shl_j, shl_k];

                    // prepare cache and buffer
                    let thread_idx = rayon::current_thread_index().unwrap_or(0);
                    let mut cache = thread_cache[thread_idx].lock().unwrap();
                    let mut buf = thread_buffer[thread_idx].lock().unwrap();

                    // call integral function
                    unsafe {
                        self.integral_block(integrator, &mut buf, &shls, &[], cint_opt, &mut cache);
                    }

                    // copy buffer to output slice
                    let out = unsafe { cast_mut_slice(&*out) };
                    let out_offsets = [cgto_loc_i, cgto_loc_j, cgto_loc_k, 0];
                    let buf_shape = [cgto_i, cgto_j, cgto_k, n_comp];
                    copy_4d_s1(out, &out_offsets, &out_shape, &buf, &buf_shape);
                });
            },
            4 => {
                let out_shape =
                    [cgto_shape[I], cgto_shape[J], cgto_shape[K], cgto_shape[L], n_comp];

                let iter_layout = [
                    (shls_slice[I][1] - shls_slice[I][0]) as usize,
                    (shls_slice[J][1] - shls_slice[J][0]) as usize,
                    (shls_slice[K][1] - shls_slice[K][0]) as usize,
                    (shls_slice[L][1] - shls_slice[L][0]) as usize,
                ]
                .f();
                let iter_indices = IndexedIterLayout::new(&iter_layout, ColMajor).unwrap();

                iter_indices.into_par_iter().for_each(|([idx_i, idx_j, idx_k, idx_l], _)| {
                    let shl_i = idx_i as c_int + shls_slice[I][0];
                    let shl_j = idx_j as c_int + shls_slice[J][0];
                    let shl_k = idx_k as c_int + shls_slice[K][0];
                    let shl_l = idx_l as c_int + shls_slice[L][0];
                    let cgto_loc_i = cgto_locs[I][idx_i];
                    let cgto_loc_j = cgto_locs[J][idx_j];
                    let cgto_loc_k = cgto_locs[K][idx_k];
                    let cgto_loc_l = cgto_locs[L][idx_l];
                    let cgto_i = cgto_locs[I][idx_i + 1] - cgto_loc_i;
                    let cgto_j = cgto_locs[J][idx_j + 1] - cgto_loc_j;
                    let cgto_k = cgto_locs[K][idx_k + 1] - cgto_loc_k;
                    let cgto_l = cgto_locs[L][idx_l + 1] - cgto_loc_l;

                    let shls = [shl_i, shl_j, shl_k, shl_l];

                    // prepare cache and buffer
                    let thread_idx = rayon::current_thread_index().unwrap_or(0);
                    let mut cache = thread_cache[thread_idx].lock().unwrap();
                    let mut buf = thread_buffer[thread_idx].lock().unwrap();

                    // call integral function
                    unsafe {
                        self.integral_block(integrator, &mut buf, &shls, &[], cint_opt, &mut cache);
                    }

                    // copy buffer to output slice
                    let out = unsafe { cast_mut_slice(&*out) };
                    let out_offsets = [cgto_loc_i, cgto_loc_j, cgto_loc_k, cgto_loc_l, 0];
                    let buf_shape = [cgto_i, cgto_j, cgto_k, cgto_l, n_comp];
                    copy_5d_s1(out, &out_offsets, &out_shape, &buf, &buf_shape);
                });
            },
            _ => unreachable!(),
        }

        /* #endregion */

        Ok(())
    }

    #[allow(non_snake_case)]
    pub fn integral_s2ij_inplace<F>(
        &self,
        integrator: &dyn Integrator,
        out: &mut [F],
        shls_slice: &[[c_int; 2]],
        cint_opt: Option<&CIntOptimizer>,
    ) -> Result<(), CIntError>
    where
        F: ComplexFloat + Send + Sync,
    {
        /* #region sanity check and preparation */

        self.check_float_type::<F>()?;
        self.check_shls_slice(integrator, shls_slice, CIntSymm::S1)?;
        if let Some(cint_opt) = cint_opt {
            self.check_optimizer(integrator, cint_opt)?;
        }

        // dimensions

        let n_comp = integrator.n_comp(); // number of components for intor
        let n_center = integrator.n_center(); // atom center number for intor
        let cgto_shape = self.cgto_shape_s2ij(shls_slice); // AO shape, without intor component
        let cgto_locs = self.cgto_locs(shls_slice); // AO relative locations mapped to shells, 0-indexed

        // cache (thread local)
        let cache_size = self.max_cache_size(integrator, shls_slice);
        let buffer_size = self.max_buffer_size(integrator, shls_slice);
        let thread_cache = (0..rayon::current_num_threads())
            .map(|_| Mutex::new(unsafe { aligned_uninitialized_vec::<f64>(cache_size) }))
            .collect_vec();
        let thread_buffer = (0..rayon::current_num_threads())
            .map(|_| Mutex::new(unsafe { aligned_uninitialized_vec::<F>(buffer_size) }))
            .collect_vec();

        /* #endregion */

        /* #region parallel integration generation */

        const I: usize = 0; // index of first shell
        const J: usize = 1; // index of second shell
        const K: usize = 2; // index of third shell
        const L: usize = 3; // index of fourth shell

        let nidx_i = (shls_slice[I][1] - shls_slice[I][0]) as usize;
        let nidx_ij = nidx_i * (nidx_i + 1) / 2;

        match n_center {
            2 => {
                let out_shape = [cgto_shape[0], n_comp];

                (0..nidx_ij).into_par_iter().for_each(|idx_ij| {
                    let [idx_i, idx_j] = unravel_s2_indices(idx_ij);

                    let shl_i = idx_i as c_int + shls_slice[I][0];
                    let shl_j = idx_j as c_int + shls_slice[J][0];
                    let cgto_loc_i = cgto_locs[I][idx_i];
                    let cgto_loc_j = cgto_locs[J][idx_j];
                    let cgto_i = cgto_locs[I][idx_i + 1] - cgto_loc_i;
                    let cgto_j = cgto_locs[J][idx_j + 1] - cgto_loc_j;

                    let shls = [shl_i, shl_j];

                    // prepare cache and buffer
                    let thread_idx = rayon::current_thread_index().unwrap_or(0);
                    let mut cache = thread_cache[thread_idx].lock().unwrap();
                    let mut buf = thread_buffer[thread_idx].lock().unwrap();

                    // call integral function
                    unsafe {
                        self.integral_block(integrator, &mut buf, &shls, &[], cint_opt, &mut cache);
                    }

                    // copy buffer to output slice
                    let out = unsafe { cast_mut_slice(&*out) };
                    let out_offsets = [cgto_loc_i, cgto_loc_j, 0];
                    let buf_shape = [cgto_i, cgto_j, n_comp];
                    copy_3d_s2ij(out, &out_offsets, &out_shape, &buf, &buf_shape);
                });
            },
            3 => {
                let out_shape = [cgto_shape[0], cgto_shape[1], n_comp];

                let iter_layout = [nidx_ij, (shls_slice[K][1] - shls_slice[K][0]) as usize].f();
                let iter_indices = IndexedIterLayout::new(&iter_layout, ColMajor).unwrap();

                iter_indices.into_par_iter().for_each(|([idx_ij, idx_k], _)| {
                    let [idx_i, idx_j] = unravel_s2_indices(idx_ij);

                    let shl_i = idx_i as c_int + shls_slice[I][0];
                    let shl_j = idx_j as c_int + shls_slice[J][0];
                    let shl_k = idx_k as c_int + shls_slice[K][0];
                    let cgto_loc_i = cgto_locs[I][idx_i];
                    let cgto_loc_j = cgto_locs[J][idx_j];
                    let cgto_loc_k = cgto_locs[K][idx_k];
                    let cgto_i = cgto_locs[I][idx_i + 1] - cgto_loc_i;
                    let cgto_j = cgto_locs[J][idx_j + 1] - cgto_loc_j;
                    let cgto_k = cgto_locs[K][idx_k + 1] - cgto_loc_k;

                    let shls = [shl_i, shl_j, shl_k];

                    // prepare cache and buffer
                    let thread_idx = rayon::current_thread_index().unwrap_or(0);
                    let mut cache = thread_cache[thread_idx].lock().unwrap();
                    let mut buf = thread_buffer[thread_idx].lock().unwrap();

                    // call integral function
                    unsafe {
                        self.integral_block(integrator, &mut buf, &shls, &[], cint_opt, &mut cache);
                    }

                    // copy buffer to output slice
                    let out = unsafe { cast_mut_slice(&*out) };
                    let out_offsets = [cgto_loc_i, cgto_loc_j, cgto_loc_k, 0];
                    let buf_shape = [cgto_i, cgto_j, cgto_k, n_comp];
                    copy_4d_s2ij(out, &out_offsets, &out_shape, &buf, &buf_shape);
                });
            },
            4 => {
                let out_shape = [cgto_shape[0], cgto_shape[1], cgto_shape[2], n_comp];

                let iter_layout = [
                    nidx_ij,
                    (shls_slice[K][1] - shls_slice[K][0]) as usize,
                    (shls_slice[L][1] - shls_slice[L][0]) as usize,
                ]
                .f();
                let iter_indices = IndexedIterLayout::new(&iter_layout, ColMajor).unwrap();

                iter_indices.into_par_iter().for_each(|([idx_ij, idx_k, idx_l], _)| {
                    let [idx_i, idx_j] = unravel_s2_indices(idx_ij);

                    let shl_i = idx_i as c_int + shls_slice[I][0];
                    let shl_j = idx_j as c_int + shls_slice[J][0];
                    let shl_k = idx_k as c_int + shls_slice[K][0];
                    let shl_l = idx_l as c_int + shls_slice[L][0];
                    let cgto_loc_i = cgto_locs[I][idx_i];
                    let cgto_loc_j = cgto_locs[J][idx_j];
                    let cgto_loc_k = cgto_locs[K][idx_k];
                    let cgto_loc_l = cgto_locs[L][idx_l];
                    let cgto_i = cgto_locs[I][idx_i + 1] - cgto_loc_i;
                    let cgto_j = cgto_locs[J][idx_j + 1] - cgto_loc_j;
                    let cgto_k = cgto_locs[K][idx_k + 1] - cgto_loc_k;
                    let cgto_l = cgto_locs[L][idx_l + 1] - cgto_loc_l;

                    let shls = [shl_i, shl_j, shl_k, shl_l];

                    // prepare cache and buffer
                    let thread_idx = rayon::current_thread_index().unwrap_or(0);
                    let mut cache = thread_cache[thread_idx].lock().unwrap();
                    let mut buf = thread_buffer[thread_idx].lock().unwrap();

                    // call integral function
                    unsafe {
                        self.integral_block(integrator, &mut buf, &shls, &[], cint_opt, &mut cache);
                    }

                    // copy buffer to output slice
                    let out = unsafe { cast_mut_slice(&*out) };
                    let out_offsets = [cgto_loc_i, cgto_loc_j, cgto_loc_k, cgto_loc_l, 0];
                    let buf_shape = [cgto_i, cgto_j, cgto_k, cgto_l, n_comp];
                    copy_5d_s2ij(out, &out_offsets, &out_shape, &buf, &buf_shape);
                });
            },
            _ => unreachable!(),
        }

        /* #endregion */

        Ok(())
    }
}

/* #region util */

#[allow(clippy::mut_from_ref)]
pub(crate) unsafe fn cast_mut_slice<T>(slc: &[T]) -> &mut [T] {
    let len = slc.len();
    let ptr = slc.as_ptr() as *mut T;
    unsafe { std::slice::from_raw_parts_mut(ptr, len) }
}

#[inline(always)]
pub(crate) fn unravel_s2_indices(x: usize) -> [usize; 2] {
    let j = (((x * 2 + 1) as f64).sqrt() - 0.5) as usize;
    let i = x - j * (j + 1) / 2;
    [i, j]
}

#[inline(always)]
pub(crate) fn get_f_index_3d(indices: &[usize; 3], shape: &[usize; 3]) -> usize {
    indices[0] + shape[0] * (indices[1] + shape[1] * (indices[2]))
}

#[inline(always)]
pub(crate) fn get_f_index_4d(indices: &[usize; 4], shape: &[usize; 4]) -> usize {
    indices[0] + shape[0] * (indices[1] + shape[1] * (indices[2] + shape[2] * (indices[3])))
}

#[inline(always)]
pub(crate) fn get_f_index_5d(indices: &[usize; 5], shape: &[usize; 5]) -> usize {
    indices[0]
        + shape[0]
            * (indices[1]
                + shape[1] * (indices[2] + shape[2] * (indices[3] + shape[3] * (indices[4]))))
}

#[inline(always)]
pub(crate) fn get_f_index_3d_s2ij(indices: &[usize; 3], shape: &[usize; 2]) -> usize {
    indices[0] + indices[1] * (indices[1] + 1) / 2 + shape[0] * (indices[2])
}

#[inline(always)]
pub(crate) fn get_f_index_4d_s2ij(indices: &[usize; 4], shape: &[usize; 3]) -> usize {
    indices[0]
        + indices[1] * (indices[1] + 1) / 2
        + shape[0] * (indices[2] + shape[1] * (indices[3]))
}

#[inline(always)]
pub(crate) fn get_f_index_5d_s2ij(indices: &[usize; 5], shape: &[usize; 4]) -> usize {
    indices[0]
        + indices[1] * (indices[1] + 1) / 2
        + shape[0] * (indices[2] + shape[1] * (indices[3] + shape[2] * (indices[4])))
}

#[inline(always)]
pub(crate) fn copy_3d_s1<T>(
    out: &mut [T],
    out_offsets: &[usize; 3],
    out_shape: &[usize; 3],
    buf: &[T],
    buf_shape: &[usize; 3],
) where
    T: Copy,
{
    for c in 0..buf_shape[2] {
        for j in 0..buf_shape[1] {
            for i in 0..buf_shape[0] {
                let out_indices = [out_offsets[0] + i, out_offsets[1] + j, out_offsets[2] + c];
                let buf_indices = [i, j, c];
                let out_index = get_f_index_3d(&out_indices, out_shape);
                let buf_index = get_f_index_3d(&buf_indices, buf_shape);
                out[out_index] = buf[buf_index];
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_4d_s1<T>(
    out: &mut [T],
    out_offsets: &[usize; 4],
    out_shape: &[usize; 4],
    buf: &[T],
    buf_shape: &[usize; 4],
) where
    T: Copy,
{
    for c in 0..buf_shape[3] {
        for k in 0..buf_shape[2] {
            for j in 0..buf_shape[1] {
                for i in 0..buf_shape[0] {
                    let out_indices = [
                        out_offsets[0] + i,
                        out_offsets[1] + j,
                        out_offsets[2] + k,
                        out_offsets[3] + c,
                    ];
                    let buf_indices = [i, j, k, c];
                    let out_index = get_f_index_4d(&out_indices, out_shape);
                    let buf_index = get_f_index_4d(&buf_indices, buf_shape);
                    out[out_index] = buf[buf_index];
                }
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_5d_s1<T>(
    out: &mut [T],
    out_offsets: &[usize; 5],
    out_shape: &[usize; 5],
    buf: &[T],
    buf_shape: &[usize; 5],
) where
    T: Copy,
{
    for c in 0..buf_shape[4] {
        for l in 0..buf_shape[3] {
            for k in 0..buf_shape[2] {
                for j in 0..buf_shape[1] {
                    for i in 0..buf_shape[0] {
                        let out_indices = [
                            out_offsets[0] + i,
                            out_offsets[1] + j,
                            out_offsets[2] + k,
                            out_offsets[3] + l,
                            out_offsets[4] + c,
                        ];
                        let buf_indices = [i, j, k, l, c];
                        let out_index = get_f_index_5d(&out_indices, out_shape);
                        let buf_index = get_f_index_5d(&buf_indices, buf_shape);
                        out[out_index] = buf[buf_index];
                    }
                }
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_3d_s2ij<T>(
    out: &mut [T],
    out_offsets: &[usize; 3],
    out_shape: &[usize; 2],
    buf: &[T],
    buf_shape: &[usize; 3],
) where
    T: Copy,
{
    if out_offsets[0] != out_offsets[1] {
        // offdiag parts
        for c in 0..buf_shape[2] {
            for j in 0..buf_shape[1] {
                for i in 0..buf_shape[0] {
                    let out_indices = [out_offsets[0] + i, out_offsets[1] + j, out_offsets[2] + c];
                    let buf_indices = [i, j, c];
                    let out_index = get_f_index_3d_s2ij(&out_indices, out_shape);
                    let buf_index = get_f_index_3d(&buf_indices, buf_shape);
                    out[out_index] = buf[buf_index];
                }
            }
        }
    } else {
        // diag parts
        for c in 0..buf_shape[2] {
            for j in 0..buf_shape[1] {
                for i in 0..(j + 1) {
                    let out_indices = [out_offsets[0] + i, out_offsets[1] + j, out_offsets[2] + c];
                    let buf_indices = [i, j, c];
                    let out_index = get_f_index_3d_s2ij(&out_indices, out_shape);
                    let buf_index = get_f_index_3d(&buf_indices, buf_shape);
                    out[out_index] = buf[buf_index];
                }
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_4d_s2ij<T>(
    out: &mut [T],
    out_offsets: &[usize; 4],
    out_shape: &[usize; 3],
    buf: &[T],
    buf_shape: &[usize; 4],
) where
    T: Copy,
{
    if out_offsets[0] != out_offsets[1] {
        // offdiag parts
        for c in 0..buf_shape[3] {
            for k in 0..buf_shape[2] {
                for j in 0..buf_shape[1] {
                    for i in 0..buf_shape[0] {
                        let out_indices = [
                            out_offsets[0] + i,
                            out_offsets[1] + j,
                            out_offsets[2] + k,
                            out_offsets[3] + c,
                        ];
                        let buf_indices = [i, j, k, c];
                        let out_index = get_f_index_4d_s2ij(&out_indices, out_shape);
                        let buf_index = get_f_index_4d(&buf_indices, buf_shape);
                        out[out_index] = buf[buf_index];
                    }
                }
            }
        }
    } else {
        // diag parts
        for c in 0..buf_shape[3] {
            for k in 0..buf_shape[2] {
                for j in 0..buf_shape[1] {
                    for i in 0..(j + 1) {
                        let out_indices = [
                            out_offsets[0] + i,
                            out_offsets[1] + j,
                            out_offsets[2] + k,
                            out_offsets[3] + c,
                        ];
                        let buf_indices = [i, j, k, c];
                        let out_index = get_f_index_4d_s2ij(&out_indices, out_shape);
                        let buf_index = get_f_index_4d(&buf_indices, buf_shape);
                        out[out_index] = buf[buf_index];
                    }
                }
            }
        }
    }
}

#[inline(always)]
pub(crate) fn copy_5d_s2ij<T>(
    out: &mut [T],
    out_offsets: &[usize; 5],
    out_shape: &[usize; 4],
    buf: &[T],
    buf_shape: &[usize; 5],
) where
    T: Copy,
{
    if out_offsets[0] != out_offsets[1] {
        // offdiag parts
        for c in 0..buf_shape[4] {
            for l in 0..buf_shape[3] {
                for k in 0..buf_shape[2] {
                    for j in 0..buf_shape[1] {
                        for i in 0..buf_shape[0] {
                            let out_indices = [
                                out_offsets[0] + i,
                                out_offsets[1] + j,
                                out_offsets[2] + k,
                                out_offsets[3] + l,
                                out_offsets[4] + c,
                            ];
                            let buf_indices = [i, j, k, l, c];
                            let out_index = get_f_index_5d_s2ij(&out_indices, out_shape);
                            let buf_index = get_f_index_5d(&buf_indices, buf_shape);
                            out[out_index] = buf[buf_index];
                        }
                    }
                }
            }
        }
    } else {
        // diag parts
        for c in 0..buf_shape[4] {
            for l in 0..buf_shape[3] {
                for k in 0..buf_shape[2] {
                    for j in 0..buf_shape[1] {
                        for i in 0..(j + 1) {
                            let out_indices = [
                                out_offsets[0] + i,
                                out_offsets[1] + j,
                                out_offsets[2] + k,
                                out_offsets[3] + l,
                                out_offsets[4] + c,
                            ];
                            let buf_indices = [i, j, k, l, c];
                            let out_index = get_f_index_5d_s2ij(&out_indices, out_shape);
                            let buf_index = get_f_index_5d(&buf_indices, buf_shape);
                            out[out_index] = buf[buf_index];
                        }
                    }
                }
            }
        }
    }
}

/// Create an unaligned uninitialized vector with the given size.
///
/// # Safety
///
/// Caller must ensure that the vector is properly initialized before using it.
///
/// This is not a very good function, since `set_len` on uninitialized memory is
/// undefined-behavior (UB).
/// Nevertheless, if `T` is some type of `MaybeUninit`, then this will not UB.
#[allow(clippy::uninit_vec)]
#[inline]
pub unsafe fn unaligned_uninitialized_vec<T>(size: usize) -> Vec<T> {
    let mut v: Vec<T> = vec![];
    v.try_reserve_exact(size).unwrap();
    unsafe { v.set_len(size) };
    v
}

/// Create an uninitialized vector with the given size and alignment.
///
/// - None: if the size is 0 or allocation fails.
/// - Some: pointer to the allocated memory.
///
/// https://users.rust-lang.org/t/how-can-i-allocate-aligned-memory-in-rust/33293
#[inline]
pub fn aligned_alloc(numbytes: usize, alignment: usize) -> Option<NonNull<()>> {
    extern crate alloc;

    if numbytes == 0 {
        return None;
    }
    let layout = alloc::alloc::Layout::from_size_align(numbytes, alignment).unwrap();
    NonNull::new(unsafe { alloc::alloc::alloc(layout) }).map(|p| p.cast::<()>())
}

/// Create an conditionally aligned uninitialized vector with the given size.
///
/// # Safety
///
/// Caller must ensure that the vector is properly initialized before using it.
///
/// This is not a very good function, since `set_len` on uninitialized memory is
/// undefined-behavior (UB).
/// Nevertheless, if `T` is some type of `MaybeUninit`, then this will not UB.
#[inline]
pub unsafe fn aligned_uninitialized_vec<T>(size: usize) -> Vec<T> {
    const MIN_ALIGN: usize = 64; // minimal number of elements in vector to be aligned
    const ALIGNMENT: usize = 64; // 64 bytes alignment (minimal requirement for AVX-512)

    if size == 0 {
        vec![]
    } else if size < MIN_ALIGN {
        unsafe { unaligned_uninitialized_vec(size) }
    } else {
        let sizeof = core::mem::size_of::<T>();
        let pointer = aligned_alloc(size * sizeof, ALIGNMENT);
        if let Some(pointer) = pointer {
            let mut v = unsafe { Vec::from_raw_parts(pointer.as_ptr() as *mut T, size, size) };
            unsafe { v.set_len(size) };
            return v;
        } else {
            panic!("Allocation failed (probably due to out-of-memory)")
        }
    }
}

/* #endregion */

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
