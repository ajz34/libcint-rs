//! Implementation of grids implementation at lower API level.

use crate::prelude::*;
use rayon::prelude::*;
use rstsr_common::layout::IndexedIterLayout;
use rstsr_common::prelude::*;

const BLKSIZE: usize = 312;

/// Implementation of grids implementation at lower API level
impl CInt {
    pub fn max_cache_size_grids(&self, integrator: &dyn Integrator, shls_slice: &[[c_int; 2]]) -> usize {
        let natm = self.natm() as c_int;
        let nbas = self.nbas() as c_int;
        let shls_min = shls_slice.iter().map(|x| x[0]).min().unwrap_or(0);
        let shls_max = shls_slice.iter().map(|x| x[1]).max().unwrap_or(nbas);

        (shls_min..shls_max)
            .into_par_iter()
            .map(|shl| unsafe {
                let shls = [shl, shl, 0, BLKSIZE as i32];
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

    #[doc(hidden)]
    pub fn integral_grids_inplace<F>(
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

        let n_comp = match self.cint_type {
            Spheric | Cartesian => integrator.n_comp(),
            Spinor => integrator.n_spinor_comp(),
        }; // number of components for intor

        // n_center must be 2 for int1e_grids
        let cgto_shape = self.cgto_shape_s1(shls_slice); // AO shape, without intor component
        let cgto_locs = self.cgto_locs(shls_slice); // AO relative locations mapped to shells, 0-indexed
        let n_grids = self.ngrids();
        let mut grid_locs = (0..n_grids).step_by(BLKSIZE).collect_vec();
        grid_locs.push(n_grids);

        // cache (thread local)
        let cache_size = self.max_cache_size_grids(integrator, shls_slice);
        let buffer_size = BLKSIZE * self.max_buffer_size(integrator, shls_slice);
        let thread_cache = (0..rayon::current_num_threads()).map(|_| unsafe { aligned_uninitialized_vec::<f64>(cache_size) }).collect_vec();
        let thread_buffer = (0..rayon::current_num_threads()).map(|_| unsafe { aligned_uninitialized_vec::<F>(buffer_size) }).collect_vec();

        /* #endregion */

        /* #region parallel integration generation */

        let out_shape = [n_grids, cgto_shape[0], cgto_shape[1], n_comp];

        const I: usize = 0; // index of first shell
        const J: usize = 1; // index of second shell

        let nidx_i = (shls_slice[I][1] - shls_slice[I][0]) as usize;
        let nidx_j = (shls_slice[J][1] - shls_slice[J][0]) as usize;
        let nidx_g = grid_locs.len() - 1;

        let iter_layout = [nidx_g, nidx_i, nidx_j].f();
        let iter_indices = IndexedIterLayout::new(&iter_layout, ColMajor).unwrap();

        iter_indices.into_iter().for_each(|([idx_g, idx_i, idx_j], _)| {
            let shl_i = idx_i as c_int + shls_slice[I][0];
            let shl_j = idx_j as c_int + shls_slice[J][0];
            let cgto_loc_i = cgto_locs[I][idx_i];
            let cgto_loc_j = cgto_locs[J][idx_j];
            let cgto_loc_g = grid_locs[idx_g];
            let cgto_i = cgto_locs[I][idx_i + 1] - cgto_loc_i;
            let cgto_j = cgto_locs[J][idx_j + 1] - cgto_loc_j;
            let cgto_g = grid_locs[idx_g + 1] - cgto_loc_g;

            let shls = [shl_i, shl_j, cgto_loc_g as c_int, (cgto_loc_g + cgto_g) as c_int];

            // prepare cache and buffer
            let thread_idx = rayon::current_thread_index().unwrap_or(0);
            let cache = unsafe { cast_mut_slice(&thread_cache[thread_idx]) };
            let buf = unsafe { cast_mut_slice(&thread_buffer[thread_idx]) };

            // call integral function
            unsafe { self.integral_block(integrator, buf, &shls, &[], cint_opt, cache) };

            // copy buffer to output slice
            let out = unsafe { cast_mut_slice(&*out) };
            let out_offsets = [cgto_loc_g, cgto_loc_i, cgto_loc_j, 0];
            let buf_shape = [cgto_g, cgto_i, cgto_j, n_comp];
            copy_f_4d_s1(out, &out_offsets, &out_shape, buf, &buf_shape);
        });

        /* #endregion */

        Ok(())
    }

    #[doc(hidden)]
    pub fn integral_grids_row_major_inplace<F>(
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

        let n_comp = match self.cint_type {
            Spheric | Cartesian => integrator.n_comp(),
            Spinor => integrator.n_spinor_comp(),
        }; // number of components for intor

        // n_center must be 2 for int1e_grids
        let cgto_shape = self.cgto_shape_s1(shls_slice); // AO shape, without intor component
        let cgto_locs = self.cgto_locs(shls_slice); // AO relative locations mapped to shells, 0-indexed
        let n_grids = self.ngrids();
        let mut grid_locs = (0..n_grids).step_by(BLKSIZE).collect_vec();
        grid_locs.push(n_grids);

        // cache (thread local)
        let cache_size = self.max_cache_size_grids(integrator, shls_slice);
        let buffer_size = BLKSIZE * self.max_buffer_size(integrator, shls_slice);
        let thread_cache = (0..rayon::current_num_threads()).map(|_| unsafe { aligned_uninitialized_vec::<f64>(cache_size) }).collect_vec();
        let thread_buffer = (0..rayon::current_num_threads()).map(|_| unsafe { aligned_uninitialized_vec::<F>(buffer_size) }).collect_vec();

        /* #endregion */

        /* #region parallel integration generation */

        let out_shape = [n_comp, n_grids, cgto_shape[0], cgto_shape[1]];

        const I: usize = 0; // index of first shell
        const J: usize = 1; // index of second shell

        let nidx_i = (shls_slice[I][1] - shls_slice[I][0]) as usize;
        let nidx_j = (shls_slice[J][1] - shls_slice[J][0]) as usize;
        let nidx_g = grid_locs.len() - 1;

        let iter_layout = [nidx_g, nidx_i, nidx_j].c();
        let iter_indices = IndexedIterLayout::new(&iter_layout, RowMajor).unwrap();

        iter_indices.into_iter().for_each(|([idx_g, idx_i, idx_j], _)| {
            let shl_i = idx_i as c_int + shls_slice[I][0];
            let shl_j = idx_j as c_int + shls_slice[J][0];
            let cgto_loc_i = cgto_locs[I][idx_i];
            let cgto_loc_j = cgto_locs[J][idx_j];
            let cgto_loc_g = grid_locs[idx_g];
            let cgto_i = cgto_locs[I][idx_i + 1] - cgto_loc_i;
            let cgto_j = cgto_locs[J][idx_j + 1] - cgto_loc_j;
            let cgto_g = grid_locs[idx_g + 1] - cgto_loc_g;

            let shls = [shl_i, shl_j, cgto_loc_g as c_int, (cgto_loc_g + cgto_g) as c_int];

            // prepare cache and buffer
            let thread_idx = rayon::current_thread_index().unwrap_or(0);
            let cache = unsafe { cast_mut_slice(&thread_cache[thread_idx]) };
            let buf = unsafe { cast_mut_slice(&thread_buffer[thread_idx]) };

            // call integral function
            unsafe { self.integral_block(integrator, buf, &shls, &[], cint_opt, cache) };

            // copy buffer to output slice
            let out = unsafe { cast_mut_slice(&*out) };
            let out_offsets = [0, cgto_loc_g, cgto_loc_i, cgto_loc_j];
            let buf_shape = [cgto_g, cgto_i, cgto_j, n_comp];
            copy_c_4d_s1(out, &out_offsets, &out_shape, buf, &buf_shape);
        });

        /* #endregion */

        Ok(())
    }
}
