//! Implementation of integral at lower API level (crafting, row-major).

use crate::prelude::*;
use rayon::prelude::*;
use rstsr_common::layout::IndexedIterLayout;
use rstsr_common::prelude::*;

/// Implementation of integral at lower API level (crafting, row-major).
impl CInt {
    #[doc(hidden)]
    pub fn integral_row_major_inplace<F>(
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
            CIntSymm::S1 => self.integral_s1_row_major_inplace(integrator, out, shls_slice, cint_opt),
            CIntSymm::S2ij => self.integral_s2ij_row_major_inplace(integrator, out, shls_slice, cint_opt),
            CIntSymm::S2kl => self.integral_s2kl_row_major_inplace(integrator, out, shls_slice, cint_opt),
            CIntSymm::S4 => self.integral_s4_row_major_inplace(integrator, out, shls_slice, cint_opt),
            CIntSymm::S8 => self.integral_s8_row_major_inplace(integrator, out, shls_slice, cint_opt),
        }
    }

    #[doc(hidden)]
    pub fn integral_s1_row_major_inplace<F>(
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
        let n_center = integrator.n_center(); // atom center number for intor
        let cgto_shape = self.cgto_shape_s1(shls_slice); // AO shape, without intor component
        let cgto_locs = self.cgto_locs(shls_slice); // AO relative locations mapped to shells, 0-indexed

        // cache (thread local)
        let cache_size = self.max_cache_size(integrator, shls_slice);
        let buffer_size = self.max_buffer_size(integrator, shls_slice);
        let thread_cache = (0..rayon::current_num_threads()).map(|_| unsafe { aligned_uninitialized_vec::<f64>(cache_size) }).collect_vec();
        let thread_buffer = (0..rayon::current_num_threads()).map(|_| unsafe { aligned_uninitialized_vec::<F>(buffer_size) }).collect_vec();

        /* #endregion */

        /* #region parallel integration generation */

        const I: usize = 0; // index of first shell
        const J: usize = 1; // index of second shell
        const K: usize = 2; // index of third shell
        const L: usize = 3; // index of fourth shell

        let nidx_i = (shls_slice[I][1] - shls_slice[I][0]) as usize;
        let nidx_j = (shls_slice[J][1] - shls_slice[J][0]) as usize;

        match n_center {
            2 => {
                let out_shape = [n_comp, cgto_shape[I], cgto_shape[J]];
                let iter_layout = [nidx_i, nidx_j].c();
                let iter_indices = IndexedIterLayout::new(&iter_layout, RowMajor).unwrap();

                iter_indices.into_par_iter().for_each(|([idx_i, idx_j], _)| {
                    let shl_i = idx_i as c_int + shls_slice[I][0];
                    let shl_j = idx_j as c_int + shls_slice[J][0];
                    let cgto_loc_i = cgto_locs[I][idx_i];
                    let cgto_loc_j = cgto_locs[J][idx_j];
                    let cgto_i = cgto_locs[I][idx_i + 1] - cgto_loc_i;
                    let cgto_j = cgto_locs[J][idx_j + 1] - cgto_loc_j;

                    let shls = [shl_i, shl_j];

                    // prepare cache and buffer
                    let thread_idx = rayon::current_thread_index().unwrap_or(0);
                    let cache = unsafe { cast_mut_slice(&thread_cache[thread_idx]) };
                    let buf = unsafe { cast_mut_slice(&thread_buffer[thread_idx]) };

                    // call integral function
                    unsafe { self.integral_block(integrator, buf, &shls, &[], cint_opt, cache) };

                    // copy buffer to output slice
                    let out = unsafe { cast_mut_slice(&*out) };
                    let out_offsets = [0, cgto_loc_i, cgto_loc_j];
                    let buf_shape = [cgto_i, cgto_j, n_comp];
                    copy_c_3d_s1(out, &out_offsets, &out_shape, buf, &buf_shape);
                });
            },
            3 => {
                let out_shape = [n_comp, cgto_shape[I], cgto_shape[J], cgto_shape[K]];

                let nidx_k = (shls_slice[K][1] - shls_slice[K][0]) as usize;
                let iter_layout = [nidx_i, nidx_j, nidx_k].c();
                let iter_indices = IndexedIterLayout::new(&iter_layout, RowMajor).unwrap();

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
                    let cache = unsafe { cast_mut_slice(&thread_cache[thread_idx]) };
                    let buf = unsafe { cast_mut_slice(&thread_buffer[thread_idx]) };

                    // call integral function
                    unsafe { self.integral_block(integrator, buf, &shls, &[], cint_opt, cache) };

                    // copy buffer to output slice
                    let out = unsafe { cast_mut_slice(&*out) };
                    let out_offsets = [0, cgto_loc_i, cgto_loc_j, cgto_loc_k];
                    let buf_shape = [cgto_i, cgto_j, cgto_k, n_comp];
                    copy_c_4d_s1(out, &out_offsets, &out_shape, buf, &buf_shape);
                });
            },
            4 => {
                let out_shape = [n_comp, cgto_shape[I], cgto_shape[J], cgto_shape[K], cgto_shape[L]];

                let nidx_k = (shls_slice[K][1] - shls_slice[K][0]) as usize;
                let nidx_l = (shls_slice[L][1] - shls_slice[L][0]) as usize;
                let iter_layout = [nidx_i, nidx_j, nidx_k, nidx_l].c();
                let iter_indices = IndexedIterLayout::new(&iter_layout, RowMajor).unwrap();

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
                    let cache = unsafe { cast_mut_slice(&thread_cache[thread_idx]) };
                    let buf = unsafe { cast_mut_slice(&thread_buffer[thread_idx]) };

                    // call integral function
                    unsafe { self.integral_block(integrator, buf, &shls, &[], cint_opt, cache) };

                    // copy buffer
                    let out = unsafe { cast_mut_slice(&*out) };
                    let out_offsets = [0, cgto_loc_i, cgto_loc_j, cgto_loc_k, cgto_loc_l];
                    let buf_shape = [cgto_i, cgto_j, cgto_k, cgto_l, n_comp];
                    copy_c_5d_s1(out, &out_offsets, &out_shape, buf, &buf_shape);
                });
            },
            _ => unreachable!(),
        }

        /* #endregion */

        Ok(())
    }

    #[doc(hidden)]
    pub fn integral_s2ij_row_major_inplace<F>(
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
        let n_center = integrator.n_center(); // atom center number for intor
        let cgto_shape = self.cgto_shape_s2ij(shls_slice); // AO shape, without intor component
        let cgto_locs = self.cgto_locs(shls_slice); // AO relative locations mapped to shells, 0-indexed

        // cache (thread local)
        let cache_size = self.max_cache_size(integrator, shls_slice);
        let buffer_size = self.max_buffer_size(integrator, shls_slice);
        let thread_cache = (0..rayon::current_num_threads()).map(|_| unsafe { aligned_uninitialized_vec::<f64>(cache_size) }).collect_vec();
        let thread_buffer = (0..rayon::current_num_threads()).map(|_| unsafe { aligned_uninitialized_vec::<F>(buffer_size) }).collect_vec();

        /* #endregion */

        /* #region parallel integration generation */

        const I: usize = 0; // index of first shell
        const J: usize = 1; // index of second shell
        const K: usize = 2; // index of third shell
        const L: usize = 3; // index of fourth shell

        let nidx_i = (shls_slice[I][1] - shls_slice[I][0]) as usize;
        let nidx_ij = nidx_i * (nidx_i + 1) / 2; // number of unique (i, j) pairs

        match n_center {
            2 => {
                let out_shape = [n_comp, cgto_shape[0]];

                (0..nidx_ij).for_each(|idx_ij| {
                    let [idx_j, idx_i] = unravel_s2_indices(idx_ij);

                    let shl_i = idx_i as c_int + shls_slice[I][0];
                    let shl_j = idx_j as c_int + shls_slice[J][0];
                    let cgto_loc_i = cgto_locs[I][idx_i];
                    let cgto_loc_j = cgto_locs[J][idx_j];
                    let cgto_i = cgto_locs[I][idx_i + 1] - cgto_loc_i;
                    let cgto_j = cgto_locs[J][idx_j + 1] - cgto_loc_j;

                    let shls = [shl_i, shl_j];

                    // prepare cache and buffer
                    let thread_idx = rayon::current_thread_index().unwrap_or(0);
                    let cache = unsafe { cast_mut_slice(&thread_cache[thread_idx]) };
                    let buf = unsafe { cast_mut_slice(&thread_buffer[thread_idx]) };

                    // call integral function
                    unsafe { self.integral_block(integrator, buf, &shls, &[], cint_opt, cache) };

                    // copy buffer to output slice
                    let out = unsafe { cast_mut_slice(&*out) };
                    let out_offsets = [0, cgto_loc_i, cgto_loc_j];
                    let buf_shape = [cgto_i, cgto_j, n_comp];
                    copy_c_3d_s2ij(out, &out_offsets, &out_shape, buf, &buf_shape);
                });
            },
            3 => {
                let out_shape = [n_comp, cgto_shape[0], cgto_shape[1]];

                let nidx_k = (shls_slice[K][1] - shls_slice[K][0]) as usize;
                let iter_layout = [nidx_ij, nidx_k].c();
                let iter_indices = IndexedIterLayout::new(&iter_layout, RowMajor).unwrap();

                iter_indices.into_par_iter().for_each(|([idx_ij, idx_k], _)| {
                    let [idx_j, idx_i] = unravel_s2_indices(idx_ij);

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
                    let cache = unsafe { cast_mut_slice(&thread_cache[thread_idx]) };
                    let buf = unsafe { cast_mut_slice(&thread_buffer[thread_idx]) };

                    // call integral function
                    unsafe { self.integral_block(integrator, buf, &shls, &[], cint_opt, cache) };

                    // copy buffer to output slice
                    let out = unsafe { cast_mut_slice(&*out) };
                    let out_offsets = [0, cgto_loc_i, cgto_loc_j, cgto_loc_k];
                    let buf_shape = [cgto_i, cgto_j, cgto_k, n_comp];
                    copy_c_4d_s2ij(out, &out_offsets, &out_shape, buf, &buf_shape);
                });
            },
            4 => {
                let out_shape = [n_comp, cgto_shape[0], cgto_shape[1], cgto_shape[2]];

                let nidx_k = (shls_slice[K][1] - shls_slice[K][0]) as usize;
                let nidx_l = (shls_slice[L][1] - shls_slice[L][0]) as usize;
                let iter_layout = [nidx_ij, nidx_k, nidx_l].c();
                let iter_indices = IndexedIterLayout::new(&iter_layout, RowMajor).unwrap();

                iter_indices.into_par_iter().for_each(|([idx_ij, idx_k, idx_l], _)| {
                    let [idx_j, idx_i] = unravel_s2_indices(idx_ij);

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
                    let cache = unsafe { cast_mut_slice(&thread_cache[thread_idx]) };
                    let buf = unsafe { cast_mut_slice(&thread_buffer[thread_idx]) };

                    // call integral function
                    unsafe { self.integral_block(integrator, buf, &shls, &[], cint_opt, cache) };

                    // copy buffer
                    let out = unsafe { cast_mut_slice(&*out) };
                    let out_offsets = [0, cgto_loc_i, cgto_loc_j, cgto_loc_k, cgto_loc_l];
                    let buf_shape = [cgto_i, cgto_j, cgto_k, cgto_l, n_comp];
                    copy_c_5d_s2ij(out, &out_offsets, &out_shape, buf, &buf_shape);
                });
            },
            _ => unreachable!(),
        }

        /* #endregion */

        Ok(())
    }

    #[doc(hidden)]
    pub fn integral_s2kl_row_major_inplace<F>(
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
        // n_center must be 4 for S2kl symmetry, checked in `check_shls_slice`
        let cgto_shape = self.cgto_shape_s2kl(shls_slice); // AO shape, without intor component
        let cgto_locs = self.cgto_locs(shls_slice); // AO relative locations mapped to shells, 0-indexed

        // cache (thread local)
        let cache_size = self.max_cache_size(integrator, shls_slice);
        let buffer_size = self.max_buffer_size(integrator, shls_slice);
        let thread_cache = (0..rayon::current_num_threads()).map(|_| unsafe { aligned_uninitialized_vec::<f64>(cache_size) }).collect_vec();
        let thread_buffer = (0..rayon::current_num_threads()).map(|_| unsafe { aligned_uninitialized_vec::<F>(buffer_size) }).collect_vec();

        /* #endregion */

        /* #region parallel integration generation */

        let out_shape = [n_comp, cgto_shape[0], cgto_shape[1], cgto_shape[2]];

        const I: usize = 0; // index of first shell
        const J: usize = 1; // index of second shell
        const K: usize = 2; // index of third shell
        const L: usize = 3; // index of fourth shell

        let nidx_i = (shls_slice[I][1] - shls_slice[I][0]) as usize;
        let nidx_j = (shls_slice[J][1] - shls_slice[J][0]) as usize;
        let nidx_k = (shls_slice[K][1] - shls_slice[K][0]) as usize;
        let nidx_kl = nidx_k * (nidx_k + 1) / 2; // number of unique (k, l) pairs
        let iter_layout = [nidx_i, nidx_j, nidx_kl].c();
        let iter_indices = IndexedIterLayout::new(&iter_layout, RowMajor).unwrap();

        iter_indices.into_par_iter().for_each(|([idx_i, idx_j, idx_kl], _)| {
            let [idx_l, idx_k] = unravel_s2_indices(idx_kl);

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
            let cache = unsafe { cast_mut_slice(&thread_cache[thread_idx]) };
            let buf = unsafe { cast_mut_slice(&thread_buffer[thread_idx]) };

            // call integral function
            unsafe { self.integral_block(integrator, buf, &shls, &[], cint_opt, cache) };

            // copy buffer
            let out = unsafe { cast_mut_slice(&*out) };
            let out_offsets = [0, cgto_loc_i, cgto_loc_j, cgto_loc_k, cgto_loc_l];
            let buf_shape = [cgto_i, cgto_j, cgto_k, cgto_l, n_comp];
            copy_c_5d_s2kl(out, &out_offsets, &out_shape, buf, &buf_shape);
        });

        /* #endregion */

        Ok(())
    }

    #[doc(hidden)]
    pub fn integral_s4_row_major_inplace<F>(
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
        // n_center must be 4 for S4 symmetry, checked in `check_shls_slice`
        let cgto_shape = self.cgto_shape_s4(shls_slice); // AO shape, without intor component
        let cgto_locs = self.cgto_locs(shls_slice); // AO relative locations mapped to shells, 0-indexed

        // cache (thread local)
        let cache_size = self.max_cache_size(integrator, shls_slice);
        let buffer_size = self.max_buffer_size(integrator, shls_slice);
        let thread_cache = (0..rayon::current_num_threads()).map(|_| unsafe { aligned_uninitialized_vec::<f64>(cache_size) }).collect_vec();
        let thread_buffer = (0..rayon::current_num_threads()).map(|_| unsafe { aligned_uninitialized_vec::<F>(buffer_size) }).collect_vec();

        /* #endregion */

        /* #region parallel integration generation */

        let out_shape = [n_comp, cgto_shape[0], cgto_shape[1]];

        const I: usize = 0; // index of first shell
        const J: usize = 1; // index of second shell
        const K: usize = 2; // index of third shell
        const L: usize = 3; // index of fourth shell

        let nidx_i = (shls_slice[I][1] - shls_slice[I][0]) as usize;
        let nidx_k = (shls_slice[K][1] - shls_slice[K][0]) as usize;
        let nidx_ij = nidx_i * (nidx_i + 1) / 2; // number of unique (i, j) pairs
        let nidx_kl = nidx_k * (nidx_k + 1) / 2; // number of unique (k, l) pairs
        let iter_layout = [nidx_ij, nidx_kl].c();
        let iter_indices = IndexedIterLayout::new(&iter_layout, RowMajor).unwrap();

        iter_indices.into_par_iter().for_each(|([idx_ij, idx_kl], _)| {
            let [idx_j, idx_i] = unravel_s2_indices(idx_ij);
            let [idx_l, idx_k] = unravel_s2_indices(idx_kl);

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
            let cache = unsafe { cast_mut_slice(&thread_cache[thread_idx]) };
            let buf = unsafe { cast_mut_slice(&thread_buffer[thread_idx]) };

            // call integral function
            unsafe { self.integral_block(integrator, buf, &shls, &[], cint_opt, cache) };

            // copy buffer
            let out = unsafe { cast_mut_slice(&*out) };
            let out_offsets = [0, cgto_loc_i, cgto_loc_j, cgto_loc_k, cgto_loc_l];
            let buf_shape = [cgto_i, cgto_j, cgto_k, cgto_l, n_comp];
            copy_c_5d_s4(out, &out_offsets, &out_shape, buf, &buf_shape);
        });

        /* #endregion */

        Ok(())
    }

    #[doc(hidden)]
    pub fn integral_s8_row_major_inplace<F>(
        &self,
        integrator: &dyn Integrator,
        out: &mut [F],
        shls_slice: &[[c_int; 2]],
        cint_opt: Option<&CIntOptimizer>,
    ) -> Result<(), CIntError>
    where
        F: ComplexFloat + Send + Sync,
    {
        // Row-major s8 is actually the same to col-major s8.
        //
        // Though it is very unlikely to happen, if component dimension occurs,
        // then simply transpose col-major output; otherwise, it is only one
        // dimension and row/col-major is the same.
        //
        // And since inplace-function only writes out buffer, so simply call col-major
        // should work.

        self.integral_s8_inplace(integrator, out, shls_slice, cint_opt)
    }
}
