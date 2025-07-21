//! Properties of the [`CInt`] instance (compute and memory cost is not
//! essential in this module).

use crate::prelude::*;

/// Properties of the `CInt` instance (compute and memory cost is not essential
/// in this module).
impl CInt {
    /// Whether pseudo potential is used in the system.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.has_ecp`
    ///
    /// # Examples
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// assert!(!cint_data.has_ecp());
    ///
    /// let cint_data = init_sb2me4_cc_pvtz();
    /// assert!(cint_data.has_ecp());
    /// ```
    #[inline]
    pub fn has_ecp(&self) -> bool {
        !self.ecpbas.is_empty()
    }

    /// Whether spin-orbit coupling is enabled in ECP.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.has_ecp_soc`
    ///
    /// # Examples
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_sb2me4_cc_pvtz();
    /// assert!(!cint_data.has_ecp_soc());
    /// ```
    #[inline]
    pub fn has_ecp_soc(&self) -> bool {
        const SO_TYPE_OF: usize = cecp_ffi::SO_TYPE_OF as usize;
        self.ecpbas.iter().any(|ecp| ecp[SO_TYPE_OF] == 1)
    }

    /// Number of shells in the system.
    ///
    /// This does not count ECP shells.
    ///
    /// Please note that, `self.bas.len()` may not be the actual number of
    /// shells, because in some cases, ECP shells are merged into the GTO
    /// shells.
    ///
    /// # PySCF Equivalent
    ///
    /// Attribute `Mole.nbas`
    ///
    /// # Examples
    ///
    /// ```rust
    /// use libcint::prelude::*;
    ///
    /// let cint_data = init_h2o_def2_tzvp();
    /// assert_eq!(cint_data.nbas(), 19);
    ///
    /// let cint_data = init_sb2me4_cc_pvtz();
    /// assert_eq!(cint_data.nbas(), 130);
    ///
    /// # // This is test that when ecpbas is merged, the number of shells is still correct.
    /// # let merged = cint_data.merge_ecpbas();
    /// # assert_eq!(merged.nbas(), 130);
    /// ```
    #[inline]
    pub fn nbas(&self) -> usize {
        if self.is_ecp_merged() { self.bas.len() - self.ecpbas.len() } else { self.bas.len() }
    }

    /// Number of atoms in the system.
    #[inline]
    pub fn natm(&self) -> usize {
        self.atm.len()
    }

    /// Pointer to the shell data.
    #[inline]
    pub fn bas_ptr(&self) -> *const c_int {
        self.bas.as_ptr() as *const c_int
    }

    /// Pointer to the atom data.
    #[inline]
    pub fn atm_ptr(&self) -> *const c_int {
        self.atm.as_ptr() as *const c_int
    }

    /// Pointer to the environment data.
    #[inline]
    pub fn env_ptr(&self) -> *const f64 {
        self.env.as_ptr()
    }

    /// Nuclear effective charge of the given atom id.
    ///
    /// # Note
    ///
    /// `atom_charge != charge(atom_symbol)` when ECP is enabled.
    ///
    /// Number of electrons screened by ECP can be obtained by
    /// `charge(atom_symbol) - atom_charge`.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.atom_charge`
    ///
    /// # Examples
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// assert_eq!(cint_data.atom_charge(0), 8.0);  // O  (8)
    /// assert_eq!(cint_data.atom_charge(1), 1.0);  // H  (1)
    ///
    /// let cint_data = init_sb2me4_cc_pvtz();
    /// assert_eq!(cint_data.atom_charge(0), 23.0); // Sb (51) with ECP (-36)
    /// assert_eq!(cint_data.atom_charge(2), 6.0);  // C  (6)
    /// ```
    #[inline]
    pub fn atom_charge(&self, atm_id: usize) -> f64 {
        const NUC_MOD_OF: usize = cint_ffi::NUC_MOD_OF as usize;
        const CHARGE_OF: usize = cint_ffi::CHARGE_OF as usize;
        const FRAC_CHARGE_NUC: u32 = cint_ffi::FRAC_CHARGE_NUC;
        const PTR_FRAC_CHARGE: usize = cint_ffi::PTR_FRAC_CHARGE as usize;

        if self.atm[atm_id][NUC_MOD_OF] as u32 != FRAC_CHARGE_NUC {
            // regular QM atoms
            self.atm[atm_id][CHARGE_OF] as f64
        } else {
            // MM atoms with fractional charges
            self.env[self.atm[atm_id][PTR_FRAC_CHARGE] as usize]
        }
    }

    /// List of Nuclear effective charge of all atoms in system.
    ///
    /// # Note
    ///
    /// `atom_charge != charge(atom_symbol)` when ECP is enabled.
    ///
    /// Number of electrons screened by ECP can be obtained by
    /// `charge(atom_symbol) - atom_charge`.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.atom_charges`
    ///
    /// # Examples
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// assert_eq!(cint_data.atom_charges(), vec![8., 1., 1.]);
    ///
    /// // For Sb2Me4, the first atom is Sb (51) with ECP (-36), so the charge is 23.
    /// let cint_data = init_sb2me4_cc_pvtz();
    /// assert_eq!(cint_data.atom_charges(), vec![23., 23., 6., 6., 6., 6., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1., 1.]);
    /// ```
    #[inline]
    pub fn atom_charges(&self) -> Vec<f64> {
        (0..self.atm.len()).map(|i| self.atom_charge(i)).collect()
    }

    /// The number of spinor associated with given angular momentum $l$ and
    /// kappa $\kappa$.
    ///
    /// - If $\kappa = 0$, it returns $4 l + 2$.
    /// - If $\kappa > 0$, it returns $2 l + 2$.
    /// - If $\kappa < 0$, it returns $2 l$.
    ///
    /// # PySCF Equivalent
    ///
    /// Function `gto.mole.len_spinor`
    #[inline]
    pub fn len_spinor(l: i32, kappa: i32) -> usize {
        if kappa == 0 {
            (l * 4 + 2) as usize
        } else if kappa > 0 {
            (l * 2 + 2) as usize
        } else {
            (l * 2) as usize
        }
    }

    /// The number of Cartesian function associated with given angular momentum
    /// $l$.
    ///
    /// This will gives $\frac{(l + 1) (l + 2)}{2}$
    ///
    /// # PySCF Equivalent
    ///
    /// Function `gto.mole.len_cart`
    #[inline]
    pub fn len_cart(l: i32) -> usize {
        ((l + 1) * (l + 2) / 2) as usize
    }

    /// The number of spherical function associated with given angular momentum
    /// $l$.
    ///
    /// This will gives $2 l + 1$.
    ///
    /// # PySCF Equivalent
    ///
    /// Function `gto.mole.len_sph`
    #[inline]
    pub fn len_sph(l: i32) -> usize {
        (l * 2 + 1) as usize
    }

    /// Location mapping from shell to basis.
    ///
    /// The type of integral is specified by struct field `cint_type`.
    ///
    /// Output vector is of length `nshl + 1`, where `nshl` is the number of
    /// shells.
    ///
    /// # PySCF Equivalent
    ///
    /// This implementation follows `gto.moleintor.make_loc`.
    ///
    /// For `gto.Mole` object, methods `ao_loc`, `ao_loc_nr`, `ao_loc_2c` are
    /// also relavent.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    ///
    /// let loc_sph = cint_data.make_loc();
    /// assert_eq!(loc_sph, vec![ 0,  1,  2,  3,  4,  5,  8, 11, 14, 19, 24, 31, 32, 33, 34, 37, 38, 39, 40, 43]);
    /// ```
    #[inline]
    pub fn make_loc(&self) -> Vec<usize> {
        self.make_loc_with_type(self.cint_type)
    }

    /// Location mapping from shell to basis.
    ///
    /// # See also
    ///
    /// [`CInt::make_loc`]
    #[inline]
    pub fn ao_loc(&self) -> Vec<usize> {
        self.make_loc()
    }

    /// Location mapping from shell to basis (with integral type specified).
    ///
    /// This mapping is cumulated, and can be different for sph, cart and spinor
    /// types.
    ///
    /// Output vector is of length `nshl + 1`, where `nshl` is the number of
    /// shells.
    ///
    /// # PySCF Equivalent
    ///
    /// This implementation follows `gto.moleintor.make_loc`.
    ///
    /// For `gto.Mole` object, methods `ao_loc`, `ao_loc_nr`, `ao_loc_2c` are
    /// also relavent.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    ///
    /// let loc_sph = cint_data.make_loc_with_type(CIntType::Spheric);
    /// assert_eq!(loc_sph, vec![ 0,  1,  2,  3,  4,  5,  8, 11, 14, 19, 24, 31, 32, 33, 34, 37, 38, 39, 40, 43]);
    ///
    /// let loc_cart = cint_data.make_loc_with_type(CIntType::Cartesian);
    /// assert_eq!(loc_cart, vec![ 0,  1,  2,  3,  4,  5,  8, 11, 14, 20, 26, 36, 37, 38, 39, 42, 43, 44, 45, 48]);
    ///
    /// let loc_spinor = cint_data.make_loc_with_type(CIntType::Spinor);
    /// assert_eq!(loc_spinor, vec![ 0,  2,  4,  6,  8, 10, 16, 22, 28, 38, 48, 62, 64, 66, 68, 74, 76, 78, 80, 86]);
    /// ```
    pub fn make_loc_with_type(&self, cint_type: CIntType) -> Vec<usize> {
        const ANG_OF: usize = cint_ffi::ANG_OF as usize;
        const KAPPA_OF: usize = cint_ffi::KAPPA_OF as usize;
        const NCTR_OF: usize = cint_ffi::NCTR_OF as usize;

        let mut ao_loc = vec![0];
        let nbas = self.nbas();
        for shl in 0..nbas {
            let l = self.bas[shl][ANG_OF];
            let k = self.bas[shl][KAPPA_OF];
            let nctr = self.bas[shl][NCTR_OF] as usize;
            let val = match cint_type {
                Spheric => Self::len_sph(l) * nctr,
                Cartesian => Self::len_cart(l) * nctr,
                Spinor => Self::len_spinor(l, k) * nctr,
            };
            ao_loc.push(ao_loc[shl] + val);
        }
        ao_loc
    }

    /// Get the number of basis (atomic orbitals).
    ///
    /// The type of integral is specified by struct field `cint_type`.
    ///
    /// This value is the same to the last of [`CInt::make_loc_with_type`].
    ///
    /// # PySCF Equivalent
    ///
    /// For `gto.Mole` object, methods `nao_nr`, `nao_cart`, `nao_2c` are
    /// relavent.
    ///
    ///
    /// # Examples
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    ///
    /// let nao_sph = cint_data.nao();
    /// assert_eq!(nao_sph, 43);
    ///
    /// let cint_data = init_sb2me4_cc_pvtz();
    /// let nao_sph = cint_data.nao();
    /// assert_eq!(nao_sph, 366);
    ///
    /// # // This is test that when ecpbas is merged, the number of
    /// # // atomic orbitals is still correct.
    /// # let merged = cint_data.merge_ecpbas();
    /// # let nao_cart = merged.nao();
    /// # assert_eq!(nao_cart, 366);
    /// ```
    #[inline]
    pub fn nao(&self) -> usize {
        self.nao_with_type(self.cint_type)
    }

    /// Get the number of basis (atomic orbitals, with integral type specified).
    ///
    /// This value is different for sph, cart and spinor types.
    ///
    /// This value is the same to the last of [`CInt::make_loc_with_type`].
    ///
    /// # PySCF Equivalent
    ///
    /// For `gto.Mole` object, methods `nao_nr`, `nao_cart`, `nao_2c` are
    /// relavent.
    ///
    ///
    /// # Examples
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    ///
    /// let nao_sph = cint_data.nao_with_type(CIntType::Spheric);
    /// assert_eq!(nao_sph, 43);
    ///
    /// let nao_cart = cint_data.nao_with_type(CIntType::Cartesian);
    /// assert_eq!(nao_cart, 48);
    ///
    /// let nao_spinor = cint_data.nao_with_type(CIntType::Spinor);
    /// assert_eq!(nao_spinor, 86);
    ///
    /// let cint_data = init_sb2me4_cc_pvtz();
    /// let nao_sph = cint_data.nao_with_type(CIntType::Spheric);
    /// assert_eq!(nao_sph, 366);
    ///
    /// # // This is test that when ecpbas is merged, the number of
    /// # // atomic orbitals is still correct.
    /// # let merged = cint_data.merge_ecpbas();
    /// # let nao_cart = merged.nao_with_type(CIntType::Spheric);
    /// # assert_eq!(nao_cart, 366);
    /// ```
    #[inline]
    pub fn nao_with_type(&self, cint_type: CIntType) -> usize {
        const ANG_OF: usize = cint_ffi::ANG_OF as usize;
        const KAPPA_OF: usize = cint_ffi::KAPPA_OF as usize;
        const NCTR_OF: usize = cint_ffi::NCTR_OF as usize;

        let mut nao = 0;
        let nbas = self.nbas();
        for shl in 0..nbas {
            let l = self.bas[shl][ANG_OF];
            let k = self.bas[shl][KAPPA_OF];
            let nctr = self.bas[shl][NCTR_OF];
            let val = match cint_type {
                Spheric => Self::len_sph(l) as c_int * nctr,
                Cartesian => Self::len_cart(l) as c_int * nctr,
                Spinor => Self::len_spinor(l, k) as c_int * nctr,
            };
            nao += val as usize;
        }
        nao
    }

    /// Coordinates fo the given atom id in unit Bohr (a.u.).
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.atom_coord` with `unit="Bohr"`.
    #[inline]
    pub fn atom_coord(&self, atm_id: usize) -> [f64; 3] {
        const PTR_COORD: usize = cint_ffi::PTR_COORD as usize;

        let ptr = self.atm[atm_id][PTR_COORD] as usize;
        [self.env[ptr], self.env[ptr + 1], self.env[ptr + 2]]
    }

    /// Coordinates of all atoms in unit Bohr (a.u.).
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.atom_coords` with `unit="Bohr"`.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// # use approx::assert_relative_eq;
    /// let cint_data = init_h2o_def2_tzvp();
    /// let coords = cint_data.atom_coords();
    /// // Output:
    /// // [[  0.000000,   0.000000,   0.000000],
    /// //  [  1.776343,   0.000000,   0.000000],
    /// //  [ -0.444761,   0.000000,   1.719762]]
    /// assert_relative_eq!(cint_fingerprint(&coords), -2.4358371781626658, max_relative=1e-12);
    /// ```
    #[inline]
    pub fn atom_coords(&self) -> Vec<[f64; 3]> {
        (0..self.atm.len()).map(|i| self.atom_coord(i)).collect_vec()
    }

    /// Number of shells of the given atom
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.atom_nshells`
    ///
    /// # Examples
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// assert_eq!(cint_data.atom_nshells(0), 11);
    /// ```
    #[inline]
    pub fn atom_nshells(&self, atm_id: usize) -> usize {
        const ATOM_OF: usize = cint_ffi::ATOM_OF as usize;

        let nbas = self.nbas();
        self.bas[..nbas].iter().filter(|bas| bas[ATOM_OF] as usize == atm_id).count()
    }

    /// List of shell ids of the given atom.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.atom_shell_ids`
    ///
    /// # Examples
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// assert_eq!(cint_data.atom_shell_ids(0), vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]);
    /// ```
    #[inline]
    pub fn atom_shell_ids(&self, atm_id: usize) -> Vec<usize> {
        const ATOM_OF: usize = cint_ffi::ATOM_OF as usize;

        let nbas = self.nbas();
        self.bas[..nbas].iter().enumerate().filter_map(|(i, bas)| if bas[ATOM_OF] as usize == atm_id { Some(i) } else { None }).collect()
    }

    #[inline]
    fn check_bas_id(&self, bas_id: usize) {
        let nbas = self.nbas();
        if bas_id >= nbas {
            panic!("Invalid bas_id: {bas_id}. Number of shells is {nbas}");
        }
    }

    /// Coordinates of the given shell in unit Bohr (a.u.).
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.bas_coord`
    #[inline]
    pub fn bas_coord(&self, bas_id: usize) -> [f64; 3] {
        const PTR_COORD: usize = cint_ffi::PTR_COORD as usize;

        self.check_bas_id(bas_id);
        let atm_id = self.bas_atom(bas_id);
        let ptr = self.atm[atm_id][PTR_COORD] as usize;
        [self.env[ptr], self.env[ptr + 1], self.env[ptr + 2]]
    }

    /// The atom (0-based id) that the given shell sits on.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.bas_atom`
    #[inline]
    pub fn bas_atom(&self, bas_id: usize) -> usize {
        const ATOM_OF: usize = cint_ffi::ATOM_OF as usize;

        self.check_bas_id(bas_id);
        self.bas[bas_id][ATOM_OF] as usize
    }

    /// The angular momentum associated with the given shell.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.bas_angular`
    #[inline]
    pub fn bas_angular(&self, bas_id: usize) -> usize {
        const ANG_OF: usize = cint_ffi::ANG_OF as usize;

        self.check_bas_id(bas_id);
        self.bas[bas_id][ANG_OF] as usize
    }

    /// The number of contracted GTOs for the given shell.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.bas_nctr`
    #[inline]
    pub fn bas_nctr(&self, bas_id: usize) -> usize {
        const NCTR_OF: usize = cint_ffi::NCTR_OF as usize;

        self.check_bas_id(bas_id);
        self.bas[bas_id][NCTR_OF] as usize
    }

    /// The number of primitive GTOs for the given shell.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.bas_nprim`
    #[inline]
    pub fn bas_nprim(&self, bas_id: usize) -> usize {
        const NPRIM_OF: usize = cint_ffi::NPRIM_OF as usize;

        self.check_bas_id(bas_id);
        self.bas[bas_id][NPRIM_OF] as usize
    }

    /// Kappa (if l < j, -l-1, else l) of the given shell.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.bas_kappa`
    #[inline]
    pub fn bas_kappa(&self, bas_id: usize) -> c_int {
        const KAPPA_OF: usize = cint_ffi::KAPPA_OF as usize;

        self.check_bas_id(bas_id);
        self.bas[bas_id][KAPPA_OF]
    }

    /// Exponents of the given shell.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.bas_exp`
    #[inline]
    pub fn bas_exp(&self, bas_id: usize) -> &[f64] {
        const PTR_EXP: usize = cint_ffi::PTR_EXP as usize;

        self.check_bas_id(bas_id);
        let ptr = self.bas[bas_id][PTR_EXP] as usize;
        &self.env[ptr..ptr + self.bas_nprim(bas_id)]
    }

    /// Exponents of all shells.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.bas_exps`
    #[inline]
    pub fn bas_exps(&self) -> Vec<Vec<f64>> {
        (0..self.nbas()).map(|i| self.bas_exp(i).to_vec()).collect()
    }

    /// Shell and basis (atomic orbitals) offsets for each atom.
    ///
    /// This will give a list of `[shl_start, shl_end, ao_start, ao_end]` for
    /// each atom.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.aoslice_by_atom`. This implementation does not allow
    /// aribitary `ao_loc` (which can be obtained by [`CInt::make_loc`]).
    /// However, user can provide `cint_type` to specify if you wish to use
    /// sph, cart or spinor.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// let aoslice = cint_data.aoslice_by_atom();
    /// assert_eq!(aoslice, vec![
    ///     [ 0, 11,  0, 31],
    ///     [11, 15, 31, 37],
    ///     [15, 19, 37, 43],
    /// ]);
    /// ```
    #[inline]
    pub fn aoslice_by_atom(&self) -> Vec<[usize; 4]> {
        self.aoslice_by_atom_with_type(self.cint_type)
    }

    /// Shell and basis (atomic orbitals) offsets for each atom (with integral
    /// type specified).
    ///
    /// This will give a list of `[shl_start, shl_end, ao_start, ao_end]` for
    /// each atom.
    ///
    /// # PySCF Equivalent
    ///
    /// Method `Mole.aoslice_by_atom`. This implementation does not allow
    /// aribitary `ao_loc` (which can be obtained by [`CInt::make_loc`]).
    /// However, user can provide `cint_type` to specify if you wish to use
    /// sph, cart or spinor.
    ///
    /// # Examples
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// let aoslice = cint_data.aoslice_by_atom_with_type(CIntType::Spheric);
    /// assert_eq!(aoslice, vec![
    ///     [ 0, 11,  0, 31],
    ///     [11, 15, 31, 37],
    ///     [15, 19, 37, 43],
    /// ]);
    /// ```
    pub fn aoslice_by_atom_with_type(&self, cint_type: CIntType) -> Vec<[usize; 4]> {
        const ATOM_OF: usize = cint_ffi::ATOM_OF as usize;

        let nbas = self.nbas();
        let ao_loc = self.make_loc_with_type(cint_type);

        // category atoms with basis index
        // [atm_id, shl_start, shl_end]
        let mut category = vec![];
        if nbas > 0 {
            let mut current_atm = self.bas[0][ATOM_OF];
            let mut shl_start = 0;
            for shl in 0..self.nbas() {
                let atm_id = self.bas[shl][ATOM_OF];
                if atm_id != current_atm {
                    category.push([current_atm as usize, shl_start, shl]);
                    current_atm = atm_id;
                    shl_start = shl;
                }
            }
            category.push([current_atm as usize, shl_start, nbas]);
        }

        // fill result based on category
        let mut result = vec![];
        let mut ptr_category = 0;
        for atm_id in 0..self.atm.len() {
            if ptr_category < category.len() && atm_id == category[ptr_category][0] {
                // atom exists in bas
                let shl_start = category[ptr_category][1];
                let shl_end = category[ptr_category][2];
                let ao_start = ao_loc[shl_start];
                let ao_end = ao_loc[shl_end];
                result.push([shl_start, shl_end, ao_start, ao_end]);
                ptr_category += 1;
            } else {
                // atom missing in bas
                // following unwrap occurs for first atom
                let &[_, shl_end, _, ao_end] = result.last().unwrap_or(&[0; 4]);
                result.push([shl_end, shl_end, ao_end, ao_end]);
            }
        }

        // panic if I am wrong
        // check result length
        assert_eq!(result.len(), self.atm.len());
        // check total number of basis
        assert_eq!(result.last().unwrap()[3], self.nao_with_type(cint_type));
        result
    }

    /// Number of grids that is used in `int1e_grids`.
    ///
    /// Note that this function is mostly for internal usage. User should
    /// generally not call this function.
    pub fn ngrids(&self) -> usize {
        const NGRIDS: usize = crate::ffi::cint_ffi::NGRIDS as usize;
        self.env[NGRIDS] as usize
    }

    /// Balances the partition of shells into blocks of a given size.
    ///
    /// This function can be useful if the full bulk of integrals is not
    /// available in memory, and you have to generate them batch by batch.
    ///
    /// The outputs are
    /// - `shl_start`: the starting index of the shell in the partition
    ///   (included)
    /// - `shl_end`: the ending index of the shell in the partition (not
    ///   included)
    /// - `nbatch_ao`: the number of AOs in the batch.
    ///
    /// # PySCF equivalent
    ///
    /// `ao2mo.outcore.balance_partition`
    ///
    /// # Example
    ///
    /// ```
    /// use libcint::prelude::*;
    /// let cint_data = init_h2o_def2_tzvp();
    /// let partition = cint_data.balance_partition(20);
    /// assert_eq!(partition, vec![[0, 9, 19], [9, 17, 20], [17, 19, 4]]);
    /// ```
    pub fn balance_partition(&self, block_size: usize) -> Vec<[usize; 3]> {
        crate::util::balance_partition(self, block_size, None, None)
    }

    /// Balances the partition of shells into blocks of a given size.
    ///
    /// This function is similar to [`CInt::balance_partition`], but allows
    /// user to specify the start and end id of the shells.
    pub fn balance_partition_advanced(&self, block_size: usize, start_id: Option<usize>, end_id: Option<usize>) -> Vec<[usize; 3]> {
        crate::util::balance_partition(self, block_size, start_id, end_id)
    }
}

#[cfg(test)]
mod test {
    use crate::prelude::*;

    #[test]
    fn playground() {
        let cint_data = init_h2o_def2_tzvp();
        println!("{:?}", cint_data.aoslice_by_atom());
    }
}
