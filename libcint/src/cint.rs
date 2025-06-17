//! Main [`CInt`] with related struct definition, and impl of
//! [`CInt::integrate`].

use crate::prelude::*;

use crate::ffi::cecp_wrapper::get_ecp_integrator;
use crate::ffi::cint_wrapper::get_cint_integrator;

/* #region structs */

/// CInt data structure, which contains **all the necessary information** for
/// GTO electronic integral evaluation, and almost **all methods to evaluate**
/// GTO electronic integrals for basic usage.
///
/// Most users may wish to use the [`integrate`](Self::integrate) and
/// [`integrate_spinor`](Self::integrate_spinor) methods, which correspond to
/// PySCF's `mol.intor` method. We refer to those functions for more
/// documentations.
///
/// Documentation and code conventions:
///
/// | code name | description |
/// |------------|-------------|
/// | `atm` | atom |
/// | `bas` | shell |
/// | `shl` | shell |
/// | `ao` | atomic orbital basis (for users, absolute where starting shell is alway the first shell from `bas`) |
/// | `cgto` | atomic orbital basis (for developers, relative to specified starting shell by `shls_slice`) |
/// | `cint`, `CInt` | libcint's [`CInt`] structure or instances |
/// | `c_int` | [`std::ffi::c_int`] type ([`i32`] in most cases) |
#[derive(Debug, Clone, PartialEq)]
pub struct CInt {
    /// Slots of atoms.
    ///
    /// Names of this field can be retrived in [`cint_ffi`].
    ///
    /// | Index | Name          | Description | Related getter | Related setter |
    /// |-------|---------------|-------------|----------------|----------------|
    /// | 0 | `CHARGE_OF`       | atomic charge (ECP core electrons excluded)       | [`atom_charge`](Self::atom_charge) <br/> [`atom_charges`](Self::atom_charges) |
    /// | 1 | `PTR_COORD`       | coordinates location (pointer to `env` field)     | [`atom_coord`](Self::atom_coord) <br/> [`atom_coords`](Self::atom_coords) | [`set_geom`](Self::set_geom) |
    /// | 2 | `NUC_MOD_OF`      | nuclear mode of the atom                          | | [`set_nuc_mod`](Self::set_nuc_mod) |
    /// | 3 | `PTR_ZETA`        | pointer to zeta values (pointer to `env` field)   | | [`set_rinv_at_nucleus`](Self::set_rinv_at_nucleus) |
    /// | 4 | `PTR_FRAC_CHARGE` | fractional charge (pointer to `env` field)        | [`atom_charge`](Self::atom_charge) <br/> [`atom_charges`](Self::atom_charges) |
    /// | 5 | `RESERVE_ATMSLOT` |
    /// | 6 | `ATM_SLOTS`       |
    ///
    /// - 2 - `NUC_MOD_OF` can be three kinds:
    ///     - 1 - `POINT_NUC` point charge, which is the most common case.
    ///     - 2 - `GAUSSIAN_NUC` nuclear charge is represented by a gaussian
    ///       distribution, with a zeta value (related to `PTR_ZETA`).
    ///     - 3 - `MM_NUC` is a point charge with fractional charge, which is
    ///       used in molecular mechanics (MM) calculations, with a fractional
    ///       value (related to `PTR_FRAC_CHARGE`).
    /// - 3 - Zeta value is atomic charge that is gaussian-distributed.
    ///   (`NUC_MOD_OF` refers to `GAUSSIAN_NUC`). In other cases, this value
    ///   should be zero.
    pub atm: Vec<[c_int; 6]>,

    /// Slot of shells (as minimal building block of contracted GTO).
    ///
    /// | Index | Name          | Description | Related getter |
    /// |-------|---------------|-------------|----------------|
    /// | 0 | `ATOM_OF`        | 0-based index of corresponding atom                            |
    /// | 1 | `ANG_OF`         | angular momentum                                               | [`bas_angular`](Self::bas_angular) |
    /// | 2 | `NPRIM_OF`       | number of primitive GTO in basis                               | [`bas_nprim`](Self::bas_nprim) |
    /// | 3 | `NCTR_OF`        | number of contracted GTO in basis                              | [`bas_nctr`](Self::bas_nctr) |
    /// | 4 | `KAPPA_OF`       | kappa for spinor GTO                                           | [`bas_kappa`](Self::bas_kappa) |
    /// | 5 | `PTR_EXP`        | exponents of primitive GTOs (pointer to `env` field)           |
    /// | 6 | `PTR_COEFF`      | column-major contraction coefficients (pointer to `env` field) |
    /// | 7 | `RESERVE_BASLOT` |
    /// | 8 | `BAS_SLOTS`      |
    pub bas: Vec<[c_int; 8]>,

    /// Slot of shells for ECP (effective core potential).
    ///
    /// Other entries are the same as `bas` field, but with some additional:
    ///
    /// | Index | Name         | Description |
    /// |-------|--------------|-------------|
    /// | 1 | `ANG_OF`         | angular momentum (refers to angular of $U_L$) |
    /// | 3 | `RADI_POWER`     | number of maximum radial power to be summed   |
    /// | 4 | `SO_TYPE_OF`     | spin-orb type                                 |
    pub ecpbas: Vec<[c_int; 8]>,

    /// Slot of floating-point variables.
    ///
    /// | Index | Name          | Description | Related getter | Related setter |
    /// |-------|---------------|-------------|----------------|----------------|
    /// |  0 | `PTR_EXPCUTOFF`     | Overall cutoff for integral prescreening, value needs to be ln (threshold) | [`get_integral_screen`](Self::get_integral_screen) | [`set_integral_screen`](Self::set_integral_screen) |
    /// |  1 | `PTR_COMMON_ORIG`   | $R_C$ (common origin) of $r-R_C$ in dipole, GIAO operators                 | [`get_common_origin`](Self::get_common_origin) | [`set_common_origin`](Self::set_common_origin) |
    /// |  4 | `PTR_RINV_ORIG`     | $R_O$ (rinv origin) of $1/(r-R_O)$                                         | [`get_rinv_origin`](Self::get_rinv_origin) | [`set_rinv_origin`](Self::set_rinv_origin) |
    /// |  7 | `PTR_RINV_ZETA`     | ZETA parameter for Gaussian charge distribution (Gaussian nuclear model)   | [`get_rinv_zeta`](Self::get_rinv_zeta) | [`set_rinv_zeta`](Self::set_rinv_zeta) |
    /// |  8 | `PTR_RANGE_OMEGA`   | omega parameter in range-separated coulomb operator (>0 for LR, <0 for SR) | [`get_range_coulomb`](Self::get_range_coulomb) <br/> [`omega`](`Self::omega`) | [`set_range_coulomb`](Self::set_range_coulomb) <br/> [`set_omega`](Self::set_omega) |
    /// |  9 | `PTR_F12_ZETA`      | Yukawa potential and Slater-type geminal $e^{-\zeta r}$                    |
    /// | 10 | `PTR_GTG_ZETA`      | Gaussian type geminal $e^{-\zeta r^2}$                                     |
    /// | 11 | `NGRIDS`            | For hybrid integrals with grids                                            |
    /// | 12 | `PTR_GRIDS`         | Location of grids for `int1e_grids` or variants                            |
    /// | 17 | `AS_RINV_ORIG_ATOM` | Position of atom to be used as origin in $1/(r-R_O)$ for ECP derivatives   | [`get_rinv_origin_atom`](Self::get_rinv_origin_atom) <br/> [`set_rinv_origin_atom`](Self::set_rinv_origin_atom) |
    /// | 18 | `AS_ECPBAS_OFFSET`  | Offset of ECP basis in `bas` field, used for ECP integrals                 | | [`merge_ecpbas`](Self::merge_ecpbas) <br/> [`decopule_ecpbas`](Self::decopule_ecpbas) |
    /// | 19 | `AS_NECPBAS`        | Number of ECP shells in `ecpbas` field                                     | | [`merge_ecpbas`](Self::merge_ecpbas) <br/> [`decopule_ecpbas`](Self::decopule_ecpbas) |
    /// | 20 | `PTR_ENV_START`     | Start of data                                                              |
    pub env: Vec<f64>,

    /// Type of integral.
    ///
    /// - `Spheric`: spherical harmonics, which is the most common case.
    /// - `Cartesian`: cartesian coordinates.
    /// - `Spinor`: spinor integrals, which are used in relativistic
    ///   calculations.
    pub cint_type: CIntType,
}

/// Error for [`CInt`] operations.
#[derive(Debug, Clone)]
pub enum CIntError {
    IntegratorNotFound(String),
    IntegratorNotAvailable(String),
    RuntimeError(String),
    InvalidValue(String),
    Miscellaneous(String),
    UninitializedFieldError(UninitializedFieldError),
}

/// Distinguish between general integral or ECP (effective core potential)
/// integral.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CIntKind {
    /// General integral type, which is used for most of the integrals.
    Int,
    /// ECP integral type, which is used for effective core potential integrals.
    Ecp,
}

/// Integral types.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CIntType {
    #[default]
    /// Spherical integrals, output double float `f64`.
    Spheric,
    /// Cartesian integrals, output double float `f64`.
    Cartesian,
    /// Spinor integrals, output `Complex<f64>`, and should be called by
    /// functions with `_spinor` suffix.
    Spinor,
}

/// Symmetry of the integral.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum CIntSymm {
    #[default]
    /// No symmetry.
    S1,
    /// Symmetric for the first two indices (in F-contiguous).
    S2ij,
    /// Symmetric for the last two indices in 4-center integrals (in
    /// F-contiguous).
    S2kl,
    /// Symmetric for the first and last indices in 4-center integrals (in
    /// F-contiguous).
    S4,
    /// Symmetric for all indices in 4-center integrals (in F-contiguous).
    S8,
}

/// Optimizer to be used in integral evaluation (not intended for user).
#[derive(Debug, PartialEq)]
pub enum CIntOptimizer {
    Int(*mut CINTOpt),
    Ecp(*mut ECPOpt),
}

/// Output of the integral evaluation (output buffer with col-major shape).
///
/// You can use `.into()` to convert it to a tuple of `(Vec<F>, Vec<usize>)`,
/// where `F` is `f64` for sph and cart integrals, and `Complex<f64>` for spinor
/// integrals.
pub struct CIntOutput<F> {
    /// Output buffer.
    ///
    /// It can be `None` when user provided the output buffer for integral: the
    /// integral will be written directly to the buffer user provided, and this
    /// crate will not allocate an output buffer. Note that in this case,
    /// `.into()` will panic.
    pub out: Option<Vec<F>>,

    /// Column-major shape of the output buffer. Please see [`CInt::integrate`]
    /// for more information.
    pub shape: Vec<usize>,
}

/// Integrate arguments for function [`CInt::integrate_with_args`] and
/// [`CInt::integrate_with_args_spinor`].
///
/// Builder of this struct can be retrived by [`CInt::integrate_args_builder`]
/// and [`CInt::integrate_args_builder_spinor`].
#[derive(Builder, Debug)]
#[builder(pattern = "owned", build_fn(error = "CIntError"))]
pub struct IntegrateArgs<'l, F> {
    /// Integrator name, such as "int1e_ovlp", "int2e", etc.
    #[builder(setter(into))]
    pub intor: &'l str,

    /// Shell slices for evaluating sub-tensor of the total integral.
    #[builder(default = "&[]")]
    pub shls_slice: &'l [[usize; 2]],

    /// Symmetry of the integral. This will affect shape of output tensor.
    #[builder(default, setter(into))]
    pub aosym: CIntSymm,

    /// Output buffer for the integral evaluation.
    ///
    /// If this is `None`, the integral engine will generate a buffer for
    /// output.
    #[builder(default, setter(strip_option))]
    pub out: Option<&'l mut [F]>,

    /// Col-major or row-major output buffer (col-major default).
    ///
    /// - Col-major: the same **data layout** to PySCF's 2/3-center integrals;
    ///   shape presented in col-major.
    /// - Row-major: the same **shape** to PySCF's integrals; shape presented in
    ///   row-major.
    #[builder(default = false)]
    pub row_major: bool,
}

/// Integrate arguments for function [`CInt::integrate_cross_with_args`] and
/// [`CInt::integrate_cross_with_args_spinor`].
///
/// Builder of this struct can be retrived by
/// [`CInt::integrate_cross_args_builder`]
/// and [`CInt::integrate_cross_args_builder_spinor`].
#[derive(Builder, Debug)]
#[builder(pattern = "owned", build_fn(error = "CIntError"))]
pub struct IntorCrossArgs<'l, F> {
    /// Integrator name, such as "int1e_ovlp", "int2e", etc.
    #[builder(setter(into))]
    pub intor: &'l str,

    /// Molecules for which the integral will be evaluated.
    ///
    /// The length of molecules must be accordance with the `shls_slice` and
    /// number of centers of integrator.
    pub mols: &'l [&'l CInt],

    /// Shell slices for evaluating sub-tensor of the total integral.
    #[builder(default = "&[]")]
    pub shls_slice: &'l [[usize; 2]],

    /// Symmetry of the integral. This will affect shape of output tensor.
    #[builder(default, setter(into))]
    pub aosym: CIntSymm,

    /// Output buffer for the integral evaluation.
    ///
    /// If this is `None`, the integral engine will generate a buffer for
    /// output.
    #[builder(default, setter(strip_option))]
    pub out: Option<&'l mut [F]>,
}

/* #endregion */

/* #region get_integrator */

/// Obtain integrator by name.
///
/// # Panics
///
/// - integrator not found
pub fn get_integrator(intor: &str) -> Box<dyn Integrator> {
    get_integrator_f(intor).unwrap()
}

pub fn get_integrator_f(intor: &str) -> Result<Box<dyn Integrator>, CIntError> {
    // explicitly check cint and ecp differently
    let intor = if let Some(intor) = get_cint_integrator(&intor.to_lowercase()) {
        Ok(intor)
    } else if let Some(intor) = get_ecp_integrator(&intor.to_lowercase()) {
        Ok(intor)
    } else {
        Err(CIntError::IntegratorNotFound(intor.to_string()))
    }?;

    // should have checked anything, gracefully returns
    Ok(intor)
}

/* #endregion */

/* #region CIntError impl */

impl Display for CIntError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        Debug::fmt(self, f)
    }
}

impl Error for CIntError {}

impl From<UninitializedFieldError> for CIntError {
    fn from(err: UninitializedFieldError) -> Self {
        CIntError::UninitializedFieldError(err)
    }
}

/* #endregion */

/* #region CIntType impl */

impl From<&str> for CIntType {
    #[inline]
    fn from(cint_type: &str) -> Self {
        match cint_type.to_lowercase().as_str() {
            "sph" | "spheric" | "spherical" => CIntType::Spheric,
            "cart" | "cartesian" => CIntType::Cartesian,
            "spinor" => CIntType::Spinor,
            _ => panic!("Unknown integral type: {cint_type}"),
        }
    }
}

/* #endregion */

/* #region CIntSymm impl */

impl From<&str> for CIntSymm {
    #[inline]
    fn from(aosym: &str) -> Self {
        match aosym {
            "S1" | "s1" => CIntSymm::S1,
            "S2ij" | "s2ij" | "S2uv" | "s2uv" => CIntSymm::S2ij,
            "S2kl" | "s2kl" => CIntSymm::S2kl,
            "S4" | "s4" => CIntSymm::S4,
            "S8" | "s8" => CIntSymm::S8,
            _ => panic!("Unknown symmetry: {aosym}"),
        }
    }
}

impl From<Option<&str>> for CIntSymm {
    #[inline]
    fn from(aosym: Option<&str>) -> Self {
        aosym.unwrap_or("s1").into()
    }
}

/* #endregion */

/* #region CIntOptimizer impl */

unsafe impl Send for CIntOptimizer {}
unsafe impl Sync for CIntOptimizer {}

impl Drop for CIntOptimizer {
    fn drop(&mut self) {
        match self {
            Self::Int(opt) => unsafe { cint_ffi::CINTdel_optimizer(opt) },
            Self::Ecp(opt) => unsafe { cecp_ffi::ECPdel_optimizer(opt) },
        }
    }
}

impl CIntOptimizer {
    #[inline]
    pub fn as_ptr(&self) -> *const c_void {
        match self {
            Self::Int(opt) => *opt as *const c_void,
            Self::Ecp(opt) => *opt as *const c_void,
        }
    }
}

/* #endregion */

/* #region CIntOutput impl */

impl<F> From<CIntOutput<F>> for (Vec<F>, Vec<usize>) {
    fn from(output: CIntOutput<F>) -> Self {
        if output.out.is_none() {
            panic!(
                "You have called integral with output buffer, and should not call `into` to get the integral buffer, but use your mutable reference instead."
            );
        }
        (output.out.unwrap_or_default(), output.shape)
    }
}

/* #endregion */

/// Implementation of integral at higher API level (for basic user usage).
impl CInt {
    /// Main electronic integral driver (for non-spinor type) in column-major.
    ///
    /// # PySCF equivalent
    ///
    /// `mol.intor(intor, aosym, shls_slice)`
    ///
    /// # Arguments
    ///
    /// - `intor`: name of the integral to be evaluated, such as `"int1e_ovlp"`,
    ///   `"int2e"`, `"ECPscalar_iprinvip"`, `"int2e_giao_sa10sp1spsp2"`, etc.
    /// - `aosym`: symmetry of the integral, you can specify `"s1"`, `"s2ij"`,
    ///   `"s2kl"`, `"s4"`, `"s8"`.
    /// - `shls_slice`: shell slices for evaluating sub-tensor of the total
    ///   integral. Please note that this is either be `None`, or a slice of
    ///   type `&[[usize; 2]]` or something similar (vectors, arrays) **with
    ///   length the same to the number of components in integral**.
    ///
    /// <div class="warning">
    ///
    /// **Does not check whether AO symmetry is correct to the correspoding
    /// integrator**
    ///
    /// For example,
    /// - `int2e_ip2` $(\mu \nu | \kappa \nabla \lambda)$ is meaningful for
    ///   `"s1"` and `"s2ij"` symmetry, not meaningful for other kinds of
    ///   symmetry.
    /// - `int2e_ip1` $(\nabla \mu \nu | \kappa \lambda)$ is meaningful for
    ///   `"s1"` and `"s2kl"` symmetry, not meaningful for other kinds of
    ///   symmetry.
    /// - `int2e` is meaningful for all symmetries.
    /// - `int1e_ovlp` is meaningful for `"s1"` and `"s2ij"` symmetry.
    ///
    /// However, this wrapper will generally not panic, nor raise error, if you
    /// use an incorrect symmetry for the integral. An exception is that if
    /// you specify `"s2kl"`, `"s4"`, `"s8"` on non-4-center integrals, it will
    /// raise error.
    ///
    /// </div>
    ///
    /// # Outputs
    ///
    /// Returns a [`CIntOutput<f64>`] for the integral evaluation, which
    /// contains
    /// - `out` (`Vec<f64>`): the output buffer for the integral evaluation.
    /// - `shape` (`Vec<usize>`): the shape of the output buffer.
    ///
    /// You can use `.into()` to convert it to a tuple of `(Vec<f64>,
    /// Vec<usize>)`:
    ///
    /// ```rust
    /// # use libcint::prelude::*;
    /// # let cint_data = init_h2o_def2_tzvp();
    /// let (out, shape) = cint_data.integrate("int1e_ovlp", "s2ij", None).into();
    /// ```
    ///
    /// # Examples
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
    /// // aosym default: "s1"
    /// // shls_slice default: None (evaluating all shells)
    /// let (out, shape) = cint_data.integrate("int1e_ovlp", None, None).into();
    /// assert_eq!(shape, vec![43, 43]);
    ///
    /// // utilize "s2ij" symmetry
    /// let (out, shape) = cint_data.integrate("int3c2e_ip2", "s2ij", None).into();
    /// assert_eq!(shape, vec![946, 43, 3]);
    ///
    /// // calculate sub-tensor of the integral
    /// let nbas = cint_data.nbas();
    /// let shl_slices = [[0, nbas], [0, nbas], [3, 12]]; // 3-centers for int3c2e_ip2
    /// let (out, shape) = cint_data.integrate("int3c2e_ip2", "s2ij", shl_slices).into();
    /// assert_eq!(shape, vec![946, 29, 3]);
    ///
    /// let shls_slice = [[0, nbas], [0, nbas], [3, 12], [4, 15]]; // 4-centers for int2e_ip2
    /// let (out, shape) = cint_data.integrate("int2e_ip2", "s2ij", shls_slice).into();
    /// assert_eq!(shape, vec![946, 29, 33, 3]);
    /// ```
    ///
    /// ECP (effective core potential) integrals are also supported.
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let cint_data = init_sb2me4_cc_pvtz();
    ///
    /// let (out, shape) = cint_data.integrate("ECPscalar_iprinvip", "s1", None).into();
    /// assert_eq!(shape, vec![366, 366, 9]);
    /// ```
    ///
    /// # Column-major convention
    ///
    /// <div class="warning">
    ///
    /// [`CInt`] always returns column-major shape with its output buffer.
    ///
    /// **Please note that [`CInt`] is somehow different to PySCF's
    /// convention.** The user should be very clear about the convention
    /// before using this function, especially when you are evaluating
    /// integrals with derivatives or using `shls_slice` option.
    ///
    /// </div>
    ///
    /// The following table summarizes difference between PySCF and [`CInt`]
    /// conventions.
    /// - "Strides" are indicated by numbers. Recall that a tensor is defined by
    ///   its data buffer in RAM and layout; the layouot contains shape, strides
    ///   and offset.
    ///
    ///   Those numbers are not the actual stride size, but actually indicates
    ///   the relative size of that axis in  tensor. For row-major tensors of
    ///   3-dimensional tensor, it is $[2, 1, 0]$; for col-major tensors, it is
    ///   $[0, 1, 2]$.
    /// - $t$ indicates the number of components in the integral.
    /// - $\mu, \nu, \kappa, \lambda$ indicate the indices of the atomic
    ///   orbitals (basis).
    /// - $\mathrm{tp}(\mu \nu)$ means triangular-packed shape of the $(\mu,
    ///   \nu)$ indices pair. This is always packed by upper-triangular indices
    ///   for col-major, or equilvalently lower-triangular indices for
    ///   row-major.
    /// - "same data layout" means the underlying data is the same, regardless
    ///   of how shape and strides are defined. If same data layout, it means
    ///   [`CInt`] data is the transposed PySCF's data, without any additional
    ///   copy or memory allocation.
    ///
    /// | centers | comp | symm   | PySCF shape | PySCF strides | [`CInt`] shape | [`CInt`] strides | same data layout |
    /// |---------|------|--------|-------------|---------------|----------------|------------------|------------------|
    /// | 2       | ✘    | `s1`   | $[\mu, \nu]$                                   | $[0, 1]$          | $[\mu, \nu]$                                  | $[0, 1]$          | ✔ |
    /// | 2       | ✔    | `s1`   | $[t, \mu, \nu]$                                | $[2, 0, 1]$       | $[\mu, \nu, t]$                               | $[0, 1, 2]$       | ✔ |
    /// | 2       | ✘    | `s2ij` | Not supported                                  | -                 | $[\mathrm{tp}(\mu \nu)]$                      | $\[0\]$           | - |
    /// | 2       | ✔    | `s2ij` | Not supported                                  | -                 | $[\mathrm{tp}(\mu \nu), t]$                   | $[0, 1]$          | - |
    /// | 3       | ✘    | `s1`   | $[\mu, \nu, \kappa]$                           | $[0, 1, 2]$       | $[\mu, \nu, \kappa]$                          | $[0, 1, 2]$       | ✔ |
    /// | 3       | ✔    | `s1`   | $[t, \mu, \nu, \kappa]$                        | $[3, 0, 1, 2]$    | $[\mu, \nu, \kappa, t]$                       | $[0, 1, 2, 3]$    | ✔ |
    /// | 3       | ✘    | `s2ij` | $[\mathrm{tp}(\mu \nu), \kappa]$               | $[0, 1]$          | $[\mathrm{tp}(\mu \nu), \kappa]$              | $[0, 1]$          | ✔ |
    /// | 3       | ✔    | `s2ij` | $[t, \mathrm{tp}(\mu \nu), \kappa]$            | $[2, 0, 1]$       | $[\mathrm{tp}(\mu \nu), \kappa, t]$           | $[0, 1, 2]$       | ✔ |
    /// | 4       | ✘    | `s1`   | $[\mu, \nu, \kappa, \lambda]$                  | $[3, 2, 1, 0]$    | $[\mu, \nu, \kappa, \lambda]$                 | $[0, 1, 2, 3]$    | ✘ |
    /// | 4       | ✔    | `s1`   | $[t, \mu, \nu, \kappa, \lambda]$               | $[4, 3, 2, 1, 0]$ | $[\mu, \nu, \kappa, \lambda, t]$              | $[0, 1, 2, 3, 4]$ | ✘ |
    /// | 4       | ✘    | `s2ij` | $[\mathrm{tp}(\mu \nu), \kappa, \lambda]$      | $[2, 1, 0]$       | $[\mathrm{tp}(\mu \nu), \kappa, \lambda]$     | $[0, 1, 2]$       | ✘ |
    /// | 4       | ✔    | `s2ij` | $[t, \mathrm{tp}(\mu \nu), \kappa, \lambda]$   | $[3, 2, 1, 0]$    | $[\mathrm{tp}(\mu \nu), \kappa, \lambda, t]$  | $[0, 1, 2, 3]$    | ✘ |
    /// | 4       | ✘    | `s2kl` | $[\mu, \nu, \mathrm{tp}(\kappa \lambda)]$      | $[2, 1, 0]$       | $[\mu, \nu, \mathrm{tp}(\kappa \lambda)]$     | $[0, 1, 2]$       | ✘ |
    /// | 4       | ✔    | `s2kl` | $[t, \mu, \nu, \mathrm{tp}(\kappa \lambda)]$   | $[3, 2, 1, 0]$    | $[\mu, \nu, \mathrm{tp}(\kappa \lambda), t]$  | $[0, 1, 2, 3]$    | ✘ |
    /// | 4       | ✘    | `s4`   | $[\mathrm{tp}(\mu \nu), \mathrm{tp}(\kappa \lambda)]$ | $[1, 0]$ | $[\mathrm{tp}(\mu \nu), \mathrm{tp}(\kappa \lambda)]$ | $[0, 1]$ | ✘ |
    /// | 4       | ✘    | `s8`   | $[\mathrm{tp}(\mathrm{tp}(\mu \nu) \mathrm{tp}(\kappa \lambda))]$ | $\[0\]$ | $[\mathrm{tp}(\mathrm{tp}(\mu \nu) \mathrm{tp}(\kappa \lambda))]$ | $\[0\]$ | ✔ |
    ///
    /// # Spheric or Cartesian
    ///
    /// The `CInt` supports both spherical and cartesian integrals. You can
    /// either set [`CInt::cint_type`] or use `_sph` or `_cart` suffix in
    /// integrator name argument `intor`:
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let mut cint_data = init_h2o_def2_tzvp();
    ///
    /// // spherical integral
    /// let (out, shape) = cint_data.integrate("int1e_ovlp", None, None).into();
    /// assert_eq!(shape, vec![43, 43]);
    ///
    /// // cartesian integral using with-clause
    /// let (out, shape) = cint_data.with_cint_type("cart", |cint_data| {
    ///     cint_data.integrate("int1e_ovlp", None, None).into()
    /// });
    /// assert_eq!(shape, vec![48, 48]);
    ///
    /// // cartesian integral using suffix
    /// let (out, shape) = cint_data.integrate("int1e_ovlp_cart", None, None).into();
    /// assert_eq!(shape, vec![48, 48]);
    /// ```
    ///
    /// <div class="warning">
    ///
    /// **This function cannot handle spinor integral**
    ///
    /// Due to rust's strict type system, it is very difficult to handle both
    /// spinor and spheric/cartesian integrals in the same function. So spinor
    /// integral should be called by
    /// [`integrate_spinor`](Self::integrate_spinor), which output is
    /// [`Complex<f64>`].
    ///
    /// </div>
    ///
    /// ```should_panic
    /// # use libcint::prelude::*;
    /// # let cint_data = init_h2o_def2_tzvp();
    /// let (out, shape) = cint_data.integrate("int1e_ovlp_spinor", None, None).into();
    /// // panics: Expected float type size 16 bytes, but got 8 bytes.
    /// ```
    ///
    /// # See also
    ///
    /// - [`integrate_f`](Self::integrate_f) for fallible counterpart.
    /// - [`integrate_spinor`](Self::integrate_spinor) for spinor type
    ///   integrals.
    /// - [`integrate_cross`](Self::integrate_cross) for cross integrals between
    ///   multiple molecules.
    /// - [`integrate_with_args`](Self::integrate_with_args) and
    ///   [`integrate_cross_with_args`](Self::integrate_cross_with_args) for
    ///   more advanced usage (full arguments that this crate supports).
    pub fn integrate(&self, intor: &str, aosym: impl Into<CIntSymm>, shls_slice: impl Into<ShlsSlice>) -> CIntOutput<f64> {
        self.integrate_f(intor, aosym, shls_slice.into()).unwrap()
    }

    /// Main electronic integral driver (for spinor type) in column-major.
    ///
    /// We refer most documentation to [`integrate`](Self::integrate).
    ///
    /// # PySCF equivalent
    ///
    /// `mol.intor(f"{intor}_spinor", aosym, shls_slice)`
    ///
    /// # Examples
    ///
    /// The following example uses a pre-defined water molecule with def2-TZVP
    /// basis set.
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let mut cint_data = init_h2o_def2_tzvp();
    ///
    /// // spinor integral using with-clause
    /// let (out, shape) = cint_data.with_cint_type("spinor", |cint_data| {
    ///    cint_data.integrate_spinor("int1e_ovlp", None, None).into()
    /// });
    /// assert_eq!(shape, vec![86, 86]);
    ///
    /// // spinor integral using suffix
    /// let (out, shape) = cint_data.integrate_spinor("int1e_ovlp_spinor", None, None).into();
    /// assert_eq!(shape, vec![86, 86]);
    /// ```
    ///
    /// <div class="warning">
    ///
    /// **This function cannot handle spheric or cartesian integral**
    ///
    /// Due to rust's strict type system, it is very difficult to handle both
    /// spinor and spheric/cartesian integrals in the same function. So spheric
    /// and cartesian integral should be called by
    /// [`integrate`](Self::integrate), which output is `f64`.
    ///
    /// </div>
    ///
    /// ```should_panic
    /// # use libcint::prelude::*;
    /// # let cint_data = init_h2o_def2_tzvp();
    /// let (out, shape) = cint_data.integrate_spinor("int1e_ovlp_sph", None, None).into();
    /// // panics: Expected float type size 8 bytes, but got 16 bytes.
    /// ```
    ///
    /// # See also
    ///
    /// - [`integrate`](Self::integrate) for non-spinor type integrals.
    /// - [`integrate_spinor_f`](Self::integrate_spinor_f) for fallible
    ///   counterpart.
    /// - [`integrate_cross_spinor`](Self::integrate_cross_spinor) for cross
    ///   integrals between multiple molecules.
    /// - [`integrate_with_args_spinor`](Self::integrate_with_args_spinor) and
    ///   [`integrate_cross_with_args_spinor`](Self::integrate_cross_with_args_spinor)
    ///   for more advanced usage (full arguments that this crate supports).
    pub fn integrate_spinor(&self, intor: &str, aosym: impl Into<CIntSymm>, shls_slice: impl Into<ShlsSlice>) -> CIntOutput<Complex<f64>> {
        self.integrate_spinor_f(intor, aosym, shls_slice.into()).unwrap()
    }

    /// Main electronic integral driver (for non-spinor type) in column-major.
    ///
    /// This function is fallible.
    ///
    /// # See also
    ///
    /// - [`integrate`](Self::integrate) for non-fallible counterpart.
    pub fn integrate_f(
        &self,
        intor: &str,
        aosym: impl Into<CIntSymm>,
        shls_slice: impl Into<ShlsSlice>,
    ) -> Result<CIntOutput<f64>, CIntError> {
        let shls_slice = shls_slice.into();
        let integrate_args = self.integrate_args_builder().intor(intor).aosym(aosym).shls_slice(shls_slice.as_ref()).build()?;
        self.integrate_with_args_inner(integrate_args)
    }

    /// Main electronic integral driver (for spinor type) in column-major.
    ///
    /// This function is fallible.
    ///
    /// # See also
    ///
    /// - [`integrate_spinor`](Self::integrate_spinor) for non-fallible
    ///   counterpart.
    pub fn integrate_spinor_f(
        &self,
        intor: &str,
        aosym: impl Into<CIntSymm>,
        shls_slice: impl Into<ShlsSlice>,
    ) -> Result<CIntOutput<Complex<f64>>, CIntError> {
        let shls_slice = shls_slice.into();
        let integrate_args = self.integrate_args_builder_spinor().intor(intor).aosym(aosym).shls_slice(shls_slice.as_ref()).build()?;
        self.integrate_with_args_inner(integrate_args)
    }

    pub fn integrate_cross<'l>(
        intor: &str,
        mols: impl AsRef<[&'l CInt]>,
        aosym: impl Into<CIntSymm>,
        shls_slice: impl Into<ShlsSlice>,
    ) -> CIntOutput<f64> {
        CInt::integrate_cross_f(intor, mols, aosym, shls_slice).unwrap()
    }

    pub fn integrate_cross_f<'l>(
        intor: &str,
        mols: impl AsRef<[&'l CInt]>,
        aosym: impl Into<CIntSymm>,
        shls_slice: impl Into<ShlsSlice>,
    ) -> Result<CIntOutput<f64>, CIntError> {
        let shls_slice = shls_slice.into();
        let args = IntorCrossArgsBuilder::default()
            .intor(intor)
            .mols(mols.as_ref())
            .shls_slice(shls_slice.as_ref())
            .aosym(aosym.into())
            .build()?;
        CInt::integrate_cross_with_args_inner(args)
    }

    pub fn integrate_cross_spinor<'l>(
        intor: &str,
        mols: impl AsRef<[&'l CInt]>,
        aosym: impl Into<CIntSymm>,
        shls_slice: impl Into<ShlsSlice>,
    ) -> CIntOutput<Complex<f64>> {
        CInt::integrate_cross_spinor_f(intor, mols, aosym, shls_slice).unwrap()
    }

    pub fn integrate_cross_spinor_f<'l>(
        intor: &str,
        mols: impl AsRef<[&'l CInt]>,
        aosym: impl Into<CIntSymm>,
        shls_slice: impl Into<ShlsSlice>,
    ) -> Result<CIntOutput<Complex<f64>>, CIntError> {
        let shls_slice = shls_slice.into();
        let args = IntorCrossArgsBuilder::default()
            .intor(intor)
            .mols(mols.as_ref())
            .shls_slice(shls_slice.as_ref())
            .aosym(aosym.into())
            .build()?;
        CInt::integrate_cross_with_args_inner(args)
    }
}

impl CInt {
    /// Main electronic integral driver (for non-spinor type) in row-major.
    pub fn integrate_row_major(&self, intor: &str, aosym: impl Into<CIntSymm>, shls_slice: impl Into<ShlsSlice>) -> CIntOutput<f64> {
        self.integrate_row_major_f(intor, aosym, shls_slice.into()).unwrap()
    }

    /// Main electronic integral driver (for non-spinor type) in row-major.
    ///
    /// This function is fallible.
    pub fn integrate_row_major_f(
        &self,
        intor: &str,
        aosym: impl Into<CIntSymm>,
        shls_slice: impl Into<ShlsSlice>,
    ) -> Result<CIntOutput<f64>, CIntError> {
        let shls_slice = shls_slice.into();
        let integrate_args =
            self.integrate_args_builder().intor(intor).aosym(aosym).shls_slice(shls_slice.as_ref()).row_major(true).build()?;
        self.integrate_with_args_inner(integrate_args)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cint_data() {
        let cint_data = init_h2o_def2_tzvp();
        let opt = cint_data.optimizer("int2e");
        println!("int2e opt: {opt:?}");
        if let CIntOptimizer::Int(opt) = opt {
            println!("int2e opt inner {:?}", unsafe { opt.read() });
        }
    }
}
