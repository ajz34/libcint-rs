//! Getter/setter and with-clause temporary changes to [`CInt`] instance.

use crate::prelude::*;
use core::ops::Add;

/// Concating two molecules for integral.
impl Add for &CInt {
    type Output = CInt;

    fn add(self, mol2: &CInt) -> CInt {
        let mol1 = self.decopule_ecpbas();
        let mol2 = mol2.decopule_ecpbas();

        if mol1.cint_type != mol2.cint_type {
            panic!("Cannot concatenate CInt with different cint_type");
        }

        const PTR_COORD: usize = cint_ffi::PTR_COORD as usize;
        const PTR_ZETA: usize = cint_ffi::PTR_ZETA as usize;
        const ATOM_OF: usize = cint_ffi::ATOM_OF as usize;
        const PTR_EXP: usize = cint_ffi::PTR_EXP as usize;
        const PTR_COEFF: usize = cint_ffi::PTR_COEFF as usize;

        let offset_env = mol1.env.len() as c_int;
        let offset_atm = mol1.atm.len() as c_int;

        let atm2 = mol2
            .atm
            .iter()
            .map(|v| {
                let mut v = *v;
                v[PTR_COORD] += offset_env;
                v[PTR_ZETA] += offset_env;
                v
            })
            .collect_vec();
        let bas2 = mol2
            .bas
            .iter()
            .map(|v| {
                let mut v = *v;
                v[ATOM_OF] += offset_atm;
                v[PTR_EXP] += offset_env;
                v[PTR_COEFF] += offset_env;
                v
            })
            .collect_vec();
        let ecpbas2 = mol2
            .ecpbas
            .iter()
            .map(|v| {
                let mut v = *v;
                v[ATOM_OF] += offset_atm;
                v[PTR_EXP] += offset_env;
                v[PTR_COEFF] += offset_env;
                v
            })
            .collect_vec();

        let mut result = mol1;
        result.atm.extend(atm2);
        result.bas.extend(bas2);
        result.ecpbas.extend(ecpbas2);
        result.env.extend_from_slice(&mol2.env);

        result
    }
}

/* #region fakemol */

/// Create a fake molecule for gaussian distributed charges.
///
/// # See also
///
/// [`CInt::fakemol_for_charges`]
pub trait FakeMolForChargeArg {
    fn fakemol_for_charges(coords: &[[f64; 3]], arg: Self) -> CInt;
}

/// Create a fake molecule for gaussian distributed charges.
///
/// # See also
///
/// [`CInt::fakemol_for_charges`]
pub fn fakemol_for_charges<T>(coords: &[[f64; 3]], arg: T) -> CInt
where
    T: FakeMolForChargeArg,
{
    T::fakemol_for_charges(coords, arg)
}

impl FakeMolForChargeArg for f64 {
    fn fakemol_for_charges(coords: &[[f64; 3]], exponent: Self) -> CInt {
        const ATM_SLOTS: usize = cint_ffi::ATM_SLOTS as usize;
        const BAS_SLOTS: usize = cint_ffi::BAS_SLOTS as usize;
        const PTR_ENV_START: usize = cint_ffi::PTR_ENV_START as usize;
        const PTR_COORD: usize = cint_ffi::PTR_COORD as usize;
        const ATOM_OF: usize = cint_ffi::ATOM_OF as usize;
        const NPRIM_OF: usize = cint_ffi::NPRIM_OF as usize;
        const NCTR_OF: usize = cint_ffi::NCTR_OF as usize;
        const PTR_EXP: usize = cint_ffi::PTR_EXP as usize;
        const PTR_COEFF: usize = cint_ffi::PTR_COEFF as usize;
        const PI: f64 = std::f64::consts::PI;

        let nbas = coords.len();
        let mut fake_atm = vec![];
        let mut fake_bas = vec![];
        let mut fake_env = vec![0.0; PTR_ENV_START];
        let mut ptr = PTR_ENV_START;

        // atm
        coords.iter().for_each(|&coord| {
            let mut atm = [0; ATM_SLOTS];
            atm[PTR_COORD] = ptr as c_int;
            fake_atm.push(atm);
            fake_env.extend_from_slice(&coord);
            ptr += 3;
        });

        // bas
        (0..nbas).for_each(|i| {
            let mut bas = [0; BAS_SLOTS];
            bas[ATOM_OF] = i as c_int;
            bas[NPRIM_OF] = 1;
            bas[NCTR_OF] = 1;
            bas[PTR_EXP] = ptr as c_int;
            bas[PTR_COEFF] = (ptr + 1) as c_int;
            fake_bas.push(bas);
        });

        // env
        let coef = 1.0 / (2.0 * PI.sqrt() * gaussian_int(2.0, exponent));
        fake_env.extend_from_slice(&[exponent, coef]);

        CInt { atm: fake_atm, bas: fake_bas, ecpbas: vec![], env: fake_env, cint_type: CIntType::default() }
    }
}

impl FakeMolForChargeArg for &[f64] {
    fn fakemol_for_charges(coords: &[[f64; 3]], exponents: Self) -> CInt {
        const ATM_SLOTS: usize = cint_ffi::ATM_SLOTS as usize;
        const BAS_SLOTS: usize = cint_ffi::BAS_SLOTS as usize;
        const PTR_ENV_START: usize = cint_ffi::PTR_ENV_START as usize;
        const PTR_COORD: usize = cint_ffi::PTR_COORD as usize;
        const ATOM_OF: usize = cint_ffi::ATOM_OF as usize;
        const NPRIM_OF: usize = cint_ffi::NPRIM_OF as usize;
        const NCTR_OF: usize = cint_ffi::NCTR_OF as usize;
        const PTR_EXP: usize = cint_ffi::PTR_EXP as usize;
        const PTR_COEFF: usize = cint_ffi::PTR_COEFF as usize;
        const PI: f64 = std::f64::consts::PI;

        let nbas = coords.len();
        if nbas != exponents.len() {
            panic!("Number of coordinates must match number of exponents");
        }

        let mut fake_atm = vec![];
        let mut fake_bas = vec![];
        let mut fake_env = vec![0.0; PTR_ENV_START];
        let mut ptr = PTR_ENV_START;

        // atm
        coords.iter().for_each(|&coord| {
            let mut atm = [0; ATM_SLOTS];
            atm[PTR_COORD] = ptr as c_int;
            fake_atm.push(atm);
            fake_env.extend_from_slice(&coord);
            ptr += 3;
        });

        // bas, env
        exponents.iter().enumerate().for_each(|(i, &exponent)| {
            let mut bas = [0; BAS_SLOTS];
            bas[ATOM_OF] = i as c_int;
            bas[NPRIM_OF] = 1;
            bas[NCTR_OF] = 1;
            bas[PTR_EXP] = ptr as c_int;
            bas[PTR_COEFF] = (ptr + 1) as c_int;
            let coef = 1.0 / (2.0 * PI.sqrt() * gaussian_int(2.0, exponent));
            fake_bas.push(bas);
            fake_env.extend_from_slice(&[exponent, coef]);
            ptr += 2;
        });

        CInt { atm: fake_atm, bas: fake_bas, ecpbas: vec![], env: fake_env, cint_type: CIntType::default() }
    }
}

impl FakeMolForChargeArg for Option<f64> {
    fn fakemol_for_charges(coords: &[[f64; 3]], exponent: Self) -> CInt {
        match exponent {
            Some(exponent) => fakemol_for_charges(coords, exponent),
            None => fakemol_for_charges(coords, 1.0e+16),
        }
    }
}

/// Create a fake molecule for contracted GTO-type charge.
///
/// # See also
///
/// [`CInt::fakemol_for_cgtf_charge`]
pub fn fakemol_for_cgtf_charge(coord: [f64; 3], exponents: &[f64], coeffs: &[f64]) -> CInt {
    const ATM_SLOTS: usize = cint_ffi::ATM_SLOTS as usize;
    const BAS_SLOTS: usize = cint_ffi::BAS_SLOTS as usize;
    const PTR_ENV_START: usize = cint_ffi::PTR_ENV_START as usize;
    const PTR_COORD: usize = cint_ffi::PTR_COORD as usize;
    const ATOM_OF: usize = cint_ffi::ATOM_OF as usize;
    const NPRIM_OF: usize = cint_ffi::NPRIM_OF as usize;
    const NCTR_OF: usize = cint_ffi::NCTR_OF as usize;
    const PTR_EXP: usize = cint_ffi::PTR_EXP as usize;
    const PTR_COEFF: usize = cint_ffi::PTR_COEFF as usize;
    const PI: f64 = core::f64::consts::PI;

    let nprim = exponents.len();
    if coeffs.len() != nprim {
        panic!("Number of exponents and coefficients must match for CGTF charge");
    }

    let mut fake_atm = vec![];
    let mut fake_bas = vec![];
    let mut fake_env = vec![0.0; PTR_ENV_START];
    let mut ptr = PTR_ENV_START;

    // atm
    let mut atm = [0; ATM_SLOTS];
    atm[PTR_COORD] = ptr as c_int;
    fake_atm.push(atm);
    fake_env.extend_from_slice(&coord);
    ptr += 3;

    // bas, env
    let mut bas = [0; BAS_SLOTS];
    bas[ATOM_OF] = 0; // single atom
    bas[NPRIM_OF] = nprim as c_int;
    bas[NCTR_OF] = 1; // single charge
    bas[PTR_EXP] = ptr as c_int;
    bas[PTR_COEFF] = (ptr + nprim) as c_int;

    let coeffs_converted = coeffs.iter().zip(exponents).map(|(&c, &e)| c / (2.0 * PI.sqrt() * gaussian_int(2.0, e))).collect_vec();

    fake_bas.push(bas);

    fake_env.extend_from_slice(exponents);
    fake_env.extend_from_slice(&coeffs_converted);

    CInt { atm: fake_atm, bas: fake_bas, ecpbas: vec![], env: fake_env, cint_type: CIntType::default() }
}

impl CInt {
    /// Create a fake molecule for gaussian distributed charges.
    ///
    /// You may also fake point charges by setting an extremely large exponent.
    ///
    /// # PySCF equivalent
    ///
    /// `gto.fakemol_for_charges(coords, arg)`
    ///
    /// # Argument overloads
    ///
    /// This function passes coordinates in Bohr (a.u.) unit, with `&[[f64; 3]]`
    /// type.
    ///
    /// For the second argument,
    /// - `None`: use a large exponent (1e+16) to mimic point-charge.
    /// - `f64`: use a single exponent for all coordinates.
    /// - `&[f64]`: use a list of exponents, one for each coordinate.
    /// - `(&[f64], &[f64])`: use a list of exponents and contraction
    ///   coefficients, one for each coordinate.
    ///
    /// # Examples
    ///
    /// All examples are initialized with by
    /// ```rust
    /// use libcint::prelude::*;
    /// let data_tzvp = init_h2o_def2_tzvp();
    /// let coords_chg = [[0., 1., 2.], [2., 0., 1.], [0., 2., 1.]];
    /// ```
    ///
    /// ## Mimic point-charge by large exponent
    ///
    /// This crate implementation:
    ///
    /// ```rust
    /// # use libcint::prelude::*;
    /// # let data_tzvp = init_h2o_def2_tzvp();
    /// # let coords_chg = [[0., 1., 2.], [2., 0., 1.], [0., 2., 1.]];
    /// let data_chg = CInt::fakemol_for_charges(&coords_chg, None);
    /// let (out, shape) = CInt::integrate_cross_row_major(
    ///     "int1e_ovlp", [&data_tzvp, &data_chg], None, None).into();
    /// assert_eq!(shape, [43, 3]);
    /// assert!((cint_fp(&out) - -0.12650556883004238).abs() < 1e-10);
    /// ```
    ///
    /// PySCF equivalent:
    ///
    /// ```python
    /// coords_chg = np.asarray([[0, 1, 2], [2, 0, 1], [0, 2, 1]])
    /// mol_chg = gto.fakemol_for_charges(coords_chg)
    /// out = gto.intor_cross("int1e_ovlp", mol_tzvp, mol_chg)
    /// lib.fp(out), out.shape
    /// ```
    ///
    /// ## Universal exponent for all coordinates
    ///
    /// ```rust
    /// # use libcint::prelude::*;
    /// # let data_tzvp = init_h2o_def2_tzvp();
    /// # let coords_chg = [[0., 1., 2.], [2., 0., 1.], [0., 2., 1.]];
    /// let exp_chg = 1.0;
    /// let data_chg = CInt::fakemol_for_charges(&coords_chg, exp_chg);
    /// let (out, shape) = CInt::integrate_cross_row_major(
    ///     "int1e_ovlp", [&data_tzvp, &data_chg], None, None).into();
    /// assert_eq!(shape, [43, 3]);
    /// assert!((cint_fp(&out) - 0.0707265752256318).abs() < 1e-10);
    /// ```
    ///
    /// PySCF equivalent:
    ///
    /// ```python
    /// exp_chg = 1.0
    /// mol_chg = gto.fakemol_for_charges(coords_chg, exp_chg)
    /// out = gto.intor_cross("int1e_ovlp", mol_tzvp, mol_chg)
    /// lib.fp(out), out.shape
    /// ```
    ///
    /// ## Exponents for each coordinate
    ///
    /// This crate implementation:
    ///
    /// ```rust
    /// # use libcint::prelude::*;
    /// # let data_tzvp = init_h2o_def2_tzvp();
    /// # let coords_chg = [[0., 1., 2.], [2., 0., 1.], [0., 2., 1.]];
    /// let exp_chg = [1.0, 2.5, 4.9];
    /// let data_chg = CInt::fakemol_for_charges(&coords_chg, exp_chg.as_slice());
    /// let (out, shape) = CInt::integrate_cross_row_major(
    ///     "int1e_ovlp", [&data_tzvp, &data_chg], None, None).into();
    /// assert_eq!(shape, [43, 3]);
    /// assert!((cint_fp(&out) - 0.0424629237780389).abs() < 1e-10);
    /// ```
    ///
    /// PySCF equivalent:
    ///
    /// ```python
    /// exp_chg = [1.0, 2.5, 4.9]
    /// mol_chg = gto.fakemol_for_charges(coords_chg, exp_chg)
    /// out = gto.intor_cross("int1e_ovlp", mol_tzvp, mol_chg)
    /// lib.fp(out), out.shape
    /// ```
    pub fn fakemol_for_charges<T>(coords: &[[f64; 3]], arg: T) -> CInt
    where
        T: FakeMolForChargeArg,
    {
        T::fakemol_for_charges(coords, arg)
    }

    /// Create a fake molecule for contracted GTO-type charge.
    ///
    /// # PySCF equivalent
    ///
    /// `gto.fakemol_for_cgtf_charge(coord, exponents, coeffs)`
    ///
    /// # Examples
    ///
    /// This crate implementation:
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let data_tzvp = init_h2o_def2_tzvp();
    ///
    /// let coord = [0., 1., 2.];
    /// let exp_chg = [1.0, 2.5, 4.9];
    /// let coef_chg = [2.8, 3.3, 0.7];
    /// let data_chg = CInt::fakemol_for_cgtf_charge(coord, &exp_chg, &coef_chg);
    /// let (out, shape) = CInt::integrate_cross_row_major(
    ///     "int1e_ovlp", [&data_tzvp, &data_chg], None, None).into();
    /// assert_eq!(shape, [43, 1]);
    /// assert!((cint_fp(&out) - -0.054460537334674264).abs() < 1e-10);
    /// ```
    ///
    /// PySCF equivalent:
    ///
    /// ```python
    /// mol_tzvp = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="def2-TZVP").build()
    /// coord = np.asarray([[0., 1., 2.]])
    /// exp_chg = [1.0, 2.5, 4.9]
    /// coef_chg = [2.8, 3.3, 0.7]
    /// mol_chg = gto.fakemol_for_cgtf_charge(coord, exp_chg, coef_chg)
    /// out = gto.intor_cross("int1e_ovlp", mol_tzvp, mol_chg)
    /// lib.fp(out), out.shape
    /// ```
    pub fn fakemol_for_cgtf_charge(coord: [f64; 3], exponents: &[f64], coeffs: &[f64]) -> CInt {
        fakemol_for_cgtf_charge(coord, exponents, coeffs)
    }
}

#[test]
fn playground() {
    use crate::prelude::*;
    let data_tzvp = init_h2o_def2_tzvp();

    let coord = [0., 1., 2.];
    let exp_chg = [1.0, 2.5, 4.9];
    let coef_chg = [2.8, 3.3, 0.7];
    let data_chg = CInt::fakemol_for_cgtf_charge(coord, &exp_chg, &coef_chg);
    let (out, shape) = CInt::integrate_cross_row_major("int1e_ovlp", [&data_tzvp, &data_chg], None, None).into();
    assert_eq!(shape, [43, 1]);
    // println!("out: {out:?}");
    println!("env {:?}", data_chg.env);
    assert!((cint_fp(&out) - -0.054460537334674264).abs() < 1e-10);
}

/* #endregion */

/// Usual cases for mutating `CInt` instance temporarily (something similar to
/// with-clause in Python).
impl CInt {
    /* #region cint_type */

    /// Setter of default integral type.
    ///
    /// # See also
    ///
    /// - [`with_cint_type`](CInt::with_cint_type)
    pub fn set_cint_type(&mut self, cint_type: impl Into<CIntType>) -> &mut Self {
        self.cint_type = cint_type.into();
        self
    }

    /// With-clause of default integral type.
    ///
    /// This will change the type of integrals computed by default:
    ///
    /// ```rust
    /// use libcint::prelude::*;
    ///
    /// // this test mol is Spheric
    /// let mut cint_data = init_h2o_def2_tzvp();
    /// let (_, shape) = cint_data.integrate("int1e_ovlp", None, None).into();
    /// assert_eq!(cint_data.nao(), 43);
    /// assert_eq!(shape, [43, 43]);
    ///
    /// // change to Cartesian type
    /// cint_data.with_cint_type("cart", |data| {
    ///     let (_, shape) = data.integrate("int1e_ovlp", None, None).into();
    ///     assert_eq!(data.nao(), 48);
    ///     assert_eq!(shape, [48, 48]);
    /// });
    ///
    /// // change to spinor type
    /// cint_data.with_cint_type("spinor", |data| {
    ///     let (_, shape) = data.integrate_spinor("int1e_ovlp", None, None).into();
    ///     assert_eq!(data.nao(), 86);
    ///     assert_eq!(shape, [86, 86]);
    /// });
    /// ```
    ///
    /// # See also
    ///
    /// - [`set_cint_type`](CInt::set_cint_type)
    pub fn with_cint_type<R>(&mut self, cint_type: impl Into<CIntType>, func: impl FnOnce(&mut Self) -> R) -> R {
        let old_type = self.cint_type;
        self.set_cint_type(cint_type);
        let result = func(self);
        self.set_cint_type(old_type);
        result
    }

    /* #endregion */

    /* #region common_origin */

    /// Getter of common origin for integrals (in unit Bohr).
    ///
    /// # See also
    ///
    /// - [`with_common_origin`](CInt::with_common_origin)
    /// - [`set_common_origin`](CInt::set_common_origin)
    pub fn get_common_origin(&self) -> [f64; 3] {
        const PTR_COMMON_ORIG: usize = cint_ffi::PTR_COMMON_ORIG as usize;
        self.env[PTR_COMMON_ORIG..PTR_COMMON_ORIG + 3].try_into().unwrap()
    }

    /// Set the common origin for integrals (in unit Bohr).
    ///
    /// # See also
    ///
    /// - [`with_common_origin`](CInt::with_common_origin)
    /// - [`get_common_origin`](CInt::get_common_origin)
    pub fn set_common_origin(&mut self, origin: [f64; 3]) -> &mut Self {
        const PTR_COMMON_ORIG: usize = cint_ffi::PTR_COMMON_ORIG as usize;
        self.env[PTR_COMMON_ORIG..PTR_COMMON_ORIG + 3].copy_from_slice(&origin);
        self
    }

    /// Temporarily set the common origin for integrals (in unit Bohr).
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let mut cint_data = init_h2o_def2_tzvp();
    ///
    /// let (out, _) = cint_data.integrate_row_major("int1e_r", None, None).into();
    /// assert!((cint_fp(&out) - -0.7587292491644675).abs() < 1e-10);
    ///
    /// // set common origin to [0.0, 1.0, 2.0]
    /// cint_data.with_common_origin([0.0, 1.0, 2.0], |data| {
    ///     let (out, _) = data.integrate_row_major("int1e_r", None, None).into();
    ///     assert!((cint_fp(&out) - 71.88577867872883).abs() < 1e-10);
    /// });
    /// ```
    ///
    /// PySCF equivalent:
    ///
    /// ```python
    /// mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="def2-TZVP").build()
    /// assert abs(lib.fp(mol.intor("int1e_r")) - -0.7587292491644675) < 1e-10
    /// with mol.with_common_orig([0, 1, 2]):
    ///     assert abs(lib.fp(mol.intor("int1e_r")) - 71.88577867872883) < 1e-10
    /// ```
    ///
    /// # See also
    ///
    /// - [`set_common_origin`](CInt::set_common_origin)
    /// - [`get_common_origin`](CInt::get_common_origin)
    pub fn with_common_origin<R>(&mut self, origin: [f64; 3], func: impl FnOnce(&mut Self) -> R) -> R {
        let old_origin = self.get_common_origin();
        self.set_common_origin(origin);
        let result = func(self);
        self.set_common_origin(old_origin);
        result
    }

    /* #endregion */

    /* #region rinv_origin */

    /// Getter of origin in evaluating $1/r$ for integrals (in unit Bohr).
    ///
    /// # See also
    ///
    /// - [`with_rinv_origin`](CInt::with_rinv_origin)
    /// - [`set_rinv_origin`](CInt::set_rinv_origin)
    pub fn get_rinv_origin(&self) -> [f64; 3] {
        const PTR_RINV_ORIG: usize = cint_ffi::PTR_RINV_ORIG as usize;
        self.env[PTR_RINV_ORIG..PTR_RINV_ORIG + 3].try_into().unwrap()
    }

    /// Set the origin in evaluating $1/r$ for integrals (in unit Bohr).
    ///
    /// # See also
    ///
    /// - [`with_rinv_origin`](CInt::with_rinv_origin)
    /// - [`get_rinv_origin`](CInt::get_rinv_origin)
    pub fn set_rinv_origin(&mut self, origin: [f64; 3]) -> &mut Self {
        const PTR_RINV_ORIG: usize = cint_ffi::PTR_RINV_ORIG as usize;
        self.env[PTR_RINV_ORIG..PTR_RINV_ORIG + 3].copy_from_slice(&origin);
        self
    }

    /// Temporarily set the origin in evaluating $1/r$ for integrals (in unit
    /// Bohr).
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let mut cint_data = init_h2o_def2_tzvp();
    ///
    /// let (out, _) = cint_data.integrate_row_major("int1e_rinv", None, None).into();
    /// assert!((cint_fp(&out) - 51.806443495904794).abs() < 1e-10);
    ///
    /// // set rinv origin to [0.0, 1.0, 2.0]
    /// cint_data.with_rinv_origin([0.0, 1.0, 2.0], |data| {
    ///     let (out, _) = data.integrate_row_major("int1e_rinv", None, None).into();
    ///     assert!((cint_fp(&out) - 15.72929399764994).abs() < 1e-10);
    /// });
    /// ```
    ///
    /// PySCF equivalent:
    ///
    /// ```python
    /// mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="def2-TZVP").build()
    /// assert abs(lib.fp(mol.intor("int1e_rinv")) - 51.806443495904794) < 1e-10
    /// with mol.with_rinv_orig([0, 1, 2]):
    ///     assert abs(lib.fp(mol.intor("int1e_rinv")) - 15.72929399764994) < 1e-10
    /// ```
    ///
    /// # See also
    ///
    /// - [`set_rinv_origin`](CInt::set_rinv_origin)
    /// - [`get_rinv_origin`](CInt::get_rinv_origin)
    pub fn with_rinv_origin<R>(&mut self, origin: [f64; 3], func: impl FnOnce(&mut Self) -> R) -> R {
        let old_origin = self.get_rinv_origin();
        self.set_rinv_origin(origin);
        let result = func(self);
        self.set_rinv_origin(old_origin);
        result
    }

    /* #endregion */

    /* #region rinv_orig_atom */

    /// Getter of the origin atom in evaluating $1/r$ for ECP integrals.
    ///
    /// # See also
    ///
    /// - [`with_rinv_origin_atom`](CInt::with_rinv_origin_atom)
    /// - [`set_rinv_origin_atom`](CInt::set_rinv_origin_atom)
    pub fn get_rinv_origin_atom(&self) -> usize {
        const AS_RINV_ORIG_ATOM: usize = cecp_ffi::AS_RINV_ORIG_ATOM as usize;
        self.env[AS_RINV_ORIG_ATOM] as usize
    }

    /// Setter of the origin atom in evaluating $1/r$ for ECP integrals.
    ///
    /// # See also
    ///
    /// - [`with_rinv_origin_atom`](CInt::with_rinv_origin_atom)
    /// - [`get_rinv_origin_atom`](CInt::get_rinv_origin_atom)
    pub fn set_rinv_origin_atom(&mut self, atm_id: usize) -> &mut Self {
        const AS_RINV_ORIG_ATOM: usize = cecp_ffi::AS_RINV_ORIG_ATOM as usize;
        self.env[AS_RINV_ORIG_ATOM] = atm_id as f64;
        self
    }

    /// Temporarily set the origin atom in evaluating $1/r$ for ECP integrals.
    ///
    /// This option alone is not useful. See
    /// [`with_rinv_at_nucleus`](Self::with_rinv_at_nucleus) for more
    /// information.
    ///
    /// # See also
    ///
    /// - [`get_rinv_origin_atom`](CInt::get_rinv_origin_atom)
    /// - [`set_rinv_origin_atom`](CInt::set_rinv_origin_atom)
    pub fn with_rinv_origin_atom<R>(&mut self, atm_id: usize, func: impl FnOnce(&mut Self) -> R) -> R {
        let old_atm_id = self.get_rinv_origin_atom();
        self.set_rinv_origin_atom(atm_id);
        let result = func(self);
        self.set_rinv_origin_atom(old_atm_id);
        result
    }

    /* #endregion */

    /* #region range_coulomb */

    /// Getter of the range separation parameter $\omega$ for Coulomb integrals.
    ///
    /// # See also
    ///
    /// - [`with_range_coulomb`](CInt::with_range_coulomb)
    /// - [`set_range_coulomb`](CInt::set_range_coulomb)
    ///
    /// # Alias
    ///
    /// - [`omega`](CInt::omega)
    pub fn get_range_coulomb(&self) -> f64 {
        const PTR_RANGE_OMEGA: usize = cint_ffi::PTR_RANGE_OMEGA as usize;
        self.env[PTR_RANGE_OMEGA]
    }

    /// Setter the range separation parameter $\omega$ for Coulomb integrals.
    ///
    /// # See also
    ///
    /// - [`with_range_coulomb`](CInt::with_range_coulomb)
    /// - [`get_range_coulomb`](CInt::get_range_coulomb)
    ///
    /// # Alias
    ///
    /// - [`set_omega`](CInt::set_omega)
    pub fn set_range_coulomb(&mut self, range: f64) -> &mut Self {
        const PTR_RANGE_OMEGA: usize = cint_ffi::PTR_RANGE_OMEGA as usize;
        self.env[PTR_RANGE_OMEGA] = range;
        self
    }

    /// Temporarily set the range separation parameter $\omega$ for Coulomb
    /// integrals.
    ///
    /// Common operator of ERI (for example, `int2e`) is $1 / r$; however, in
    /// range separation scheme, the long-range operator is $\mathrm{erf}
    /// (\omega r) / r$, and the short-range operator is $\mathrm{erfc} (\omega
    /// r) / r$.
    ///
    /// For libcint convention,
    ///
    /// - Positive omega indicates long-range Coulomb operator.
    /// - Negative omega indicates short-range Coulomb operator.
    /// - Exactly zero omega indicates the original Coulomb operator $1/r$.
    ///
    /// # Example
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let mut cint_data = init_h2o_def2_tzvp();
    ///
    /// let (out, _) = cint_data.integrate_row_major("int2e", None, None).into();
    /// assert!((cint_fp(&out) - 70.00106603114841).abs() < 1e-10);
    ///
    /// // set range separation parameter to 0.5 (long-range)
    /// cint_data.with_range_coulomb(0.5, |data| {
    ///    let (out, _) = data.integrate_row_major("int2e", None, None).into();
    ///    assert!((cint_fp(&out) - 23.8282413132626).abs() < 1e-10);
    /// });
    ///
    /// // set range separation parameter to -0.5 (short-range)
    /// cint_data.with_range_coulomb(-0.5, |data| {
    ///     let (out, _) = data.integrate_row_major("int2e", None, None).into();
    ///     assert!((cint_fp(&out) - 46.17282471793578).abs() < 1e-10);
    /// });
    /// ```
    ///
    /// PySCF equivalent:
    /// ```python
    /// mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="def2-TZVP").build()
    /// assert abs(lib.fp(mol.intor("int2e")) - 70.00106603114841) < 1e-10
    /// with mol.with_range_coulomb(0.5):
    ///     assert abs(lib.fp(mol.intor("int2e")) - 23.8282413132626) < 1e-10
    /// with mol.with_range_coulomb(-0.5):
    ///     assert abs(lib.fp(mol.intor("int2e")) - 46.17282471793578) < 1e-10
    /// ```
    ///
    /// # See also
    ///
    /// - [`set_range_coulomb`](CInt::set_range_coulomb)
    /// - [`get_range_coulomb`](CInt::get_range_coulomb)
    ///
    /// # Alias
    ///
    /// - [`with_omega`](CInt::with_omega)
    pub fn with_range_coulomb<R>(&mut self, range: f64, func: impl FnOnce(&mut Self) -> R) -> R {
        let old_range = self.get_range_coulomb();
        self.set_range_coulomb(range);
        let result = func(self);
        self.set_range_coulomb(old_range);
        result
    }

    /// Getter of the range separation parameter $\omega$ for Coulomb integrals.
    ///
    /// # Alias
    ///
    /// - [`get_range_coulomb`](CInt::get_range_coulomb)
    pub fn omega(&self) -> f64 {
        self.get_range_coulomb()
    }

    /// Setter of the range separation parameter $\omega$ for Coulomb integrals.
    ///
    /// # Alias
    ///
    /// - [`set_range_coulomb`](CInt::set_range_coulomb)
    pub fn set_omega(&mut self, omega: f64) -> &mut Self {
        self.set_range_coulomb(omega)
    }

    /// Temporarily set the range separation parameter $\omega$ for Coulomb
    /// integrals.
    ///
    /// # Alias
    ///
    /// - [`with_range_coulomb`](CInt::with_range_coulomb)
    pub fn with_omega<R>(&mut self, omega: f64, func: impl FnOnce(&mut Self) -> R) -> R {
        self.with_range_coulomb(omega, func)
    }

    /// Temporarily set the range separation parameter $\omega$ for Coulomb
    /// integrals (long range).
    ///
    /// This function only accepts positive omega values.
    ///
    /// # Alias
    ///
    /// - [`with_range_coulomb`](CInt::with_range_coulomb)
    pub fn with_long_range_coulomb<R>(&mut self, omega: f64, func: impl FnOnce(&mut Self) -> R) -> R {
        if omega <= 0.0 {
            panic!("This function does not allow negative or zero omega values.");
        }
        self.with_range_coulomb(omega, func)
    }

    /// Temporarily set the range separation parameter $\omega$ for Coulomb
    /// integrals (short range).
    ///
    /// This function only accepts negative omega values, and is actually the
    /// alias to `self.with_range_coulomb(-omega, func)`.
    ///
    /// # Alias
    ///
    /// - [`with_range_coulomb`](CInt::with_range_coulomb)
    pub fn with_short_range_coulomb<R>(&mut self, omega: f64, func: impl FnOnce(&mut Self) -> R) -> R {
        if omega <= 0.0 {
            panic!("This function does not allow negative or zero omega values.");
        }
        self.with_range_coulomb(-omega, func)
    }

    /* #endregion */

    /* #region f12_zeta */

    #[cfg(feature = "with_f12")]
    pub fn get_f12_zeta(&self) -> f64 {
        const PTR_F12_ZETA: usize = cint_ffi::PTR_F12_ZETA as usize;
        self.env[PTR_F12_ZETA]
    }

    #[cfg(feature = "with_f12")]
    pub fn set_f12_zeta(&mut self, zeta: f64) -> &mut Self {
        const PTR_F12_ZETA: usize = cint_ffi::PTR_F12_ZETA as usize;
        self.env[PTR_F12_ZETA] = zeta;
        self
    }

    #[cfg(feature = "with_f12")]
    pub fn with_f12_zeta<R>(&mut self, zeta: f64, func: impl FnOnce(&mut Self) -> R) -> R {
        let old_zeta = self.get_f12_zeta();
        self.set_f12_zeta(zeta);
        let result = func(self);
        self.set_f12_zeta(old_zeta);
        result
    }

    /* #endregion */

    /* #region set_nuc_mod */

    /// Set the nuclear model for a given atom.
    pub fn set_nuc_mod(&mut self, atm_id: usize, zeta: f64) -> &mut Self {
        const PTR_ZETA: usize = cint_ffi::PTR_ZETA as usize;
        const NUC_MOD_OF: usize = cint_ffi::NUC_MOD_OF as usize;
        const POINT_NUC: c_int = cint_ffi::POINT_NUC as c_int;
        const GAUSSIAN_NUC: c_int = cint_ffi::GAUSSIAN_NUC as c_int;

        let ptr = self.atm[atm_id][PTR_ZETA] as usize;
        self.env[ptr] = zeta;
        self.atm[atm_id][NUC_MOD_OF] = if zeta == 0.0 { POINT_NUC } else { GAUSSIAN_NUC };
        self
    }

    /* #endregion */

    /* #region rinv_zeta */

    /// Getter of the $\zeta$ parameter for nucleus model.
    ///
    /// # See also
    ///
    /// - [`with_rinv_zeta`](CInt::with_rinv_zeta)
    /// - [`set_rinv_zeta`](CInt::set_rinv_zeta)
    pub fn get_rinv_zeta(&self) -> f64 {
        const PTR_RINV_ZETA: usize = cint_ffi::PTR_RINV_ZETA as usize;
        self.env[PTR_RINV_ZETA]
    }

    /// Setter of the $\zeta$ parameter for nucleus model.
    ///
    /// # See also
    ///
    /// - [`with_rinv_zeta`](CInt::with_rinv_zeta)
    /// - [`get_rinv_zeta`](CInt::get_rinv_zeta)
    pub fn set_rinv_zeta(&mut self, zeta: f64) -> &mut Self {
        const PTR_RINV_ZETA: usize = cint_ffi::PTR_RINV_ZETA as usize;
        self.env[PTR_RINV_ZETA] = zeta;
        self
    }

    /// Temporarily set the $\zeta$ parameter for nucleus model.
    ///
    /// # See also
    ///
    /// - [`get_rinv_zeta`](CInt::get_rinv_zeta)
    /// - [`set_rinv_zeta`](CInt::set_rinv_zeta)
    pub fn with_rinv_zeta<R>(&mut self, zeta: f64, func: impl FnOnce(&mut Self) -> R) -> R {
        let old_zeta = self.get_rinv_zeta();
        self.set_rinv_zeta(zeta);
        let result = func(self);
        self.set_rinv_zeta(old_zeta);
        result
    }

    /* #endregion */

    /* #region rinv_at_nucleus */

    /// Set the origin to atom coordinate in evaluating $1/r$ for integrals.
    ///
    /// # See also
    ///
    /// - [`with_rinv_at_nucleus`](CInt::with_rinv_at_nucleus)
    pub fn set_rinv_at_nucleus(&mut self, atm_id: usize) -> &mut Self {
        const PTR_ZETA: usize = cint_ffi::PTR_ZETA as usize;
        let zeta = self.env[self.atm[atm_id][PTR_ZETA] as usize];
        let rinv = self.atom_coord(atm_id);
        self.set_rinv_origin(rinv);
        self.set_rinv_origin_atom(atm_id);
        if zeta != 0.0 {
            self.set_rinv_zeta(zeta);
        }
        self
    }

    /// Temporarily set the origin to atom coordinate in evaluating $1/r$ for
    /// integrals.
    ///
    /// # Example (usual integral)
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let mut cint_data = init_h2o_def2_tzvp();
    ///
    /// let (out, _) = cint_data.integrate_row_major("int1e_rinv", None, None).into();
    /// assert!((cint_fp(&out) - 51.806443495904794).abs() < 1e-10);
    ///
    /// // set rinv origin to the second atom (first Hydrogen)
    /// cint_data.with_rinv_at_nucleus(1, |data| {
    ///     let (out, _) = data.integrate_row_major("int1e_rinv", None, None).into();
    ///     assert!((cint_fp(&out) - 20.940503856155193).abs() < 1e-10);
    /// });
    /// ```
    ///
    /// PySCF equivalent:
    ///
    /// ```python
    /// mol = gto.Mole(atom="O; H 1 0.94; H 1 0.94 2 104.5", basis="def2-TZVP").build()
    /// assert abs(lib.fp(mol.intor("int1e_rinv")) - 51.806443495904794) < 1e-10
    /// with mol.with_rinv_at_nucleus(1):
    ///     assert abs(lib.fp(mol.intor("int1e_rinv")) - 20.940503856155193) < 1e-10
    /// ```
    ///
    /// # Example (ECP integral)
    ///
    /// ```rust
    /// use libcint::prelude::*;
    /// let mut cint_data = init_sb2me4_cc_pvtz();
    ///
    /// let (out, _) = cint_data.integrate_row_major("ECPscalar_iprinvip", None, None).into();
    /// assert!((cint_fp(&out) - 324.13737563392084).abs() < 1e-10);
    ///
    /// // set rinv origin to the second atom (second Sb Antimony)
    /// cint_data.with_rinv_at_nucleus(1, |data| {
    ///     let (out, _) = data.integrate_row_major("ECPscalar_iprinvip", None, None).into();
    ///     assert!((cint_fp(&out) - 302.6772698217352).abs() < 1e-10);
    /// });
    /// ```
    ///
    /// PySCF equivalent:
    ///
    /// ```python
    /// // definition of mol, see function `init_sb2me4_cc_pvtz` in this crate
    /// assert abs(lib.fp(mol.intor("ECPscalar_iprinvip")) - 324.13737563392084) < 1e-10
    /// with mol.with_rinv_at_nucleus(1):
    ///     assert abs(lib.fp(mol.intor("ECPscalar_iprinvip")) - 302.6772698217352) < 1e-10
    /// ```
    ///
    /// # See also
    ///
    /// - [`set_rinv_at_nucleus`](CInt::set_rinv_at_nucleus)
    pub fn with_rinv_at_nucleus<R>(&mut self, atm_id: usize, func: impl FnOnce(&mut Self) -> R) -> R {
        let old_rinv = self.get_rinv_origin();
        let old_zeta = self.get_rinv_zeta();
        let old_rinv_atom = self.get_rinv_origin_atom();
        self.set_rinv_at_nucleus(atm_id);
        let result = func(self);
        self.set_rinv_origin(old_rinv);
        self.set_rinv_zeta(old_zeta);
        self.set_rinv_origin_atom(old_rinv_atom);
        result
    }

    /* #endregion */

    /* #region integral_screen */

    pub fn get_integral_screen(&self) -> f64 {
        const PTR_EXPCUTOFF: usize = cint_ffi::PTR_EXPCUTOFF as usize;
        self.env[PTR_EXPCUTOFF].exp()
    }

    pub fn set_integral_screen(&mut self, screen: f64) -> &mut Self {
        const PTR_EXPCUTOFF: usize = cint_ffi::PTR_EXPCUTOFF as usize;
        self.env[PTR_EXPCUTOFF] = screen.ln();
        self
    }

    pub fn with_integral_screen<R>(&mut self, screen: f64, func: impl FnOnce(&mut Self) -> R) -> R {
        let old_screen = self.get_integral_screen();
        self.set_integral_screen(screen);
        let result = func(self);
        self.set_integral_screen(old_screen);
        result
    }

    /* #endregion */

    /* #region set_geom */

    /// Getter of the geometry of the molecule (in unit Bohr).
    ///
    /// # See also
    ///
    /// - [`with_geom`](CInt::with_geom)
    /// - [`set_geom`](CInt::set_geom)
    ///
    /// # Alias
    ///
    /// - [`atom_coords`](CInt::atom_coords)
    pub fn get_geom(&self) -> Vec<[f64; 3]> {
        const PTR_COORD: usize = cint_ffi::PTR_COORD as usize;
        self.atm
            .iter()
            .map(|atm| {
                let ptr = atm[PTR_COORD] as usize;
                self.env[ptr..ptr + 3].try_into().unwrap()
            })
            .collect()
    }

    /// Setter of the geometry of the molecule (in unit Bohr).
    ///
    /// # See also
    ///
    /// - [`with_geom`](CInt::with_geom)
    /// - [`get_geom`](CInt::get_geom)
    pub fn set_geom(&mut self, coords: &[[f64; 3]]) -> &mut Self {
        const PTR_COORD: usize = cint_ffi::PTR_COORD as usize;

        if self.natm() != coords.len() {
            panic!("Number of coordinates ({}) does not match number of atoms ({})", coords.len(), self.natm());
        }

        coords.iter().enumerate().for_each(|(i, coord)| {
            let ptr = self.atm[i][PTR_COORD] as usize;
            self.env[ptr..ptr + 3].copy_from_slice(coord);
        });
        self
    }

    /// Temporarily set the geometry of the molecule (in unit Bohr).
    ///
    /// # See also
    ///
    /// - [`get_geom`](CInt::get_geom)
    /// - [`set_geom`](CInt::set_geom)
    pub fn with_geom<R>(&mut self, coords: &[[f64; 3]], func: impl FnOnce(&mut Self) -> R) -> R {
        let old_geom = self.get_geom();
        self.set_geom(coords);
        let result = func(self);
        self.set_geom(&old_geom);
        result
    }

    /* #endregion */

    /* #region set_grids */

    /// Getter of the grids for `int1e_grids` related integrals (in unit Bohr).
    ///
    /// # See also
    ///
    /// - [`with_grids`](CInt::with_grids)
    /// - [`set_grids`](CInt::set_grids)
    /// - [`try_clear_grids`](CInt::try_clear_grids)
    pub fn get_grids(&self) -> Vec<[f64; 3]> {
        const NGRIDS: usize = cint_ffi::NGRIDS as usize;
        const PTR_GRIDS: usize = cint_ffi::PTR_GRIDS as usize;

        let ngrids = self.env[NGRIDS] as usize;
        if ngrids == 0 {
            return vec![];
        }

        let ptr_start = self.env[PTR_GRIDS] as usize;
        let env_grids = &self.env[ptr_start..];
        (0..ngrids)
            .map(|i| {
                let ptr = ptr_start + i * 3;
                env_grids[ptr..ptr + 3].try_into().unwrap()
            })
            .collect()
    }

    /// Setter of the grids for `int1e_grids` related integrals (in unit Bohr).
    ///
    /// # See also
    ///
    /// - [`with_grids`](CInt::with_grids)
    /// - [`get_grids`](CInt::get_grids)
    /// - [`try_clear_grids`](CInt::try_clear_grids)
    pub fn set_grids(&mut self, grids: &[[f64; 3]]) -> &mut Self {
        const NGRIDS: usize = cint_ffi::NGRIDS as usize;
        const PTR_GRIDS: usize = cint_ffi::PTR_GRIDS as usize;

        let ngrids = grids.len();
        if ngrids == 0 {
            self.env[NGRIDS] = 0.0;
            self.env[PTR_GRIDS] = 0.0;
            return self;
        }

        self.env[NGRIDS] = ngrids as f64;
        let prev_env_size = self.env.len();
        self.env[PTR_GRIDS] = prev_env_size as f64;

        // cast grids to &[f64]
        let grids_slice = unsafe { core::slice::from_raw_parts(grids.as_ptr() as *const f64, ngrids * 3) };
        self.env.extend_from_slice(grids_slice);

        self
    }

    /// Try to clear the grids for `int1e_grids` related integrals.
    ///
    /// If grids are at the end of `env`, they will be cleared.
    /// Otherwise, the grids will become dangling data.
    ///
    /// # See also
    ///
    /// - [`with_grids`](CInt::with_grids)
    /// - [`set_grids`](CInt::set_grids)
    /// - [`get_grids`](CInt::get_grids)
    pub fn try_clear_grids(&mut self) -> &mut Self {
        const NGRIDS: usize = cint_ffi::NGRIDS as usize;
        const PTR_GRIDS: usize = cint_ffi::PTR_GRIDS as usize;

        let ngrids = self.env[NGRIDS] as usize;
        if ngrids != 0 {
            let ptr_start = self.env[PTR_GRIDS] as usize;
            let ptr_end = ptr_start + ngrids * 3;
            if ptr_end == self.env.len() {
                // If grids are at the end of env, we can safely clear them
                self.env.truncate(ptr_start);
                // otherwise, the grids become dangling data
            }
            self.env[NGRIDS] = 0.0;
            self.env[PTR_GRIDS] = 0.0;
        }

        self
    }

    /// Temporarily set the grids for `int1e_grids` related integrals (in unit
    /// Bohr).
    ///
    /// With this function and [`integrate`](CInt::integrate) or
    /// [`integrate_row_major`](CInt::integrate_row_major), you can perform
    /// `int1e_grids` related integrals with custom grids.
    ///
    /// # See also
    ///
    /// - [`get_grids`](CInt::get_grids)
    /// - [`set_grids`](CInt::set_grids)
    /// - [`try_clear_grids`](CInt::try_clear_grids)
    pub fn with_grids<R>(&mut self, grids: &[[f64; 3]], func: impl FnOnce(&mut Self) -> R) -> R {
        let old_grids = self.get_grids();
        self.try_clear_grids();
        self.set_grids(grids);
        let result = func(self);
        self.try_clear_grids();
        self.set_grids(&old_grids);
        result
    }

    /* #endregion */
}
