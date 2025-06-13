//! Temporary changes to `CInt` instance.

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

pub trait FakeMolForChargeArg {
    /// Create a fake molecule for gaussian distributed charges.
    ///
    /// You may also fake point charges by setting an extremely large exponent.
    fn fakemol_for_charges(coords: &[[f64; 3]], arg: Self) -> CInt;
}

/// Create a fake molecule for gaussian distributed charges.
///
/// You may also fake point charges by setting an extremely large exponent.
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

        CInt {
            atm: fake_atm,
            bas: fake_bas,
            ecpbas: vec![],
            env: fake_env,
            cint_type: CIntType::default(),
        }
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

        CInt {
            atm: fake_atm,
            bas: fake_bas,
            ecpbas: vec![],
            env: fake_env,
            cint_type: CIntType::default(),
        }
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

impl FakeMolForChargeArg for (&[f64], &[f64]) {
    fn fakemol_for_charges(coords: &[[f64; 3]], arg: Self) -> CInt {
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

        let (exponents, contr_coeffs) = arg;

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
        exponents.iter().zip(contr_coeffs.iter()).enumerate().for_each(
            |(i, (&exponent, &contr_coeff))| {
                let mut bas = [0; BAS_SLOTS];
                bas[ATOM_OF] = i as c_int;
                bas[NPRIM_OF] = 1;
                bas[NCTR_OF] = 1;
                bas[PTR_EXP] = ptr as c_int;
                bas[PTR_COEFF] = (ptr + 1) as c_int;
                let coef = contr_coeff / (2.0 * PI.sqrt() * gaussian_int(2.0, exponent));
                fake_bas.push(bas);
                fake_env.extend_from_slice(&[exponent, coef]);
                ptr += 2;
            },
        );

        CInt {
            atm: fake_atm,
            bas: fake_bas,
            ecpbas: vec![],
            env: fake_env,
            cint_type: CIntType::default(),
        }
    }
}

/* #endregion */

/// Usual cases for mutating `CInt` instance temporarily (something similar to
/// with-clause in Python).
impl CInt {
    /* #region cint_type */

    pub fn set_cint_type(&mut self, cint_type: impl Into<CIntType>) -> &mut Self {
        self.cint_type = cint_type.into();
        self
    }

    pub fn with_cint_type<R>(
        &mut self,
        cint_type: impl Into<CIntType>,
        func: impl FnOnce(&mut Self) -> R,
    ) -> R {
        let old_type = self.cint_type;
        self.set_cint_type(cint_type);
        let result = func(self);
        self.set_cint_type(old_type);
        result
    }

    /* #endregion */

    /* #region common_origin */

    pub fn get_common_origin(&self) -> [f64; 3] {
        const PTR_COMMON_ORIG: usize = cint_ffi::PTR_COMMON_ORIG as usize;
        self.env[PTR_COMMON_ORIG..PTR_COMMON_ORIG + 3].try_into().unwrap()
    }

    pub fn set_common_origin(&mut self, origin: [f64; 3]) -> &mut Self {
        const PTR_COMMON_ORIG: usize = cint_ffi::PTR_COMMON_ORIG as usize;
        self.env[PTR_COMMON_ORIG..PTR_COMMON_ORIG + 3].copy_from_slice(&origin);
        self
    }

    pub fn with_common_origin<R>(
        &mut self,
        origin: [f64; 3],
        func: impl FnOnce(&mut Self) -> R,
    ) -> R {
        let old_origin = self.get_common_origin();
        self.set_common_origin(origin);
        let result = func(self);
        self.set_common_origin(old_origin);
        result
    }

    /* #endregion */

    /* #region rinv_origin */

    pub fn get_rinv_origin(&self) -> [f64; 3] {
        const PTR_RINV_ORIG: usize = cint_ffi::PTR_RINV_ORIG as usize;
        self.env[PTR_RINV_ORIG..PTR_RINV_ORIG + 3].try_into().unwrap()
    }

    pub fn set_rinv_origin(&mut self, origin: [f64; 3]) -> &mut Self {
        const PTR_RINV_ORIG: usize = cint_ffi::PTR_RINV_ORIG as usize;
        self.env[PTR_RINV_ORIG..PTR_RINV_ORIG + 3].copy_from_slice(&origin);
        self
    }

    pub fn with_rinv_origin<R>(
        &mut self,
        origin: [f64; 3],
        func: impl FnOnce(&mut Self) -> R,
    ) -> R {
        let old_origin = self.get_rinv_origin();
        self.set_rinv_origin(origin);
        let result = func(self);
        self.set_rinv_origin(old_origin);
        result
    }

    /* #endregion */

    /* #region range_coulomb */

    pub fn get_range_coulomb(&self) -> f64 {
        const PTR_RANGE_OMEGA: usize = cint_ffi::PTR_RANGE_OMEGA as usize;
        self.env[PTR_RANGE_OMEGA]
    }

    pub fn set_range_coulomb(&mut self, range: f64) -> &mut Self {
        const PTR_RANGE_OMEGA: usize = cint_ffi::PTR_RANGE_OMEGA as usize;
        self.env[PTR_RANGE_OMEGA] = range;
        self
    }

    pub fn with_range_coulomb<R>(&mut self, range: f64, func: impl FnOnce(&mut Self) -> R) -> R {
        let old_range = self.get_range_coulomb();
        self.set_range_coulomb(range);
        let result = func(self);
        self.set_range_coulomb(old_range);
        result
    }

    pub fn omega(&self) -> f64 {
        self.get_range_coulomb()
    }

    pub fn set_omega(&mut self, omega: f64) -> &mut Self {
        self.set_range_coulomb(omega)
    }

    pub fn with_long_range_coulomb<R>(
        &mut self,
        omega: f64,
        func: impl FnOnce(&mut Self) -> R,
    ) -> R {
        self.with_range_coulomb(omega, func)
    }

    pub fn with_short_range_coulomb<R>(
        &mut self,
        omega: f64,
        func: impl FnOnce(&mut Self) -> R,
    ) -> R {
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

    pub fn get_rinv_zeta(&self) -> f64 {
        const PTR_RINV_ZETA: usize = cint_ffi::PTR_RINV_ZETA as usize;
        self.env[PTR_RINV_ZETA]
    }

    pub fn set_rinv_zeta(&mut self, zeta: f64) -> &mut Self {
        const PTR_RINV_ZETA: usize = cint_ffi::PTR_RINV_ZETA as usize;
        self.env[PTR_RINV_ZETA] = zeta;
        self
    }

    pub fn with_rinv_zeta<R>(&mut self, zeta: f64, func: impl FnOnce(&mut Self) -> R) -> R {
        let old_zeta = self.get_rinv_zeta();
        self.set_rinv_zeta(zeta);
        let result = func(self);
        self.set_rinv_zeta(old_zeta);
        result
    }

    /* #endregion */

    /* #region rinv_at_nucleus */

    pub fn set_rinv_at_nucleus(&mut self, atm_id: usize) -> &mut Self {
        const PTR_ZETA: usize = cint_ffi::PTR_ZETA as usize;
        let zeta = self.env[self.atm[atm_id][PTR_ZETA] as usize];
        let rinv = self.atom_coord(atm_id);
        self.set_rinv_origin(rinv);
        if zeta != 0.0 {
            self.set_rinv_zeta(zeta);
        }
        self
    }

    pub fn with_rinv_at_nucleus<R>(
        &mut self,
        atm_id: usize,
        func: impl FnOnce(&mut Self) -> R,
    ) -> R {
        let old_rinv = self.get_rinv_origin();
        let old_zeta = self.get_rinv_zeta();
        self.set_rinv_at_nucleus(atm_id);
        let result = func(self);
        self.set_rinv_origin(old_rinv);
        self.set_rinv_zeta(old_zeta);
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

    pub fn set_geom(&mut self, coords: &[[f64; 3]]) -> &mut Self {
        const PTR_COORD: usize = cint_ffi::PTR_COORD as usize;

        if self.natm() != coords.len() {
            panic!(
                "Number of coordinates ({}) does not match number of atoms ({})",
                coords.len(),
                self.natm()
            );
        }

        coords.iter().enumerate().for_each(|(i, coord)| {
            let ptr = self.atm[i][PTR_COORD] as usize;
            self.env[ptr..ptr + 3].copy_from_slice(coord);
        });
        self
    }

    pub fn with_geom<R>(&mut self, coords: &[[f64; 3]], func: impl FnOnce(&mut Self) -> R) -> R {
        let old_geom = self.get_geom();
        self.set_geom(coords);
        let result = func(self);
        self.set_geom(&old_geom);
        result
    }

    /* #endregion */
}
