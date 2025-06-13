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
}
