//! Temporary changes to `CInt` instance.

use crate::prelude::*;

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
