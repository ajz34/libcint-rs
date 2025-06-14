use approx::assert_relative_eq;
use libcint::prelude::*;
use num::complex::ComplexFloat;

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test() {
        let mut cint_data = init_na_ecpso();

        cint_data.with_cint_type("sph", |cint_data| {
            let (out, _) = cint_data.integrate("ECPso", "s1", None).into();
            assert_relative_eq!(16.31443178892987, cint_fp(&out), epsilon = 1e-10);
        });

        cint_data.with_cint_type("cart", |cint_data| {
            let (out, _) = cint_data.integrate("ECPso", "s1", None).into();
            assert_relative_eq!(8.274184128427798, cint_fp(&out), epsilon = 1e-10);
        });

        cint_data.with_cint_type("spinor", |cint_data| {
            let (out, _) = cint_data.integrate_spinor("ECPso", "s1", None).into();
            assert_relative_eq!(0.0, cint_fp(&out).re(), epsilon = 1e-10);
            assert_relative_eq!(30.596605690744546, cint_fp(&out).im(), epsilon = 1e-10);
        });
    }

    fn init_na_ecpso() -> CInt {
        let atm = vec![[3, 20, 4, 23, 0, 0]];
        let bas = vec![[0, 0, 1, 1, 0, 24, 25, 0], [0, 1, 1, 1, 0, 26, 27, 0], [0, 1, 1, 1, 0, 28, 29, 0], [0, 2, 1, 1, 0, 30, 31, 0]];
        let env = vec![
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.,
            0.9448630622825309,
            0.9448630622825309,
            0.,
            0.,
            1.,
            2.526475110984259,
            4.,
            16.50286631934116,
            1.,
            2.917322170855303,
            1.,
            2.6093322745198853,
            0.,
            -3.,
            -3.,
            0.,
            -3.,
            -3.,
            0.,
            -3.,
            -3.,
            0.,
            -3.,
            -3.,
        ];
        let ecpbas = vec![
            [0, 0, 1, 1, 0, 32, 33, 0],
            [0, 0, 1, 1, 1, 32, 34, 0],
            [0, 1, 1, 1, 0, 35, 36, 0],
            [0, 1, 1, 1, 1, 35, 37, 0],
            [0, 2, 1, 1, 0, 38, 39, 0],
            [0, 2, 1, 1, 1, 38, 40, 0],
            [0, 3, 1, 1, 0, 41, 42, 0],
            [0, 3, 1, 1, 1, 41, 43, 0],
        ];

        CInt { atm, bas, ecpbas, env, cint_type: CIntType::Spheric }
    }
}
