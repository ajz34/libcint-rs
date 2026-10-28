//! Reuse of libcint optimizers ([`CIntOpt`]) across integral calls.
//!
//! These tests cover the contract of `IntegrateArgs::cintopt` /
//! `IntorCrossArgs::cintopt`: a handle built once by [`CInt::get_cintopt`] must
//! give exactly the same integrals as the implicit per-call build, must be
//! reusable across shell batches, and must be rejected when it cannot match the
//! integrator, the basis size, or (for `int2c2e`) the angular momentum.

use approx::assert_relative_eq;
use libcint::prelude::*;

#[cfg(test)]
mod test {
    use super::*;

    /// Integrals with a reused optimizer must equal the implicit-build ones.
    #[test]
    fn reuse_int2c2e() {
        let mut data = init_h2o_def2_tzvp();
        let opt = data.get_cintopt("int2c2e");
        assert_eq!(opt.intor(), "int2c2e");
        assert_eq!(opt.nbas(), data.nbas());

        let (auto, shape_auto) = data.integrate("int2c2e", "s1", None).into();

        for _ in 0..2 {
            let args = data.integrate_args_builder().intor("int2c2e").cintopt(&opt).build().unwrap();
            let (reused, shape) = data.integrate_with_args(args).into();
            assert_eq!(shape, shape_auto);
            assert_eq!(reused, auto);
        }

        // the handle also works for the cartesian variant
        let cart = data.with_cint_type("cart", |data| {
            let args = data.integrate_args_builder().intor("int2c2e").cintopt(&opt).build().unwrap();
            let (out, _) = data.integrate_with_args(args).into();
            out
        });
        let cart_ref = data.with_cint_type("cart", |data| {
            let (out, _) = data.integrate("int2c2e", "s1", None).into();
            out
        });
        assert_eq!(cart, cart_ref);
    }

    /// The memory-bound J/K pattern: one optimizer, many aux-basis batches.
    #[test]
    fn reuse_int3c2e_across_batches() {
        let h2o = init_h2o_def2_tzvp();
        let jk = init_h2o_def2_jk();
        let merged = &h2o + &jk;
        let nbas_bra = h2o.nbas();
        let nbas = merged.nbas();
        let opt = merged.get_cintopt("int3c2e");
        assert_eq!(opt.nbas(), nbas);

        let mid = nbas_bra + (nbas - nbas_bra) / 2;
        for aux in [[nbas_bra, mid], [mid, nbas]] {
            let slices = [[0, nbas_bra], [0, nbas_bra], aux];
            let (auto, _) = merged.integrate("int3c2e", "s2ij", slices).into();
            let args = merged.integrate_args_builder().intor("int3c2e").aosym("s2ij").shls_slice(&slices).cintopt(&opt).build().unwrap();
            let (reused, _) = merged.integrate_with_args(args).into();
            assert_eq!(reused, auto);
        }
    }

    /// Derivatives have their own optimizer data (`ng` differs) and must not be
    /// mixed up.
    #[test]
    fn reuse_int3c2e_derivatives() {
        let h2o = init_h2o_def2_tzvp();
        let jk = init_h2o_def2_jk();
        let merged = &h2o + &jk;
        let nbas_bra = h2o.nbas();
        let nbas = merged.nbas();
        let slices = [[0, nbas_bra], [0, nbas_bra], [nbas_bra, nbas]];

        for intor in ["int3c2e_ip1", "int3c2e_ip2"] {
            let opt = merged.get_cintopt(intor);
            assert_eq!(opt.intor(), intor);
            let (auto, _) = merged.integrate(intor, "s2ij", slices).into();
            let args = merged.integrate_args_builder().intor(intor).aosym("s2ij").shls_slice(&slices).cintopt(&opt).build().unwrap();
            let (reused, _) = merged.integrate_with_args(args).into();
            assert_eq!(reused, auto);
        }
    }

    #[test]
    fn reuse_int2e_and_derivative() {
        let data = init_h2o_def2_tzvp();
        let slices = [[0, 4], [0, 4], [0, 4], [0, 4]];

        for intor in ["int2e", "int2e_ip1"] {
            let opt = data.get_cintopt(intor);
            let (auto, _) = data.integrate(intor, "s1", slices).into();
            let args = data.integrate_args_builder().intor(intor).shls_slice(&slices).cintopt(&opt).build().unwrap();
            let (reused, _) = data.integrate_with_args(args).into();
            assert_eq!(reused, auto);
        }
    }

    #[test]
    fn reuse_int2c2e_spinor() {
        let data = init_h2o_def2_tzvp();
        // suffix is accepted by `get_cintopt` as well as by `integrate_spinor`
        let opt = data.get_cintopt("int2c2e_spinor");
        let slices = [[0, 4], [0, 4]];
        let (auto, shape) = data.integrate_spinor("int2c2e_spinor", "s1", slices).into();
        let args = data.integrate_args_builder_spinor().intor("int2c2e_spinor").shls_slice(&slices).cintopt(&opt).build().unwrap();
        let (reused, shape2) = data.integrate_with_args_spinor(args).into();
        assert_eq!(shape, shape2);
        assert_eq!(reused, auto);
    }

    #[test]
    fn reuse_ecp_scalar() {
        let data = init_sb2me4_cc_pvtz();
        let opt = data.get_cintopt("ECPscalar");
        assert_eq!(opt.intor(), "ECPscalar");

        let (auto, shape) = data.integrate("ECPscalar", "s1", None).into();
        let args = data.integrate_args_builder().intor("ECPscalar").cintopt(&opt).build().unwrap();
        let (reused, shape2) = data.integrate_with_args(args).into();
        assert_eq!(shape, shape2);
        assert_eq!(reused, auto);
    }

    /// Three-center cross integral with a bra-only optimizer (PySCF `aux_e2`
    /// pattern): `int3c2e` only reads optimizer data for the first two centres,
    /// so a prefix optimizer is valid.
    #[test]
    fn reuse_cross_int3c2e_with_bra_opt() {
        let h2o = init_h2o_def2_tzvp();
        let jk = init_h2o_def2_jk();
        let opt = h2o.get_cintopt("int3c2e");
        assert_eq!(opt.nbas(), h2o.nbas());

        let mols = [&h2o, &h2o, &jk];
        let slices = [[0, h2o.nbas()], [0, h2o.nbas()], [0, jk.nbas()]];
        let (auto, shape) = CInt::integrate_cross("int3c2e", mols, "s2ij", slices).into();

        let args = h2o
            .integrate_cross_args_builder()
            .intor("int3c2e")
            .mols(&mols)
            .aosym("s2ij")
            .shls_slice(&slices)
            .cintopt(&opt)
            .build()
            .unwrap();
        let (reused, shape2) = CInt::integrate_cross_with_args(args).into();
        assert_eq!(shape, shape2);
        assert_eq!(reused, auto);
    }

    #[test]
    fn mismatched_intor_is_rejected() {
        let data = init_h2o_def2_tzvp();
        let opt = data.get_cintopt("int2c2e");
        let args = data.integrate_args_builder().intor("int3c2e").cintopt(&opt).build().unwrap();
        let err = match data.integrate_with_args_f(args) {
            Ok(_) => panic!("expected InvalidValue"),
            Err(err) => err,
        };
        assert!(matches!(err.kind, CIntErrorBase::InvalidValue(_)), "{err:?}");
    }

    #[test]
    fn insufficient_nbas_is_rejected() {
        let small = init_h2o_def2_tzvp();
        let big = init_c10h22_def2_qzvp();
        let opt = small.get_cintopt("int2e");
        let slices = [[0, 100], [0, 100], [0, 100], [0, 100]];
        let args = big.integrate_args_builder().intor("int2e").shls_slice(&slices).cintopt(&opt).build().unwrap();
        let err = match big.integrate_with_args_f(args) {
            Ok(_) => panic!("expected InvalidValue"),
            Err(err) => err,
        };
        assert!(matches!(err.kind, CIntErrorBase::InvalidValue(_)), "{err:?}");
    }

    /// `int2c2e` has no on-the-fly fallback for its angular-momentum index
    /// table, so an optimizer built with too small a maximum angular momentum
    /// must be rejected instead of being passed to C.
    #[test]
    fn int2c2e_insufficient_max_l_is_rejected() {
        let s_only = fakemol_for_charges(&[[0.0, 0.0, 0.0]], 1.0);
        let data = init_h2o_def2_tzvp();
        let opt = s_only.get_cintopt("int2c2e");
        let args = data.integrate_args_builder().intor("int2c2e").cintopt(&opt).build().unwrap();
        let err = match data.integrate_with_args_f(args) {
            Ok(_) => panic!("expected InvalidValue"),
            Err(err) => err,
        };
        assert!(matches!(err.kind, CIntErrorBase::InvalidValue(_)), "{err:?}");
    }

    /// Families where libcint ignores the optimizer take the NULL path now;
    /// values must be unchanged.
    #[test]
    fn ignored_families_still_work() {
        let data = init_h2o_def2_tzvp();
        let (out, shape) = data.integrate("int1e_kin", "s1", [[0, 12], [8, 18]]).into();
        assert_eq!(shape, [32, 26]);
        assert_relative_eq!(cint_fp(&out), 0.8338209486794069, epsilon = 1e-10);
    }

    #[test]
    #[ignore = "timing only; run with `cargo test --release -- --ignored --nocapture`"]
    fn optimizer_reuse_bench() {
        let h2o = init_h2o_def2_tzvp();
        let jk = init_h2o_def2_jk();
        let merged = &h2o + &jk;
        let nbas_bra = h2o.nbas();
        let nbas = merged.nbas();
        let nbatch = 8;
        let batches: Vec<[usize; 2]> = (0..nbatch)
            .map(|b| {
                let s = nbas_bra + (nbas - nbas_bra) * b / nbatch;
                let e = nbas_bra + (nbas - nbas_bra) * (b + 1) / nbatch;
                [s, e]
            })
            .collect();

        for intor in ["int3c2e_ip1", "int3c2e_ip2"] {
            let t = std::time::Instant::now();
            for aux in &batches {
                let slices = [[0, nbas_bra], [0, nbas_bra], *aux];
                let _ = merged.integrate(intor, "s2ij", slices);
            }
            let auto = t.elapsed().as_secs_f64();

            let opt = merged.get_cintopt(intor);
            let t = std::time::Instant::now();
            for aux in &batches {
                let slices = [[0, nbas_bra], [0, nbas_bra], *aux];
                let args = merged.integrate_args_builder().intor(intor).aosym("s2ij").shls_slice(&slices).cintopt(&opt).build().unwrap();
                let _ = merged.integrate_with_args(args);
            }
            let reused = t.elapsed().as_secs_f64();
            println!("{intor}: per-call build {auto:.3}s, reused opt {reused:.3}s, ratio {:.2}", auto / reused);
            assert!(reused < auto);
        }
    }
}
