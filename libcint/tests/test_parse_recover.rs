//! Tests for verifying parsed molecules produce correct integrals.
//!
//! These tests verify that parsing molecules with `CIntMolInput` produces
//! CInt instances that give correct integral values.

#[cfg(feature = "bse")]
mod test_parse_recover {
    use libcint::{
        parse::mol::CIntMolInputBuilder,
        prelude::{init_c10h22_def2_qzvp, init_h2o_def2_jk, init_h2o_def2_tzvp, init_sb2me4_cc_pvtz},
    };
    use std::collections::HashMap;

    #[test]
    fn test_recover_h2o_def2_tzvp() {
        let mol_input = CIntMolInputBuilder::default().atom("O; H 1 0.94; H 1 0.94 2 104.5").basis("def2-TZVP").build().unwrap();
        let mol = mol_input.create_mol();

        let cint_ref = init_h2o_def2_tzvp();
        // perform check of similarity of CInt instances
        assert_eq!(mol.cint.atm, cint_ref.atm);
        assert_eq!(mol.cint.bas, cint_ref.bas);
        assert_eq!(mol.cint.ecpbas, cint_ref.ecpbas);
        assert_eq!(mol.cint.env.len(), cint_ref.env.len());
        let env_diff = mol.cint.env.iter().zip(cint_ref.env.iter()).map(|(a, b)| (a - b).abs()).fold(0.0, |acc, x| acc + x);
        assert!(env_diff < 1e-5, "CInt env differs from reference by {}", env_diff);
    }

    /// Test H2O with def2-universal-jkfit basis.
    #[test]
    fn test_recover_h2o_def2_jk() {
        let mol = CIntMolInputBuilder::default()
            .atom("O; H 1 0.94; H 1 0.94 2 104.5")
            .basis("def2-universal-jkfit")
            .build()
            .unwrap()
            .create_mol();

        let cint_ref = init_h2o_def2_jk();
        // perform check of similarity of CInt instances
        assert_eq!(mol.cint.atm, cint_ref.atm);
        assert_eq!(mol.cint.bas, cint_ref.bas);
        assert_eq!(mol.cint.ecpbas, cint_ref.ecpbas);
        assert_eq!(mol.cint.env.len(), cint_ref.env.len());
        let env_diff = mol.cint.env.iter().zip(cint_ref.env.iter()).map(|(a, b)| (a - b).abs()).fold(0.0, |acc, x| acc + x);
        assert!(env_diff < 1e-5, "CInt env differs from reference by {}", env_diff);
    }

    /// Test Sb2Me4 with ECP (cc-pVTZ-PP for Sb).
    #[test]
    fn test_recover_sb2me4_cc_pvtz() {
        let mut basis_map = HashMap::new();
        basis_map.insert("C".to_string(), "cc-pVTZ".to_string());
        basis_map.insert("H".to_string(), "cc-pVTZ".to_string());
        basis_map.insert("Sb".to_string(), "cc-pVTZ-PP".to_string());

        let atom_str = r#"
            Sb        -1.33937843      0.44597852     -1.27279684
            Sb         1.33937843     -0.44597852     -1.27279684
            C         -1.40429524      1.10441871      0.83468205
            C         -2.16210130     -1.56132398     -0.84717555
            C          2.16210130      1.56132398     -0.84717555
            C          1.40429524     -1.10441871      0.83468205
            H         -0.69918639      1.91987631      1.00872018
            H         -1.16111477      0.29030616      1.51873028
            H         -2.40124532      1.47235562      1.08516843
            H         -2.02002046     -2.22909286     -1.69887295
            H         -1.69052287     -2.01612927      0.02577778
            H         -3.23450854     -1.49489801     -0.65423339
            H          2.02002046      2.22909286     -1.69887295
            H          3.23450854      1.49489801     -0.65423339
            H          1.69052287      2.01612927      0.02577778
            H          0.69918639     -1.91987631      1.00872018
            H          2.40124532     -1.47235562      1.08516843
            H          1.16111477     -0.29030616      1.51873028
        "#;

        let mol = CIntMolInputBuilder::default().atom(atom_str).basis(basis_map).build().unwrap().create_mol();

        let cint_ref = init_sb2me4_cc_pvtz();
        // perform check of similarity of CInt instances
        assert_eq!(mol.cint.atm, cint_ref.atm);
        assert_eq!(mol.cint.bas, cint_ref.bas);
        assert_eq!(mol.cint.ecpbas, cint_ref.ecpbas);
        assert_eq!(mol.cint.env.len(), cint_ref.env.len());
        let env_diff = mol.cint.env.iter().zip(cint_ref.env.iter()).map(|(a, b)| (a - b).abs()).fold(0.0, |acc, x| acc + x);
        assert!(env_diff < 1e-5, "CInt env differs from reference by {}", env_diff);
    }

    /// Test C10H22 with def2-QZVP basis.
    #[test]
    fn test_recover_c10h22_def2_qzvp() {
        let atom_str = r#"
            C          0.99590        0.00874        0.02912
            C          2.51497        0.01491        0.04092
            C          3.05233        0.96529        1.10755
            C          4.57887        0.97424        1.12090
            C          5.10652        1.92605        2.19141
            C          6.63201        1.93944        2.20819
            C          7.15566        2.89075        3.28062
            C          8.68124        2.90423        3.29701
            C          9.20897        3.85356        4.36970
            C         10.73527        3.86292        4.38316
            C         11.27347        4.81020        5.45174
            C         12.79282        4.81703        5.46246
            H          0.62420       -0.67624       -0.73886
            H          0.60223        1.00743       -0.18569
            H          0.59837       -0.31563        0.99622
            H          2.88337        0.31565       -0.94674
            H          2.87961       -1.00160        0.22902
            H          2.67826        0.66303        2.09347
            H          2.68091        1.98005        0.91889
            H          4.95608        1.27969        0.13728
            H          4.95349       -0.03891        1.31112
            H          4.72996        1.62033        3.17524
            H          4.73142        2.93925        2.00202
            H          7.00963        2.24697        1.22537
            H          7.00844        0.92672        2.39728
            H          6.77826        2.58280        4.26344
            H          6.77905        3.90354        3.09209
            H          9.05732        3.21220        2.31361
            H          9.05648        1.89051        3.48373
            H          8.83233        3.54509        5.35255
            H          8.83406        4.86718        4.18229
            H         11.10909        4.16750        3.39797
            H         11.10701        2.84786        4.56911
            H         10.90576        4.50686        6.43902
            H         10.90834        5.82716        5.26681
            H         13.18649        3.81699        5.67312
            H         13.16432        5.49863        6.23161
            H         13.18931        5.14445        4.49536
        "#;

        let mol = CIntMolInputBuilder::default().atom(atom_str).basis("def2-QZVP").build().unwrap().create_mol();

        let cint_ref = init_c10h22_def2_qzvp();
        // perform check of similarity of CInt instances
        assert_eq!(mol.cint.atm, cint_ref.atm);
        assert_eq!(mol.cint.bas, cint_ref.bas);
        assert_eq!(mol.cint.ecpbas, cint_ref.ecpbas);
        assert_eq!(mol.cint.env.len(), cint_ref.env.len());
        let env_diff = mol.cint.env.iter().zip(cint_ref.env.iter()).map(|(a, b)| (a - b).abs()).fold(0.0, |acc, x| acc + x);
        assert!(env_diff < 1e-5, "CInt env differs from reference by {}", env_diff);
    }
}
