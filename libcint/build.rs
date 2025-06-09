use std::path::PathBuf;

fn build_ecp() {
    cc::Build::new()
        .file("src/cecp/f2c_dgemm.c")
        .file("src/cecp/nr_ecp.c")
        .file("src/cecp/nr_ecp_deriv.c")
        .flag_if_supported("-Wno-unused-parameter")
        .flag_if_supported("-Wno-implicit-function-declaration")
        .flag_if_supported("-Wno-parentheses")
        .compile("cecp");
    println!("cargo::rerun-if-changed=src/cecp");
}

/// This build script is only for testing purposes.
///
/// It is not used in the actual library.
/// To make this library work, you either
/// - link `cint` in `build.rs` of your own project
/// - try `libcint-src` crate (which may also requires some configurations)
fn dev_link_cint() {
    println!("cargo:rerun-if-env-changed=CINT_DEV");
    if std::env::var("CINT_DEV").is_ok() {
        std::env::var("LD_LIBRARY_PATH")
            .unwrap()
            .split(":")
            .filter(|path| !path.is_empty())
            .filter(|path| PathBuf::from(path).exists())
            .for_each(|path| {
                println!("cargo:rustc-link-search=native={path}");
            });
        println!("cargo:rustc-link-lib=cint");
    }
}

fn main() {
    build_ecp();
    dev_link_cint();
}
