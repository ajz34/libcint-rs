fn main() {
    build_cint();
    build_ecp();
}

/// Builds the `libcint` library from the provided path.
fn build_cint() {
    // let dst = cmake::Config::new("external_deps").build();
    // println!("cargo:rustc-link-search=native={}/lib", dst.display());
    // println!("cargo:rustc-link-lib=cint");
    // println!("cargo:rustc-link-lib=quadmath");
}

/// Builds the `cecp` code from PySCF.
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
