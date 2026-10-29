fn main() {
    // everything should have been prepared by libcint-src/build.rs
    println!("cargo:rustc-link-lib=static=cecp");
    if cfg!(feature = "no_linkage") {
        return;
    }
    if cfg!(feature = "static") {
        println!("cargo:rustc-link-lib=static=cint");
        // It seems that quadmath is usually not linked by libcint on macOS (at
        // least for conda shipped versions).
        // Anyway, if encountered any problem, the API user should try dynamic
        // linking, feature `no_quadmath`, or manually write build.rs.
        // `cfg!(target_os)` in a build script is the host; the target is
        // `CARGO_CFG_TARGET_OS`.
        let target_os = std::env::var("CARGO_CFG_TARGET_OS").unwrap();
        if !cfg!(feature = "no_quadmath") && target_os != "macos" {
            println!("cargo:rustc-link-lib=quadmath");
        }
    } else {
        println!("cargo:rustc-link-lib=cint");
    }
    // if user specified environment variable CINT_LINK_QUADMATH, we will force
    // link
    if std::env::var("CINT_LINK_QUADMATH").is_ok() {
        println!("cargo:rustc-link-lib=quadmath");
    }
}
