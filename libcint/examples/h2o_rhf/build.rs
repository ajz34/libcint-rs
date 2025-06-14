use std::path::PathBuf;

/// This build script is only for testing purposes.
///
/// It is not used in the actual library.
/// To make this library work, you either
/// - link `cint` in `build.rs` of your own project
/// - try `libcint-src` crate (which may also requires some configurations)
fn dev_link() {
    println!("cargo:rerun-if-env-changed=CINT_DEV");
    if std::env::var("CINT_DEV").is_ok() {
        std::env::var("LD_LIBRARY_PATH").unwrap().split(":").filter(|path| !path.is_empty()).filter(|path| PathBuf::from(path).exists()).for_each(
            |path| {
                println!("cargo:rustc-link-search=native={path}");
            },
        );
        println!("cargo:rustc-link-lib=cint");
        println!("cargo:rustc-link-lib=openblas");
    }
}

fn main() {
    dev_link();
}
