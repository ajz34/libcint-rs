fn main() {
    // everything should have been prepared by libcint-src/build.rs
    println!("cargo:rustc-link-lib=static=cecp");
    if cfg!(feature = "static") {
        println!("cargo:rustc-link-lib=static=cint");
        println!("cargo:rustc-link-lib=quadmath");
    } else {
        println!("cargo:rustc-link-lib=cint");
    }
}
