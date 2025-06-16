fn main() {
    // read environment variables
    // - CINT_SRC: source of libcint (should be git repository URL or path to local
    //   source)
    // - CINT_VER: version of libcint (`v6.1.2` for example, should start with `v`)

    let cint_src = std::env::var("CINT_SRC").unwrap_or({
        if cfg!(feature = "qcint") { "https://github.com/sunqm/qcint.git".into() } else { "https://github.com/sunqm/libcint.git".into() }
    });
    let cint_ver = std::env::var("CINT_VER").unwrap_or("v6.1.2".into());

    if cfg!(feature = "build_from_source") {
        let flag_shared = if cfg!(feature = "static") { "OFF" } else { "ON" };
        let with_f12 = if cfg!(feature = "with_f12") { "ON" } else { "OFF" };
        let with_4c1e = if cfg!(feature = "with_4c1e") { "ON" } else { "OFF" };
        let dst = cmake::Config::new("external_deps")
            .define("BUILD_SHARED_LIBS", flag_shared)
            .define("CINT_SRC", cint_src)
            .define("CINT_VER", cint_ver)
            .define("WITH_F12", with_f12)
            .define("WITH_4C1E", with_4c1e)
            .build();
        println!("cargo:rustc-link-search=native={}/lib", dst.display());
    }
}
