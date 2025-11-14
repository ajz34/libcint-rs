use std::path::PathBuf;

fn build_libcint() {
    // read environment variables
    // - CINT_SRC: source of libcint (should be git repository URL or path to local
    //   source)
    // - CINT_VER: version of libcint (`v6.1.2` for example, should start with `v`)

    let cint_src = std::env::var("CINT_SRC").unwrap_or({
        if cfg!(feature = "qcint") {
            "https://github.com/sunqm/qcint.git".into()
        } else {
            "https://github.com/sunqm/libcint.git".into()
        }
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
        // CMAKE_INSTALL_LIBDIR can be lib64 on some platforms
        println!("cargo:rustc-link-search=native={}/lib", dst.display());
        println!("cargo:rustc-link-search=native={}/lib64", dst.display());
    }
}

/// Generate link search paths from a list of paths.
///
/// This allows paths like `/path/to/lib1:/path/to/lib2` to be split into
/// individual paths.
fn generate_link_search_paths(paths: &str) -> Vec<String> {
    let split_char = if cfg!(windows) { ";" } else { ":" };
    paths.split(split_char).map(|path| path.to_string()).collect()
}

/// Generate root candidates for library search paths.
///
/// Code modified from
///
/// https://github.com/coreylowman/cudarc/blob/main/build.rs
fn root_candidates(env_candidates: &[&str]) -> Vec<PathBuf> {
    let root_candidates = ["/usr", "/usr/local", "/usr/local/share", "/opt"];

    env_candidates
        .iter()
        .map(|p| p.to_string())
        .map(std::env::var)
        .filter_map(Result::ok)
        .flat_map(|path| generate_link_search_paths(&path))
        .filter(|path| !path.is_empty())
        .chain(root_candidates.into_iter().map(|p| p.to_string()))
        .map(|p| p.into())
        .collect()
}

/// Generate candidates for library search paths.
///
/// Code modified from
///
/// https://github.com/coreylowman/cudarc/blob/main/build.rs
fn lib_candidates() -> impl Iterator<Item = PathBuf> {
    let lib_candidates = [
        "",
        "lib",
        "lib/stubs",
        "lib/x64",
        "lib/Win32",
        "lib/x86_64",
        "lib/x86_64-linux-gnu",
        "lib64",
        "lib64/stubs",
        "targets/x86_64-linux",
        "targets/x86_64-linux/lib",
        "targets/x86_64-linux/lib/stubs",
    ];
    lib_candidates.into_iter().map(|p| p.into())
}

fn path_candidates(env_candidates: &[&str]) -> impl Iterator<Item = PathBuf> {
    root_candidates(env_candidates)
        .into_iter()
        .flat_map(|root| lib_candidates().map(move |lib| root.join(lib)))
        .filter(|path| path.exists())
        .map(|path| std::fs::canonicalize(path).unwrap())
}

fn link_cint() {
    let env_candidates = ["CINT_DIR", "REST_EXT_DIR", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "PATH"];
    // minimal rerun-if-env-changed to avoid unnecessary rebuilds
    println!("cargo:rerun-if-env-changed=CINT_DIR");
    for path in path_candidates(&env_candidates) {
        println!("cargo:rustc-link-search=native={}", path.display());
    }

    if cfg!(feature = "static") {
        println!("cargo:rustc-link-lib=static=cint");
        println!("cargo:rustc-link-lib=quadmath");
    } else {
        println!("cargo:rustc-link-lib=cint");
    }
}

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

fn main() {
    build_libcint();
    link_cint();
    build_ecp();
}
