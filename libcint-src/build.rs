use std::path::PathBuf;

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

fn main() {
    if cfg!(feature = "build_from_source") {
        let flag_shared = if cfg!(feature = "static") { "OFF" } else { "ON" };
        let dst = cmake::Config::new("external_deps").define("BUILD_SHARED_LIBS", flag_shared).build();
        println!("cargo:rustc-link-search=native={}/lib", dst.display());
    } else {
        let env_candidates = ["CINT_DIR", "REST_EXT_DIR", "LD_LIBRARY_PATH", "DYLD_LIBRARY_PATH", "PATH"];
        for env_var in env_candidates.iter() {
            if let Ok(path) = std::env::var(env_var) {
                println!("cargo:rerun-if-env-changed={env_var}");
                println!("cargo:rerun-if-changed={path}");
            }
        }
        for path in path_candidates(&env_candidates) {
            println!("cargo:rustc-link-search=native={}", path.display());
        }
    }

    if cfg!(feature = "static") {
        println!("cargo:rustc-link-lib=static=cint");
        println!("cargo:rustc-link-lib=quadmath");
    } else {
        println!("cargo:rustc-link-lib=cint");
    }
}
