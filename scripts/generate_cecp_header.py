# # Auto-generate some parts of `gto/nr_ecp` rust interface

# - Please let bindgen (`sudo apt install bindgen` or `cargo install bindgen-cli`) in bash `$PATH`.

# ## Preparasion

import os
import glob
import subprocess
import re

# ## `bindgen` auto-generation

with open("../src/cecp/nr_ecp.h", "r") as f:
    header = f.read()

subprocess.run([
    "bindgen", "../src/cecp/nr_ecp.h",
    "-o", "cecp.rs",
    "--no-layout-tests",  # disable redundant tests
    "--wasm-import-module-name=cint",  # to be modified to #[link(name = ...)]
    "--blocklist-function", "[_].*",  # exclude some function of complex.c
    "--blocklist-var", "_.*",  # exclude variables from complex.c and stdint.c
    "--blocklist-type", "_.*",  # exclude types from complex.c
])

with open("cecp.rs", "r") as f:
    token = f.read()

# exclude some function of complex.c
res = re.findall(r"(?:\#.*?\})", token.replace("\n", "NEWLINE"))
for r in res:
    if "pub fn c" in r and "pub fn cint" not in r:
        token = token.replace(r.replace("NEWLINE", "\n") + "\n", "")

# change to #[link(name = "cint")]
token = token.replace("wasm_import_module", "name")

# change mutability
token = token \
    .replace("dims: *mut", "dims: *const") \
    .replace("shls: *mut", "shls: *const") \
    .replace("atm: *mut", "atm: *const") \
    .replace("bas: *mut", "bas: *const") \
    .replace("env: *mut", "env: *const") \
    .replace("opt: *mut ECPOpt", "opt: *const ECPOpt")

with open("cecp.rs", "w") as f:
    f.write(token)

subprocess.run(["mv", "cecp.rs", "../src/"])
