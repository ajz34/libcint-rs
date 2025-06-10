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
    "-o", "cecp_ffi.rs",
    "--allowlist-file", "../src/cecp/nr_ecp.h",
    "--no-layout-tests",
    "--merge-extern-blocks",
])

with open("cecp_ffi.rs", "r") as f:
    token = f.read()

# change mutability
token = token \
    .replace("dims: *mut", "dims: *const") \
    .replace("shls: *mut", "shls: *const") \
    .replace("atm: *mut", "atm: *const") \
    .replace("bas: *mut", "bas: *const") \
    .replace("env: *mut", "env: *const") \
    .replace("opt: *mut ECPOpt", "opt: *const ECPOpt")

# +
# clean some `::std::` code
token = """
use core::ffi::{c_int, c_ulong};

""" + token

token = token.replace("::std::os::raw::c_int", "c_int")
token = token.replace("::std::os::raw::c_ulong", "c_ulong")
token = token.replace("::std::option::Option", "Option")
# -

with open("cecp_ffi.rs", "w") as f:
    f.write(token)

subprocess.run(["mv", "cecp_ffi.rs", "../src/ffi/cecp_ffi.rs"])

# ## Finalize with format

os.chdir("..")

subprocess.run(["cargo", "fmt"])
