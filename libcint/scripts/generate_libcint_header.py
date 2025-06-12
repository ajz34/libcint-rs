# # Auto-generate some parts of `libcint` rust interface

# This file is organized to be compatible to jupyter/jupytext, open as notebook.

# - This file must be manually modified (specify `dir_libcint_src`) to be executed.
# - Please let bindgen (`sudo apt install bindgen` or `cargo install bindgen-cli`) in bash `$PATH`.

# ## User Specifications

# -- libcint source code location --
# the same to `git clone https://github.com/sunqm/libcint`
dir_libcint_src = "/home/a/Software/source/libcint-6.1.2/"

# ## Preparasion

# ### python import

import os
import glob
import subprocess
import re

# +
with open(f"{dir_libcint_src}/include/cint.h.in", "r") as f:
    raw_cint_h_in = f.read()

with open(f"{dir_libcint_src}/include/cint_funcs.h", "r") as f:
    raw_cint_funcs_h = f.read()

with open(f"{dir_libcint_src}/CMakeLists.txt", "r") as f:
    raw_cmakelists = f.read()
# -

# ### obtain known intors

# -- configure list of integrals known to `cint_funcs.h` --
# known intor of `cint_funcs.h`
known_intor = []
for line in raw_cint_funcs_h.split("\n"):
    if line.startswith("extern CINTOptimizerFunction"):
        known_intor.append(line.strip().replace(";", "").split()[-1].replace("_optimizer", ""))
known_intor = sorted(set(known_intor))

# +
# -- configure list of integrals known to source code, autocode included --
# all intor generated in source code (autocode included)
# pattern: ALL_CINT(<intor>)
# pattern: ALL_CINT1E(<intor>)
# actual_intor = []
# for src_path in list(glob.iglob(f"{dir_libcint_src}/src/**/*.c", recursive=True)):
#     with open(src_path, "r") as f:
#         for line in f.readlines():
#             if line.startswith("ALL_CINT(") or line.startswith("ALL_CINT1E("):
#                 actual_intor.append(line.split("(")[-1].replace(")", "").replace(";", "").strip())
# actual_intor = sorted(set(actual_intor))

# +
# -- configure list of integrals known to source code, autocode included --
# all intor generated in source code (autocode included)
# pattern: 
#   CACHE_SIZE_T <intor>(...
#   FINT ng[] = {...};
# pattern rule:
#   <intor> ends with _sph, _cart, _spinor; no # in <intor>

actual_intor = {}
valid_suffix = ("_sph", "_cart", "_spinor")
for src_path in list(glob.iglob(f"{dir_libcint_src}/src/**/*.c", recursive=True)):
    with open(src_path, "r") as f:
        lines = f.readlines()
        for n, line in enumerate(lines):
            if line.startswith("CACHE_SIZE_T int"):
                # match pattern for intor
                intor = line.split("(")[0].split()[-1].strip()
                is_intor_valid = False
                for suffix in valid_suffix:
                    if intor.endswith(suffix):
                        intor = intor.replace(suffix, "")
                        is_intor_valid = True
                if "#" in intor or intor in actual_intor or not is_intor_valid:
                    continue
                # match pattern for ng
                ng = None
                for l in lines[n:n+5]:
                    if l.strip().startswith("FINT ng[]"):
                        ng = [int(i) for i in l.strip().replace(";", "").replace("}", "").split("{")[-1].split(",")]
                if ng is None:
                    raise ValueError(f"{intor}")
                actual_intor[intor] = ng
# -

# assert known intors are all included in all intors
assert len(set(known_intor) - set(actual_intor)) == 0


# ### Declare exceptions

def rule_with_f12(intor):
    return "_stg" in intor or "_yp" in intor


def rule_with_4c1e(intor):
    return "int4c1e" in intor


# ## `bindgen` auto-generation

# ### configure `cint.h`

# -- configure version of cint.h --
ver = {
    "cint_VERSION_MAJOR": None,
    "cint_VERSION_MINOR": None,
    "cint_VERSION_PATCH": None,
}
for line in raw_cmakelists.split("\n"):
    for key in ver:
        if line.startswith("set(" + key):
            ver[key] = int(line.strip().split()[1].replace(")", "").replace("\"", ""))
            break
cint_h = raw_cint_h_in \
    .replace("@cint_VERSION@", f"{ver['cint_VERSION_MAJOR']}.{ver['cint_VERSION_MINOR']}.{ver['cint_VERSION_PATCH']}") \
    .replace("@cint_SOVERSION@", f"{ver['cint_VERSION_MAJOR']}") \
    .replace("#cmakedefine I8", "/* #undef I8 */") \
    .replace("#cmakedefine CACHE_SIZE_I8", "/* #undef CACHE_SIZE_I8 */")

# ### configure `cint_funcs.h`

# necessary header (modified from libcint 6.1.1)
cint_funcs = """
#include "cint.h"

#if !defined HAVE_DEFINED_CINTINTEGRALFUNCTION
#define HAVE_DEFINED_CINTINTEGRALFUNCTION
typedef void CINTOptimizerFunction(
            CINTOpt **opt,
            const FINT *atm, FINT natm, const FINT *bas, FINT nbas, const double *env);
typedef CACHE_SIZE_T CINTIntegralFunctionReal(
            double *out, const FINT *dims, const FINT *shls,
            const FINT *atm, FINT natm, const FINT *bas, FINT nbas, const double *env,
            const CINTOpt *opt, double *cache);
typedef CACHE_SIZE_T CINTIntegralFunctionComplex(
            double complex *out, const FINT *dims, const FINT *shls,
            const FINT *atm, FINT natm, const FINT *bas, FINT nbas, const double *env,
            const CINTOpt *opt, double *cache);
#endif
"""

# for unknown intors (mostly not in autocode, but in main source code), generate signature for header
# exclusion (F12): _stg_, _yp_
token = """
extern CINTOptimizerFunction       {0:}_optimizer;
extern CINTIntegralFunctionReal    {0:}_cart;
extern CINTIntegralFunctionReal    {0:}_sph;
extern CINTIntegralFunctionComplex {0:}_spinor;
"""
cint_funcs = cint_funcs + "".join([
    token.format(intor)
    for intor in sorted(set(actual_intor))
    if not rule_with_f12(intor)])
cint_funcs = cint_funcs.replace("#include <cint.h>", "#include \"cint.h\"")

# for unknown intors (mostly not in autocode, but in main source code), generate signature for header
# exclusion (F12): _stg_, _yp_; they do not have cart and spinor
token = """
extern CINTOptimizerFunction       {0:}_optimizer;
extern CINTIntegralFunctionReal    {0:}_sph;
"""
cint_funcs = cint_funcs + "".join([
    token.format(intor)
    for intor in sorted(set(actual_intor))
    if rule_with_f12(intor)])
cint_funcs = cint_funcs.replace("#include <cint.h>", "#include \"cint.h\"")

# +
# libcint2 interface compability
# token = """
# extern CINTOptimizerFunction       c{0:}_optimizer;
# extern CINTIntegralFunctionReal    c{0:}_cart;
# extern CINTIntegralFunctionReal    c{0:}_sph;
# extern CINTIntegralFunctionComplex c{0:}_spinor;
# """
# cint_funcs = cint_funcs + "".join([token.format(intor) for intor in sorted(actual_intor) if intor not in ["int2e"]])
# -

# ### Merge `cint.h` and `cint_funcs.h`

cint_funcs = cint_funcs.replace('#include "cint.h"', cint_h)

with open("cint_funcs.h", "w") as f:
    f.write(cint_funcs)

# ### `bindgen` transpile header to rust

subprocess.run([
    "bindgen", "cint_funcs.h",
    "-o", "cint_ffi.rs",
    "--allowlist-file", "cint_funcs.h",
    "--no-layout-tests",
    "--merge-extern-blocks",
])

# ### modify bindgen auto-generated file for our purpose

with open("cint_ffi.rs", "r") as f:
    token = f.read()

# exclude some function of complex.c
res = re.findall(r"(?:\#.*?\})", token.replace("\n", "NEWLINE"))
for r in res:
    if "pub fn c" in r and "pub fn cint" not in r:
        token = token.replace(r.replace("NEWLINE", "\n") + "\n", "")

# +
# clean some `::std::` code
token = """
use core::ffi::c_int;

""" + token

token = token.replace("::std::os::raw::c_int", "c_int")
token = token.replace("::std::option::Option", "Option")
# -

with open("cint_ffi.rs", "w") as f:
    f.write(token)

# ### cleanup

subprocess.run(["mv", "cint_ffi.rs", "../src/ffi/cint_ffi.rs"])
subprocess.run(["rm", "cint_funcs.h"])

# ## Auto-generate `cint_wrapper.rs`

# ### Integrator

# +
# initial information

token_wrapper = """
/* Generated by python scripts for libcint low-level wrapper. */
/* Should not modify manually. */

use crate::cint::CIntKind;
use crate::ffi::cint_ffi::*;
use crate::ffi::wrapper_traits::Integrator;
use crate::impl_integrator;
use core::any::Any;
use core::ffi::{c_int, c_void};

#[allow(unused_imports)]
use crate::ffi::wrapper_traits::{panic_cart, panic_spinor};

"""


# +
# Integrator: macro expansion

def gen_impl_integrator(intor):
    token = f"""impl_integrator!(
    {0:},
    {0:}_optimizer,
    {0:}_sph,
    {0:}_cart,
    {0:}_spinor,
    {5:}, {6:}, {7:}
    {2:}, {3:}, {4:}, vec!{5:},
    "{1:}", "{0:}",
    CIntKind::Int);\n"""
    integrator_category = intor.split("_")[0]

    # functions and availability
    is_sph_available = "true"
    is_cart_available = "true"
    is_spinor_available = "true"
    func_integral_sph = f"{intor}_sph"
    func_integral_cart = f"{intor}_cart"
    func_integral_spinor = f"{intor}_spinor"
    if rule_with_f12(intor):
        is_cart_available = "false"
        is_spinor_available = "false"
        func_integral_cart = "panic_cart"
        func_integral_spinor = "panic_spinor"

    # counting numbers
    intor_center = intor[:5]
    if intor_center == "int2e":
        intor_center = "int4c"
    if intor_center == "int1e":
        intor_center = "int2c"
    n_center = int(intor_center[3])
    ng = actual_intor[intor]
    assert len(ng) == 8
    comp_1e, comp_2e, comp_tensor = ng[-3:]
    comp_all = max(comp_1e, 1) * max(comp_2e, 1) * comp_tensor

    # write token
    token = f"""impl_integrator!(
        {intor},
        {intor}_optimizer,
        {func_integral_sph},
        {func_integral_cart},
        {func_integral_spinor},
        {is_sph_available},
        {is_cart_available},
        {is_spinor_available},
        {comp_all}, {comp_tensor}, {n_center},
        vec!{ng},
        "{integrator_category}", "{intor}",
        CIntKind::Int);\n"""

    # cargo features for control
    if rule_with_f12(intor):
        token = '#[cfg(feature = "with_f12")]\n' + token
    if rule_with_4c1e(intor):
        token = '#[cfg(feature = "with_4c1e")]\n' + token
        
    return token

for intor in actual_intor:
    token_wrapper += gen_impl_integrator(intor)

# +
# get_integrator

token_wrapper += """

pub fn get_cint_integrator(name: &str) -> Option<Box<dyn Integrator>> {
    match name {
"""

for intor in actual_intor:
    if rule_with_f12(intor):
        token_wrapper += '#[cfg(feature = "with_f12")]\n'
    if rule_with_4c1e(intor):
        token_wrapper += '#[cfg(feature = "with_4c1e")]\n'
        
    token_wrapper += f"""
        "{intor}" => Some(Box::new({intor})),
    """.strip()

token_wrapper += """
        _ => None,
    }
}
"""
# -

with open("../src/ffi/cint_wrapper.rs", "w") as f:
    f.write(token_wrapper)

# ## Finalize with format

os.chdir("..")

subprocess.run(["cargo", "fmt"])
