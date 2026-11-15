# # Auto-generate some parts of `libcint` rust interface

# This file is organized to be compatible to jupyter/jupytext, open as notebook.

# - This file must be manually modified (specify `dir_libcint_src`) to be executed.
# - Please let bindgen (`sudo apt install bindgen` or `cargo install bindgen-cli`) in bash `$PATH`.

import os
import glob
import subprocess
import re

# ## User Specifications

# -- libcint source code location --
# the same to `git clone https://github.com/sunqm/libcint`
dir_libcint_src = "/home/a/Software/source/libcint-6.1.2/"
dir_qcint_src = "/home/a/Git-Others/qcint/"

# ## Preparasion (qcint)

# +
# only certain files are compiled by qcint, defined by CMakeLists
with open(f"{dir_qcint_src}/CMakeLists.txt", "r") as f:
    raw_cmakelists = f.read()

lines = raw_cmakelists.split("\n")
for n, l in enumerate(lines):
    if l.startswith("set(cintSrc"):
        break
files = []
for l in lines[n+1:]:
    if l.strip() == "":
        break
    files += l.strip().replace(")", "").split()
qcint_src_paths = [f"{dir_qcint_src}/{f}" for f in files]

# +
# -- configure list of integrals known to source code, autocode included --
# all intor generated in source code (autocode included)
# pattern:
#   CACHE_SIZE_T <intor>(...
#   FINT ng[] = {...};
# pattern rule:
#   <intor> ends with _sph, _cart, _spinor; no # in <intor>
#
# `BREIT0(X, ncomp_tensor)` in breit.c expands to
# `int2e_breit_X_spinor` and `int2e_breit_X_optimizer`; these names never
# appear literally in the source, so add them after the direct scan.

def scan_intor_defs(src_paths):
    actual = {}
    valid_suffix = ("_sph", "_cart", "_spinor")
    texts = []
    for src_path in src_paths:
        with open(src_path, "r") as f:
            text = f.read()
        texts.append(text)
        lines = text.splitlines()
        for n, line in enumerate(lines):
            if line.startswith("CACHE_SIZE_T int"):
                # match pattern for intor
                intor = line.split("(")[0].split()[-1].strip()
                is_intor_valid = False
                for suffix in valid_suffix:
                    if intor.endswith(suffix):
                        intor = intor.replace(suffix, "")
                        is_intor_valid = True
                if "#" in intor or intor in actual or not is_intor_valid:
                    continue
                # match pattern for ng
                ng = None
                for l in lines[n:n+5]:
                    if l.strip().startswith("FINT ng[]") or l.strip().startswith("int ng[]"):
                        ng = [int(i) for i in l.strip().replace(";", "").replace("}", "").split("{")[-1].split(",")]
                if ng is None:
                    raise ValueError(f"{intor}")
                actual[intor] = ng

    spinor_only = set()
    for text in texts:
        for line in text.splitlines():
            if line.lstrip().startswith("#define"):
                continue
            m = re.search(r"\bBREIT0\(\s*(\w+)\s*,", line)
            if m is None:
                continue
            intor = f"int2e_breit_{m.group(1)}"
            if intor in actual:
                continue
            # The Breit driver calls the gaunt/gauge components; use the gaunt
            # component's ng metadata for the wrapper.
            base = f"int2e_{m.group(1)}"
            if base not in actual:
                raise ValueError(f"BREIT0 component {base} is not an integral")
            actual[intor] = actual[base]
            spinor_only.add(intor)

    return actual, spinor_only


def scan_ssc_defs(src_paths):
    """Scan `_ssc` entries, which are outside the `_sph/_cart/_spinor` taxonomy.

    These functions are defined in the libcint sources but deliberately omitted
    from `cint_funcs.h`, so the header/bindgen path cannot see them.  The `_ssc`
    suffix means `(representation, representation | Cartesian)`: centers 0/1 use
    the representation from the driver, center 2 is Cartesian.  Each family
    shares a `<family>_ssc_optimizer`.
    """
    ssc = {}
    optimizers = set()
    for src_path in src_paths:
        with open(src_path, "r") as f:
            text = f.read()
        for m in re.finditer(r"CACHE_SIZE_T\s+(int\w*?_ssc)\s*\(([^)]*)\)\s*\{", text):
            name = m.group(1)
            signature = m.group(2)
            body = text[m.end():m.end() + 900]
            ng_match = re.search(r"FINT ng\[\]\s*=\s*\{([^}]*)\}", body)
            if ng_match is None:
                raise ValueError(f"{name}: cannot find ng[]")
            drv_match = re.search(r"return\s+(CINT\w+_drv)\s*\(", body)
            is_spinor = "complex" in signature or (drv_match is not None and "spinor" in drv_match.group(1))
            ssc[name] = {
                "ng": [int(x) for x in ng_match.group(1).split(",")],
                "repr": "Spinor" if is_spinor else "Spheric",
                "family": name.split("_")[0],
            }
        for m in re.finditer(r"void\s+(int\w*?_ssc_optimizer)\s*\(", text):
            optimizers.add(m.group(1))
    return ssc, optimizers


actual_intor_qcint, _ = scan_intor_defs(qcint_src_paths)
# -

# ## Preparasion

# ### python import

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
# pattern:
#   CACHE_SIZE_T <intor>(...
#   FINT ng[] = {...};
# pattern rule:
#   <intor> ends with _sph, _cart, _spinor; no # in <intor>

actual_intor, spinor_only_intor = scan_intor_defs(
    glob.iglob(f"{dir_libcint_src}/src/**/*.c", recursive=True)
)
# -

# `_ssc` entries are not part of the standard header; scan them separately.
ssc_intor, ssc_optimizer = scan_ssc_defs(
    sorted(glob.glob(f"{dir_libcint_src}/src/**/*.c", recursive=True))
)
assert ssc_intor, "no `_ssc` integral definitions found in the libcint sources"
for _name, _meta in ssc_intor.items():
    _family = _meta["family"]
    assert f"{_family}_ssc_optimizer" in ssc_optimizer, f"{_name}: missing {_family}_ssc_optimizer"
    assert int(_family[3]) == 3, f"{_name}: expected a three-center `_ssc` integral"

# assert known intors are all included in all intors
assert len(set(known_intor) - set(actual_intor)) == 0


# ### Declare exceptions

def rule_qcint_exclusion(intor):
    return (intor in actual_intor) and (intor not in actual_intor_qcint)


def rule_with_f12(intor):
    return "_stg" in intor or "_yp" in intor


def rule_with_4c1e(intor):
    return "int4c1e" in intor


# integrals introduced in libcint v6.2.0; gated by cargo feature `api-v6_2`
api_v6_2_intor = {
    "int1e_ipr_origj",
    "int1e_r_origj_ip",
    "int1e_iprr_origj",
    "int1e_rr_origj_ip",
    "int3c1e_ipip1",
    "int3c1e_ip1ip2",
    "int3c1e_ipvip1",
    "int3c1e_ipip2",
}


def rule_api_v6_2(intor):
    return intor in api_v6_2_intor


# integrals added to qcint in v6.2.0 but present in libcint since v6.1.x;
# available for libcint any version, but for qcint only with `api-v6_2`.
# Based on actual symbols/source, not qcint's installed header: qcint 6.1.2
# declares the first four below in cint_funcs.h without exporting them, and
# qcint 6.2.0 exports int3c2e and ~20 other integrals not declared by its header.
qcint_v6_2_intor = {
    "int1e_ipipipiprinv",
    "int1e_ipipipnuc",
    "int1e_ipipiprinv",
    "int1e_ipipiprinvip",
    "int1e_ipipnucip",
    "int1e_ipiprinvip",
    "int1e_ipiprinvipip",
    "int1e_ipiprinvrip",
    "int1e_iprinviprip",
    "int1e_iprinvr",
    "int1e_r2_origi_ip2",
    "int1e_r4_origi_ip2",
    "int1e_rinvipiprip",
    "int3c1e_ip1_r2_origk",
    "int3c1e_ip1_r4_origk",
    "int3c1e_ip1_r6_origk",
}


def rule_qcint_v6_2(intor):
    return intor in qcint_v6_2_intor


# `int3c2e*_spinor` drivers evaluate two spinor centers against one scalar
# spherical auxiliary center.  The `_ssc` variants use a Cartesian auxiliary
# center and must not be added to this list.
int3c2e_spinor_sph_aux_layout = {
    "int3c2e",
    "int3c2e_ip1",
    "int3c2e_ip2",
    "int3c2e_pvp1",
    "int3c2e_pvxp1",
    "int3c2e_ig1",
    "int3c2e_spsp1",
    "int3c2e_ipspsp1",
    "int3c2e_spsp1ip2",
    "int3c2e_ipip1",
    "int3c2e_ipip2",
    "int3c2e_ipvip1",
    "int3c2e_ip1ip2",
}

def gen_ssc_ffi_block():
    out = [
        "",
        "",
        "/* `int3c2e*_ssc` functions are exported by libcint but omitted from its",
        " * `cint_funcs.h`, so they are declared here. qcint has the internal",
        " * `is_ssc`/`c2s_*_ssc` helpers but no public `_ssc` entry points. The `_c`",
        " * suffix avoids colliding with the unit-struct integrator names in",
        " * `cint_wrapper.rs`. */",
        '#[cfg(not(feature = "qcint"))]',
        'unsafe extern "C" {',
    ]
    for optimizer in sorted(ssc_optimizer):
        out += [
            f"    pub fn {optimizer}(",
            "        opt: *mut *mut CINTOpt,",
            "        atm: *const c_int,",
            "        natm: c_int,",
            "        bas: *const c_int,",
            "        nbas: c_int,",
            "        env: *const f64,",
            "    );",
        ]
    for name, meta in ssc_intor.items():
        out_ty = "__BindgenComplex<f64>" if meta["repr"] == "Spinor" else "f64"
        out += [
            f'    #[link_name = "{name}"]',
            f"    pub fn {name}_c(",
            f"        out: *mut {out_ty},",
            "        dims: *const c_int,",
            "        shls: *const c_int,",
            "        atm: *const c_int,",
            "        natm: c_int,",
            "        bas: *const c_int,",
            "        nbas: c_int,",
            "        env: *const f64,",
            "        opt: *const CINTOpt,",
            "        cache: *mut f64,",
            "    ) -> c_int;",
        ]
    out += ["}", ""]
    return "\n".join(out)


def gen_ssc_impl_entries():
    out = []
    for name, meta in ssc_intor.items():
        ng = meta["ng"]
        comp_all = max(ng[5], 1) * max(ng[6], 1) * ng[7]
        comp_tensor = ng[7]
        family = meta["family"]
        n_center = int(family[3])
        if meta["repr"] == "Spheric":
            sph, cart, spinor = f"{name}_c", "panic_cart", "panic_spinor"
            avail = ("true", "false", "false")
            driver = "Spheric"
            centers = "&[CIntType::Spheric, CIntType::Spheric, CIntType::Cartesian]"
        else:
            sph, cart, spinor = "panic_sph", "panic_cart", f"{name}_c"
            avail = ("false", "false", "true")
            driver = "Spinor"
            centers = "&[CIntType::Spinor, CIntType::Spinor, CIntType::Cartesian]"
        out.append(f'''#[cfg(not(feature = "qcint"))]
impl_integrator!(
    {name},
    {family}_ssc_optimizer,
    {sph},
    {cart},
    {spinor},
    {avail[0]},
    {avail[1]},
    {avail[2]},
    {comp_all},
    {comp_tensor},
    {n_center},
    vec!{ng},
    "{family}",
    "{name}",
    CIntKind::Int,
    CIntIntegralLayout {{ driver_type: CIntType::{driver}, center_types: {centers} }}
);
''')
    return "".join(out)


def gen_ssc_match_arms():
    return "".join(
        f'''#[cfg(not(feature = "qcint"))]
        "{name}" => Some(Box::new({name})),
'''
        for name in ssc_intor
    )


ssc_ffi_block = gen_ssc_ffi_block()

ssc_impl_entries = gen_ssc_impl_entries()

ssc_match_arms = gen_ssc_match_arms()



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
# exclusion (F12): _stg_, _yp_; spinor-only: BREIT0
token = """
extern CINTOptimizerFunction       {0:}_optimizer;
extern CINTIntegralFunctionReal    {0:}_cart;
extern CINTIntegralFunctionReal    {0:}_sph;
extern CINTIntegralFunctionComplex {0:}_spinor;
"""
cint_funcs = cint_funcs + "".join([
    token.format(intor)
    for intor in sorted(set(actual_intor))
    if not rule_with_f12(intor) and intor not in spinor_only_intor])
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

# BREIT0 spinor-only intors: only `_optimizer` and `_spinor` are defined
token = """
extern CINTOptimizerFunction       {0:}_optimizer;
extern CINTIntegralFunctionComplex {0:}_spinor;
"""
cint_funcs = cint_funcs + "".join([
    token.format(intor)
    for intor in sorted(spinor_only_intor)])
cint_funcs = cint_funcs.replace("#include <cint.h>", "#include \"cint.h\"")

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
    f.write(token + ssc_ffi_block)

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

use crate::cint::{CIntIntegralLayout, CIntKind, CIntType};
use crate::ffi::cint_ffi::*;
use crate::ffi::wrapper_traits::Integrator;
use crate::impl_integrator;
use core::any::Any;
use core::ffi::{c_int, c_void};

#[allow(unused_imports)]
use crate::ffi::wrapper_traits::{panic_cart, panic_sph, panic_spinor};

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
    if intor in spinor_only_intor:
        is_sph_available = "false"
        is_cart_available = "false"
        func_integral_sph = "panic_sph"
        func_integral_cart = "panic_cart"

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
    center_layout = ""
    if intor in int3c2e_spinor_sph_aux_layout:
        center_layout = """,
        CIntIntegralLayout {
            driver_type: CIntType::Spinor,
            center_types: &[CIntType::Spinor, CIntType::Spinor, CIntType::Spheric],
        }"""
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
        CIntKind::Int{center_layout});\n"""

    # cargo features for control
    if rule_with_f12(intor):
        token = '#[cfg(feature = "with_f12")]\n' + token
    if rule_with_4c1e(intor):
        token = '#[cfg(feature = "with_4c1e")]\n' + token
    if rule_qcint_v6_2(intor):
        token = '#[cfg(any(not(feature = "qcint"), feature = "api-v6_2"))]\n' + token
    elif rule_qcint_exclusion(intor):
        token = '#[cfg(not(feature = "qcint"))]\n' + token
    if rule_api_v6_2(intor):
        token = '#[cfg(feature = "api-v6_2")]\n' + token

    return token

for intor in actual_intor:
    token_wrapper += gen_impl_integrator(intor)
    if intor == "int3c2e":
        token_wrapper += ssc_impl_entries

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
    if rule_qcint_v6_2(intor):
        token_wrapper += '#[cfg(any(not(feature = "qcint"), feature = "api-v6_2"))]\n'
    elif rule_qcint_exclusion(intor):
        token_wrapper += '#[cfg(not(feature = "qcint"))]\n'
    if rule_api_v6_2(intor):
        token_wrapper += '#[cfg(feature = "api-v6_2")]\n'

    token_wrapper += f"""
        "{intor}" => Some(Box::new({intor})),
    """.strip()
    if intor == "int3c2e":
        token_wrapper += "\n" + ssc_match_arms

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
