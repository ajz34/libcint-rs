/* File from pyscf/pyscf: lib/gto/nr_ecp.h

   Copyright 2014-2018 The PySCF Developers. All Rights Reserved.
  
   Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at
 
        http://www.apache.org/licenses/LICENSE-2.0
 
    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.

 *
 * Author: Qiming Sun <osirpt.sun@gmail.com>
 */

#ifndef HAVE_DEFINED_ECP_H
#define HAVE_DEFINED_ECP_H
// #include <stdint.h>
#include <complex.h>

//#define ATOM_OF         0
//#define ANG_OF          1
#define RADI_POWER      3
#define SO_TYPE_OF      4

#define ECP_LMAX        5
//#define PTR_EXP         5
//#define PTR_COEFF       6
#define SIM_ZERO        1e-50
#define EXPCUTOFF       39   // 1e-17
#define CUTOFF          460  // ~ 1e200
#define CLOSE_ENOUGH(x, y)      (fabs(x-y) < 1e-12*fabs(y) || fabs(x-y) < 1e-12)
#define SQUARE(r)       (r[0]*r[0]+r[1]*r[1]+r[2]*r[2])
#define CART_CUM        (455+1) // upto l = 12
#define K_TAYLOR_MAX    7
#define K_TAB_COL       24      // >= (7*2+1+K_TAYLOR_MAX)
#define K_TAB_ENTRIES   400
#define K_TAB_INTERVAL  (16./K_TAB_ENTRIES)    // [0,16], interval 0.04
#define MALLOC_INSTACK(var, n) \
                var = (void *)cache; \
                cache = (void *)(((uintptr_t)(var + (n)) + 7) & (-(uintptr_t)8));
#define MARK_STACK      cache0 = cache;
#define RESTORE_STACK   cache = cache0;

// Held in env, to get *ecpbas, necpbas
#define AS_RINV_ORIG_ATOM       17
#define AS_ECPBAS_OFFSET        18
#define AS_NECPBAS              19

// for radial grids
#define LEVEL0          5
//#define LEVEL_MAX       11      // 2047 points
#define LEVEL_MAX       11
#define MIN(X, Y)       ((X) < (Y) ? (X) : (Y))
#define MAX(X, Y)       ((X) > (Y) ? (X) : (Y))
#endif

typedef struct {
    double *u_ecp;
} ECPOpt;

void NPdset0(double *p, const unsigned long n);

void ECPdel_optimizer(ECPOpt **opt);

int ECPscalar_sph(
    double *out, int *dims, int *shls, int *atm, int natm,
    int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPscalar_cart(
    double *out, int *dims, int *shls, int *atm, int natm,
    int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
void ECPscalar_optimizer(ECPOpt **opt, int *atm, int natm, int *bas, int nbas, double *env);

int ECPso_sph(double *out, int *dims, int *shls, int *atm, int natm,
              int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPso_cart(double *out, int *dims, int *shls, int *atm, int natm,
               int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPso_spinor(double complex *out, int *dims, int *shls, int *atm, int natm,
                 int *bas, int nbas, double *env, ECPOpt *opt, double *cache);

int ECPscalar_ignuc_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPscalar_ignuc_cart(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
void ECPscalar_ignuc_optimizer(ECPOpt **opt, int *atm, int natm, int *bas, int nbas, double *env);

int ECPscalar_ipnuc_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPscalar_ipnuc_cart(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
void ECPscalar_ipnuc_optimizer(ECPOpt **opt, int *atm, int natm, int *bas, int nbas, double *env);

int ECPscalar_ipipnuc_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPscalar_ipipnuc_cart(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
void ECPscalar_ipipnuc_optimizer(ECPOpt **opt, int *atm, int natm, int *bas, int nbas, double *env);

int ECPscalar_ipnucip_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPscalar_ipnucip_cart(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
void ECPscalar_ipnucip_optimizer(ECPOpt **opt, int *atm, int natm, int *bas, int nbas, double *env);

int ECPscalar_iprinv_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPscalar_iprinv_cart(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
void ECPscalar_iprinv_optimizer(ECPOpt **opt, int *atm, int natm, int *bas, int nbas, double *env);

int ECPscalar_ipiprinv_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPscalar_ipiprinv_cart(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
void ECPscalar_ipiprinv_optimizer(ECPOpt **opt, int *atm, int natm, int *bas, int nbas, double *env);

int ECPscalar_iprinvip_sph(double *out, int *dims, int *shls, int *atm, int natm,
                  int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
int ECPscalar_iprinvip_cart(double *out, int *dims, int *shls, int *atm, int natm,
                   int *bas, int nbas, double *env, ECPOpt *opt, double *cache);
void ECPscalar_iprinvip_optimizer(ECPOpt **opt, int *atm, int natm, int *bas, int nbas, double *env);


