
/*
 * Copyright (C) 2013-  Qiming Sun <osirpt.sun@gmail.com>
 *
 * Parameters and function signature for libcint.
 */

#define CINT_VERSION            "6.1.1"
#define CINT_SOVERSION          6

/* #undef I8 */
#ifdef I8
#include <stdint.h>
#define FINT int64_t
#else
#define FINT int
#endif

/* #undef CACHE_SIZE_I8 */
#ifdef CACHE_SIZE_I8
#include <stdint.h>
#define CACHE_SIZE_T int64_t
#else
#define CACHE_SIZE_T FINT
#endif

// global parameters in env
// Overall cutoff for integral prescreening, value needs to be ~ln(threshold)
#define PTR_EXPCUTOFF           0
// R_C of (r-R_C) in dipole, GIAO operators
#define PTR_COMMON_ORIG         1
// R_O in 1/|r-R_O|
#define PTR_RINV_ORIG           4
// ZETA parameter for Gaussian charge distribution (Gaussian nuclear model)
#define PTR_RINV_ZETA           7
// omega parameter in range-separated coulomb operator
// LR interaction: erf(omega*r12)/r12 if omega > 0
// SR interaction: erfc(omega*r12)/r12 if omega < 0
#define PTR_RANGE_OMEGA         8
// Yukawa potential and Slater-type geminal e^{-zeta r}
#define PTR_F12_ZETA            9
// Gaussian type geminal e^{-zeta r^2}
#define PTR_GTG_ZETA            10
#define NGRIDS                  11
#define PTR_GRIDS               12
#define PTR_ENV_START           20


// slots of atm
#define CHARGE_OF       0
#define PTR_COORD       1
#define NUC_MOD_OF      2
#define PTR_ZETA        3
#define PTR_FRAC_CHARGE 4
#define RESERVE_ATMSLOT 5
#define ATM_SLOTS       6


// slots of bas
#define ATOM_OF         0
#define ANG_OF          1
#define NPRIM_OF        2
#define NCTR_OF         3
#define KAPPA_OF        4
#define PTR_EXP         5
#define PTR_COEFF       6
#define RESERVE_BASLOT  7
#define BAS_SLOTS       8

// slots of gout
#define POSX            0
#define POSY            1
#define POSZ            2
#define POS1            3
// For 2-electron integral with two spin operators
// SIGMA1X * SIGMA2X     0
// SIGMA1Y * SIGMA2X     1
// SIGMA1Z * SIGMA2X     2
// I1_2x2  * SIGMA2X     3
// SIGMA1X * SIGMA2Y     4
// SIGMA1Y * SIGMA2Y     5
// SIGMA1Z * SIGMA2Y     6
// I1_2x2  * SIGMA2Y     7
// SIGMA1X * SIGMA2Z     8
// SIGMA1Y * SIGMA2Z     9
// SIGMA1Z * SIGMA2Z     10
// I1_2x2  * SIGMA2Z     11
// SIGMA1X * I2_2x2      12
// SIGMA1Y * I2_2x2      13
// SIGMA1Z * I2_2x2      14
// I1_2x2  * I2_2x2      15
#define POSXX           0
#define POSYX           1
#define POSZX           2
#define POS1X           3
#define POSXY           4
#define POSYY           5
#define POSZY           6
#define POS1Y           7
#define POSXZ           8
#define POSYZ           9
#define POSZZ           10
#define POS1Z           11
#define POSX1           12
#define POSY1           13
#define POSZ1           14
#define POS11           15

// tensor
#define TSRX        0
#define TSRY        1
#define TSRZ        2
#define TSRXX       0
#define TSRXY       1
#define TSRXZ       2
#define TSRYX       3
#define TSRYY       4
#define TSRYZ       5
#define TSRZX       6
#define TSRZY       7
#define TSRZZ       8

// other boundaries
#define MXRYSROOTS      32 // > ANG_MAX*2+1 for 4c2e
#define ANG_MAX         15 // l = 0..15
#define LMAX1           16 // > ANG_MAX
#define CART_MAX        136 // > (ANG_MAX*(ANG_MAX+1)/2)
#define SHLS_MAX        1048576
#define NPRIM_MAX       64
#define NCTR_MAX        64

#define POINT_NUC       1
#define GAUSSIAN_NUC    2
#define FRAC_CHARGE_NUC 3

#define bas(SLOT,I)     bas[BAS_SLOTS * (I) + (SLOT)]
#define atm(SLOT,I)     atm[ATM_SLOTS * (I) + (SLOT)]

#if !defined HAVE_DEFINED_CINTOPT_H
#define HAVE_DEFINED_CINTOPT_H
typedef struct {
    double rij[3];
    double eij;
    double cceij;
} PairData;
typedef struct {
    FINT **index_xyz_array; // LMAX1**4 pointers to index_xyz
    FINT **non0ctr;
    FINT **sortedidx;
    FINT nbas;
    double **log_max_coeff;
    PairData **pairdata;  // NULL indicates not-initialized, NO_VALUE can be skipped
} CINTOpt;

// Add this macro def to make pyscf compatible with both v4 and v5
#define HAVE_DEFINED_CINTENVVARS_H
typedef struct {
        FINT *atm;
        FINT *bas;
        double *env;
        FINT *shls;
        FINT natm;
        FINT nbas;

        FINT i_l;
        FINT j_l;
        FINT k_l;
        FINT l_l;
        FINT nfi;  // number of cartesian components
        FINT nfj;
        // in int1e_grids, the grids_offset and the number of grids
        union {FINT nfk; FINT grids_offset;};
        union {FINT nfl; FINT ngrids;};
        FINT nf;  // = nfi*nfj*nfk*nfl;
        FINT rys_order; // = nrys_roots for regular ERIs. can be nrys_roots/2 for SR ERIs
        FINT x_ctr[4];

        FINT gbits;
        FINT ncomp_e1; // = 1 if spin free, = 4 when spin included, it
        FINT ncomp_e2; // corresponds to POSX,POSY,POSZ,POS1, see cint.h
        FINT ncomp_tensor; // e.g. = 3 for gradients

        /* values may diff based on the g0_2d4d algorithm */
        FINT li_ceil; // power of x, == i_l if nabla is involved, otherwise == i_l
        FINT lj_ceil;
        FINT lk_ceil;
        FINT ll_ceil;
        FINT g_stride_i; // nrys_roots * shift of (i++,k,l,j)
        FINT g_stride_k; // nrys_roots * shift of (i,k++,l,j)
        FINT g_stride_l; // nrys_roots * shift of (i,k,l++,j)
        FINT g_stride_j; // nrys_roots * shift of (i,k,l,j++)
        FINT nrys_roots;
        FINT g_size;  // ref to cint2e.c g = malloc(sizeof(double)*g_size)

        FINT g2d_ijmax;
        FINT g2d_klmax;
        double common_factor;
        double expcutoff;
        double rirj[3]; // diff by sign in different g0_2d4d algorithm
        double rkrl[3];
        double *rx_in_rijrx;
        double *rx_in_rklrx;

        double *ri;
        double *rj;
        double *rk;
        // in int2e or int3c2e, the coordinates of the fourth shell
        // in int1e_grids, the pointer for the grids coordinates
        union {double *rl; double *grids;};

        FINT (*f_g0_2e)();
        void (*f_g0_2d4d)();
        void (*f_gout)();
        CINTOpt *opt;

        /* values are assigned during calculation */
        int *idx;
        double ai[1];
        double aj[1];
        double ak[1];
        double al[1];
        double fac[1];
        double rij[3];
        double rkl[3];
} CINTEnvVars;
#endif

FINT CINTlen_cart(const FINT l);
FINT CINTlen_spinor(const FINT bas_id, const FINT *bas);

FINT CINTcgtos_cart(const FINT bas_id, const FINT *bas);
FINT CINTcgtos_spheric(const FINT bas_id, const FINT *bas);
FINT CINTcgtos_spinor(const FINT bas_id, const FINT *bas);
FINT CINTcgto_cart(const FINT bas_id, const FINT *bas);
FINT CINTcgto_spheric(const FINT bas_id, const FINT *bas);
FINT CINTcgto_spinor(const FINT bas_id, const FINT *bas);

FINT CINTtot_pgto_spheric(const FINT *bas, const FINT nbas);
FINT CINTtot_pgto_spinor(const FINT *bas, const FINT nbas);

FINT CINTtot_cgto_cart(const FINT *bas, const FINT nbas);
FINT CINTtot_cgto_spheric(const FINT *bas, const FINT nbas);
FINT CINTtot_cgto_spinor(const FINT *bas, const FINT nbas);

void CINTshells_cart_offset(FINT ao_loc[], const FINT *bas, const FINT nbas);
void CINTshells_spheric_offset(FINT ao_loc[], const FINT *bas, const FINT nbas);
void CINTshells_spinor_offset(FINT ao_loc[], const FINT *bas, const FINT nbas);

double *CINTc2s_bra_sph(double *sph, FINT nket, double *cart, FINT l);
double *CINTc2s_ket_sph(double *sph, FINT nket, double *cart, FINT l);
double *CINTc2s_ket_sph1(double *sph, double *cart, FINT lds, FINT ldc, FINT l);


double CINTgto_norm(FINT n, double a);


void CINTinit_2e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env);
void CINTinit_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                        FINT *bas, FINT nbas, double *env);
void CINTdel_2e_optimizer(CINTOpt **opt);
void CINTdel_optimizer(CINTOpt **opt);


FINT cint2e_cart(double *opijkl, FINT *shls,
                FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
                CINTOpt *opt);
void cint2e_cart_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                           FINT *bas, FINT nbas, double *env);
FINT cint2e_sph(double *opijkl, FINT *shls,
               FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
               CINTOpt *opt);
void cint2e_sph_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                          FINT *bas, FINT nbas, double *env);
FINT cint2e(double *opijkl, FINT *shls,
           FINT *atm, FINT natm, FINT *bas, FINT nbas, double *env,
           CINTOpt *opt);
void cint2e_optimizer(CINTOpt **opt, FINT *atm, FINT natm,
                      FINT *bas, FINT nbas, double *env);

#ifndef __cplusplus
#include <complex.h>

void CINTc2s_ket_spinor_sf1(double complex *gspa, double complex *gspb, double *gcart,
                            FINT lds, FINT ldc, FINT nctr, FINT l, FINT kappa);
void CINTc2s_iket_spinor_sf1(double complex *gspa, double complex *gspb, double *gcart,
                             FINT lds, FINT ldc, FINT nctr, FINT l, FINT kappa);
void CINTc2s_ket_spinor_si1(double complex *gspa, double complex *gspb, double *gcart,
                            FINT lds, FINT ldc, FINT nctr, FINT l, FINT kappa);
void CINTc2s_iket_spinor_si1(double complex *gspa, double complex *gspb, double *gcart,
                             FINT lds, FINT ldc, FINT nctr, FINT l, FINT kappa);
#endif


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

extern CINTOptimizerFunction       int1e_a01gp_optimizer;
extern CINTIntegralFunctionReal    int1e_a01gp_cart;
extern CINTIntegralFunctionReal    int1e_a01gp_sph;
extern CINTIntegralFunctionComplex int1e_a01gp_spinor;

extern CINTOptimizerFunction       int1e_cg_a11part_optimizer;
extern CINTIntegralFunctionReal    int1e_cg_a11part_cart;
extern CINTIntegralFunctionReal    int1e_cg_a11part_sph;
extern CINTIntegralFunctionComplex int1e_cg_a11part_spinor;

extern CINTOptimizerFunction       int1e_cg_irxp_optimizer;
extern CINTIntegralFunctionReal    int1e_cg_irxp_cart;
extern CINTIntegralFunctionReal    int1e_cg_irxp_sph;
extern CINTIntegralFunctionComplex int1e_cg_irxp_spinor;

extern CINTOptimizerFunction       int1e_cg_sa10nucsp_optimizer;
extern CINTIntegralFunctionReal    int1e_cg_sa10nucsp_cart;
extern CINTIntegralFunctionReal    int1e_cg_sa10nucsp_sph;
extern CINTIntegralFunctionComplex int1e_cg_sa10nucsp_spinor;

extern CINTOptimizerFunction       int1e_cg_sa10sa01_optimizer;
extern CINTIntegralFunctionReal    int1e_cg_sa10sa01_cart;
extern CINTIntegralFunctionReal    int1e_cg_sa10sa01_sph;
extern CINTIntegralFunctionComplex int1e_cg_sa10sa01_spinor;

extern CINTOptimizerFunction       int1e_cg_sa10sp_optimizer;
extern CINTIntegralFunctionReal    int1e_cg_sa10sp_cart;
extern CINTIntegralFunctionReal    int1e_cg_sa10sp_sph;
extern CINTIntegralFunctionComplex int1e_cg_sa10sp_spinor;

extern CINTOptimizerFunction       int1e_drinv_optimizer;
extern CINTIntegralFunctionReal    int1e_drinv_cart;
extern CINTIntegralFunctionReal    int1e_drinv_sph;
extern CINTIntegralFunctionComplex int1e_drinv_spinor;

extern CINTOptimizerFunction       int1e_ggkin_optimizer;
extern CINTIntegralFunctionReal    int1e_ggkin_cart;
extern CINTIntegralFunctionReal    int1e_ggkin_sph;
extern CINTIntegralFunctionComplex int1e_ggkin_spinor;

extern CINTOptimizerFunction       int1e_ggnuc_optimizer;
extern CINTIntegralFunctionReal    int1e_ggnuc_cart;
extern CINTIntegralFunctionReal    int1e_ggnuc_sph;
extern CINTIntegralFunctionComplex int1e_ggnuc_spinor;

extern CINTOptimizerFunction       int1e_ggovlp_optimizer;
extern CINTIntegralFunctionReal    int1e_ggovlp_cart;
extern CINTIntegralFunctionReal    int1e_ggovlp_sph;
extern CINTIntegralFunctionComplex int1e_ggovlp_spinor;

extern CINTOptimizerFunction       int1e_giao_a11part_optimizer;
extern CINTIntegralFunctionReal    int1e_giao_a11part_cart;
extern CINTIntegralFunctionReal    int1e_giao_a11part_sph;
extern CINTIntegralFunctionComplex int1e_giao_a11part_spinor;

extern CINTOptimizerFunction       int1e_giao_irjxp_optimizer;
extern CINTIntegralFunctionReal    int1e_giao_irjxp_cart;
extern CINTIntegralFunctionReal    int1e_giao_irjxp_sph;
extern CINTIntegralFunctionComplex int1e_giao_irjxp_spinor;

extern CINTOptimizerFunction       int1e_giao_sa10nucsp_optimizer;
extern CINTIntegralFunctionReal    int1e_giao_sa10nucsp_cart;
extern CINTIntegralFunctionReal    int1e_giao_sa10nucsp_sph;
extern CINTIntegralFunctionComplex int1e_giao_sa10nucsp_spinor;

extern CINTOptimizerFunction       int1e_giao_sa10sa01_optimizer;
extern CINTIntegralFunctionReal    int1e_giao_sa10sa01_cart;
extern CINTIntegralFunctionReal    int1e_giao_sa10sa01_sph;
extern CINTIntegralFunctionComplex int1e_giao_sa10sa01_spinor;

extern CINTOptimizerFunction       int1e_giao_sa10sp_optimizer;
extern CINTIntegralFunctionReal    int1e_giao_sa10sp_cart;
extern CINTIntegralFunctionReal    int1e_giao_sa10sp_sph;
extern CINTIntegralFunctionComplex int1e_giao_sa10sp_spinor;

extern CINTOptimizerFunction       int1e_gnuc_optimizer;
extern CINTIntegralFunctionReal    int1e_gnuc_cart;
extern CINTIntegralFunctionReal    int1e_gnuc_sph;
extern CINTIntegralFunctionComplex int1e_gnuc_spinor;

extern CINTOptimizerFunction       int1e_govlp_optimizer;
extern CINTIntegralFunctionReal    int1e_govlp_cart;
extern CINTIntegralFunctionReal    int1e_govlp_sph;
extern CINTIntegralFunctionComplex int1e_govlp_spinor;

extern CINTOptimizerFunction       int1e_grids_optimizer;
extern CINTIntegralFunctionReal    int1e_grids_cart;
extern CINTIntegralFunctionReal    int1e_grids_sph;
extern CINTIntegralFunctionComplex int1e_grids_spinor;

extern CINTOptimizerFunction       int1e_grids_ip_optimizer;
extern CINTIntegralFunctionReal    int1e_grids_ip_cart;
extern CINTIntegralFunctionReal    int1e_grids_ip_sph;
extern CINTIntegralFunctionComplex int1e_grids_ip_spinor;

extern CINTOptimizerFunction       int1e_grids_ipip_optimizer;
extern CINTIntegralFunctionReal    int1e_grids_ipip_cart;
extern CINTIntegralFunctionReal    int1e_grids_ipip_sph;
extern CINTIntegralFunctionComplex int1e_grids_ipip_spinor;

extern CINTOptimizerFunction       int1e_grids_ipvip_optimizer;
extern CINTIntegralFunctionReal    int1e_grids_ipvip_cart;
extern CINTIntegralFunctionReal    int1e_grids_ipvip_sph;
extern CINTIntegralFunctionComplex int1e_grids_ipvip_spinor;

extern CINTOptimizerFunction       int1e_grids_spvsp_optimizer;
extern CINTIntegralFunctionReal    int1e_grids_spvsp_cart;
extern CINTIntegralFunctionReal    int1e_grids_spvsp_sph;
extern CINTIntegralFunctionComplex int1e_grids_spvsp_spinor;

extern CINTOptimizerFunction       int1e_grjxp_optimizer;
extern CINTIntegralFunctionReal    int1e_grjxp_cart;
extern CINTIntegralFunctionReal    int1e_grjxp_sph;
extern CINTIntegralFunctionComplex int1e_grjxp_spinor;

extern CINTOptimizerFunction       int1e_ia01p_optimizer;
extern CINTIntegralFunctionReal    int1e_ia01p_cart;
extern CINTIntegralFunctionReal    int1e_ia01p_sph;
extern CINTIntegralFunctionComplex int1e_ia01p_spinor;

extern CINTOptimizerFunction       int1e_igkin_optimizer;
extern CINTIntegralFunctionReal    int1e_igkin_cart;
extern CINTIntegralFunctionReal    int1e_igkin_sph;
extern CINTIntegralFunctionComplex int1e_igkin_spinor;

extern CINTOptimizerFunction       int1e_ignuc_optimizer;
extern CINTIntegralFunctionReal    int1e_ignuc_cart;
extern CINTIntegralFunctionReal    int1e_ignuc_sph;
extern CINTIntegralFunctionComplex int1e_ignuc_spinor;

extern CINTOptimizerFunction       int1e_igovlp_optimizer;
extern CINTIntegralFunctionReal    int1e_igovlp_cart;
extern CINTIntegralFunctionReal    int1e_igovlp_sph;
extern CINTIntegralFunctionComplex int1e_igovlp_spinor;

extern CINTOptimizerFunction       int1e_inuc_rcxp_optimizer;
extern CINTIntegralFunctionReal    int1e_inuc_rcxp_cart;
extern CINTIntegralFunctionReal    int1e_inuc_rcxp_sph;
extern CINTIntegralFunctionComplex int1e_inuc_rcxp_spinor;

extern CINTOptimizerFunction       int1e_inuc_rxp_optimizer;
extern CINTIntegralFunctionReal    int1e_inuc_rxp_cart;
extern CINTIntegralFunctionReal    int1e_inuc_rxp_sph;
extern CINTIntegralFunctionComplex int1e_inuc_rxp_spinor;

extern CINTOptimizerFunction       int1e_ipipipiprinv_optimizer;
extern CINTIntegralFunctionReal    int1e_ipipipiprinv_cart;
extern CINTIntegralFunctionReal    int1e_ipipipiprinv_sph;
extern CINTIntegralFunctionComplex int1e_ipipipiprinv_spinor;

extern CINTOptimizerFunction       int1e_ipipipnuc_optimizer;
extern CINTIntegralFunctionReal    int1e_ipipipnuc_cart;
extern CINTIntegralFunctionReal    int1e_ipipipnuc_sph;
extern CINTIntegralFunctionComplex int1e_ipipipnuc_spinor;

extern CINTOptimizerFunction       int1e_ipipiprinv_optimizer;
extern CINTIntegralFunctionReal    int1e_ipipiprinv_cart;
extern CINTIntegralFunctionReal    int1e_ipipiprinv_sph;
extern CINTIntegralFunctionComplex int1e_ipipiprinv_spinor;

extern CINTOptimizerFunction       int1e_ipipiprinvip_optimizer;
extern CINTIntegralFunctionReal    int1e_ipipiprinvip_cart;
extern CINTIntegralFunctionReal    int1e_ipipiprinvip_sph;
extern CINTIntegralFunctionComplex int1e_ipipiprinvip_spinor;

extern CINTOptimizerFunction       int1e_ipipkin_optimizer;
extern CINTIntegralFunctionReal    int1e_ipipkin_cart;
extern CINTIntegralFunctionReal    int1e_ipipkin_sph;
extern CINTIntegralFunctionComplex int1e_ipipkin_spinor;

extern CINTOptimizerFunction       int1e_ipipnuc_optimizer;
extern CINTIntegralFunctionReal    int1e_ipipnuc_cart;
extern CINTIntegralFunctionReal    int1e_ipipnuc_sph;
extern CINTIntegralFunctionComplex int1e_ipipnuc_spinor;

extern CINTOptimizerFunction       int1e_ipipnucip_optimizer;
extern CINTIntegralFunctionReal    int1e_ipipnucip_cart;
extern CINTIntegralFunctionReal    int1e_ipipnucip_sph;
extern CINTIntegralFunctionComplex int1e_ipipnucip_spinor;

extern CINTOptimizerFunction       int1e_ipipovlp_optimizer;
extern CINTIntegralFunctionReal    int1e_ipipovlp_cart;
extern CINTIntegralFunctionReal    int1e_ipipovlp_sph;
extern CINTIntegralFunctionComplex int1e_ipipovlp_spinor;

extern CINTOptimizerFunction       int1e_ipippnucp_optimizer;
extern CINTIntegralFunctionReal    int1e_ipippnucp_cart;
extern CINTIntegralFunctionReal    int1e_ipippnucp_sph;
extern CINTIntegralFunctionComplex int1e_ipippnucp_spinor;

extern CINTOptimizerFunction       int1e_ipipprinvp_optimizer;
extern CINTIntegralFunctionReal    int1e_ipipprinvp_cart;
extern CINTIntegralFunctionReal    int1e_ipipprinvp_sph;
extern CINTIntegralFunctionComplex int1e_ipipprinvp_spinor;

extern CINTOptimizerFunction       int1e_ipipr_optimizer;
extern CINTIntegralFunctionReal    int1e_ipipr_cart;
extern CINTIntegralFunctionReal    int1e_ipipr_sph;
extern CINTIntegralFunctionComplex int1e_ipipr_spinor;

extern CINTOptimizerFunction       int1e_ipiprinv_optimizer;
extern CINTIntegralFunctionReal    int1e_ipiprinv_cart;
extern CINTIntegralFunctionReal    int1e_ipiprinv_sph;
extern CINTIntegralFunctionComplex int1e_ipiprinv_spinor;

extern CINTOptimizerFunction       int1e_ipiprinvip_optimizer;
extern CINTIntegralFunctionReal    int1e_ipiprinvip_cart;
extern CINTIntegralFunctionReal    int1e_ipiprinvip_sph;
extern CINTIntegralFunctionComplex int1e_ipiprinvip_spinor;

extern CINTOptimizerFunction       int1e_ipiprinvipip_optimizer;
extern CINTIntegralFunctionReal    int1e_ipiprinvipip_cart;
extern CINTIntegralFunctionReal    int1e_ipiprinvipip_sph;
extern CINTIntegralFunctionComplex int1e_ipiprinvipip_spinor;

extern CINTOptimizerFunction       int1e_ipiprinvrip_optimizer;
extern CINTIntegralFunctionReal    int1e_ipiprinvrip_cart;
extern CINTIntegralFunctionReal    int1e_ipiprinvrip_sph;
extern CINTIntegralFunctionComplex int1e_ipiprinvrip_spinor;

extern CINTOptimizerFunction       int1e_ipipspnucsp_optimizer;
extern CINTIntegralFunctionReal    int1e_ipipspnucsp_cart;
extern CINTIntegralFunctionReal    int1e_ipipspnucsp_sph;
extern CINTIntegralFunctionComplex int1e_ipipspnucsp_spinor;

extern CINTOptimizerFunction       int1e_ipipsprinvsp_optimizer;
extern CINTIntegralFunctionReal    int1e_ipipsprinvsp_cart;
extern CINTIntegralFunctionReal    int1e_ipipsprinvsp_sph;
extern CINTIntegralFunctionComplex int1e_ipipsprinvsp_spinor;

extern CINTOptimizerFunction       int1e_ipkin_optimizer;
extern CINTIntegralFunctionReal    int1e_ipkin_cart;
extern CINTIntegralFunctionReal    int1e_ipkin_sph;
extern CINTIntegralFunctionComplex int1e_ipkin_spinor;

extern CINTOptimizerFunction       int1e_ipkinip_optimizer;
extern CINTIntegralFunctionReal    int1e_ipkinip_cart;
extern CINTIntegralFunctionReal    int1e_ipkinip_sph;
extern CINTIntegralFunctionComplex int1e_ipkinip_spinor;

extern CINTOptimizerFunction       int1e_ipnuc_optimizer;
extern CINTIntegralFunctionReal    int1e_ipnuc_cart;
extern CINTIntegralFunctionReal    int1e_ipnuc_sph;
extern CINTIntegralFunctionComplex int1e_ipnuc_spinor;

extern CINTOptimizerFunction       int1e_ipnucip_optimizer;
extern CINTIntegralFunctionReal    int1e_ipnucip_cart;
extern CINTIntegralFunctionReal    int1e_ipnucip_sph;
extern CINTIntegralFunctionComplex int1e_ipnucip_spinor;

extern CINTOptimizerFunction       int1e_ipovlp_optimizer;
extern CINTIntegralFunctionReal    int1e_ipovlp_cart;
extern CINTIntegralFunctionReal    int1e_ipovlp_sph;
extern CINTIntegralFunctionComplex int1e_ipovlp_spinor;

extern CINTOptimizerFunction       int1e_ipovlpip_optimizer;
extern CINTIntegralFunctionReal    int1e_ipovlpip_cart;
extern CINTIntegralFunctionReal    int1e_ipovlpip_sph;
extern CINTIntegralFunctionComplex int1e_ipovlpip_spinor;

extern CINTOptimizerFunction       int1e_ippnucp_optimizer;
extern CINTIntegralFunctionReal    int1e_ippnucp_cart;
extern CINTIntegralFunctionReal    int1e_ippnucp_sph;
extern CINTIntegralFunctionComplex int1e_ippnucp_spinor;

extern CINTOptimizerFunction       int1e_ippnucpip_optimizer;
extern CINTIntegralFunctionReal    int1e_ippnucpip_cart;
extern CINTIntegralFunctionReal    int1e_ippnucpip_sph;
extern CINTIntegralFunctionComplex int1e_ippnucpip_spinor;

extern CINTOptimizerFunction       int1e_ipprinvp_optimizer;
extern CINTIntegralFunctionReal    int1e_ipprinvp_cart;
extern CINTIntegralFunctionReal    int1e_ipprinvp_sph;
extern CINTIntegralFunctionComplex int1e_ipprinvp_spinor;

extern CINTOptimizerFunction       int1e_ipprinvpip_optimizer;
extern CINTIntegralFunctionReal    int1e_ipprinvpip_cart;
extern CINTIntegralFunctionReal    int1e_ipprinvpip_sph;
extern CINTIntegralFunctionComplex int1e_ipprinvpip_spinor;

extern CINTOptimizerFunction       int1e_iprinv_optimizer;
extern CINTIntegralFunctionReal    int1e_iprinv_cart;
extern CINTIntegralFunctionReal    int1e_iprinv_sph;
extern CINTIntegralFunctionComplex int1e_iprinv_spinor;

extern CINTOptimizerFunction       int1e_iprinvip_optimizer;
extern CINTIntegralFunctionReal    int1e_iprinvip_cart;
extern CINTIntegralFunctionReal    int1e_iprinvip_sph;
extern CINTIntegralFunctionComplex int1e_iprinvip_spinor;

extern CINTOptimizerFunction       int1e_iprinviprip_optimizer;
extern CINTIntegralFunctionReal    int1e_iprinviprip_cart;
extern CINTIntegralFunctionReal    int1e_iprinviprip_sph;
extern CINTIntegralFunctionComplex int1e_iprinviprip_spinor;

extern CINTOptimizerFunction       int1e_iprinvr_optimizer;
extern CINTIntegralFunctionReal    int1e_iprinvr_cart;
extern CINTIntegralFunctionReal    int1e_iprinvr_sph;
extern CINTIntegralFunctionComplex int1e_iprinvr_spinor;

extern CINTOptimizerFunction       int1e_iprip_optimizer;
extern CINTIntegralFunctionReal    int1e_iprip_cart;
extern CINTIntegralFunctionReal    int1e_iprip_sph;
extern CINTIntegralFunctionComplex int1e_iprip_spinor;

extern CINTOptimizerFunction       int1e_ipspnucsp_optimizer;
extern CINTIntegralFunctionReal    int1e_ipspnucsp_cart;
extern CINTIntegralFunctionReal    int1e_ipspnucsp_sph;
extern CINTIntegralFunctionComplex int1e_ipspnucsp_spinor;

extern CINTOptimizerFunction       int1e_ipspnucspip_optimizer;
extern CINTIntegralFunctionReal    int1e_ipspnucspip_cart;
extern CINTIntegralFunctionReal    int1e_ipspnucspip_sph;
extern CINTIntegralFunctionComplex int1e_ipspnucspip_spinor;

extern CINTOptimizerFunction       int1e_ipsprinvsp_optimizer;
extern CINTIntegralFunctionReal    int1e_ipsprinvsp_cart;
extern CINTIntegralFunctionReal    int1e_ipsprinvsp_sph;
extern CINTIntegralFunctionComplex int1e_ipsprinvsp_spinor;

extern CINTOptimizerFunction       int1e_ipsprinvspip_optimizer;
extern CINTIntegralFunctionReal    int1e_ipsprinvspip_cart;
extern CINTIntegralFunctionReal    int1e_ipsprinvspip_sph;
extern CINTIntegralFunctionComplex int1e_ipsprinvspip_spinor;

extern CINTOptimizerFunction       int1e_irp_optimizer;
extern CINTIntegralFunctionReal    int1e_irp_cart;
extern CINTIntegralFunctionReal    int1e_irp_sph;
extern CINTIntegralFunctionComplex int1e_irp_spinor;

extern CINTOptimizerFunction       int1e_irpr_optimizer;
extern CINTIntegralFunctionReal    int1e_irpr_cart;
extern CINTIntegralFunctionReal    int1e_irpr_sph;
extern CINTIntegralFunctionComplex int1e_irpr_spinor;

extern CINTOptimizerFunction       int1e_irrp_optimizer;
extern CINTIntegralFunctionReal    int1e_irrp_cart;
extern CINTIntegralFunctionReal    int1e_irrp_sph;
extern CINTIntegralFunctionComplex int1e_irrp_spinor;

extern CINTOptimizerFunction       int1e_kin_optimizer;
extern CINTIntegralFunctionReal    int1e_kin_cart;
extern CINTIntegralFunctionReal    int1e_kin_sph;
extern CINTIntegralFunctionComplex int1e_kin_spinor;

extern CINTOptimizerFunction       int1e_kinip_optimizer;
extern CINTIntegralFunctionReal    int1e_kinip_cart;
extern CINTIntegralFunctionReal    int1e_kinip_sph;
extern CINTIntegralFunctionComplex int1e_kinip_spinor;

extern CINTOptimizerFunction       int1e_nuc_optimizer;
extern CINTIntegralFunctionReal    int1e_nuc_cart;
extern CINTIntegralFunctionReal    int1e_nuc_sph;
extern CINTIntegralFunctionComplex int1e_nuc_spinor;

extern CINTOptimizerFunction       int1e_ovlp_optimizer;
extern CINTIntegralFunctionReal    int1e_ovlp_cart;
extern CINTIntegralFunctionReal    int1e_ovlp_sph;
extern CINTIntegralFunctionComplex int1e_ovlp_spinor;

extern CINTOptimizerFunction       int1e_ovlpip_optimizer;
extern CINTIntegralFunctionReal    int1e_ovlpip_cart;
extern CINTIntegralFunctionReal    int1e_ovlpip_sph;
extern CINTIntegralFunctionComplex int1e_ovlpip_spinor;

extern CINTOptimizerFunction       int1e_p4_optimizer;
extern CINTIntegralFunctionReal    int1e_p4_cart;
extern CINTIntegralFunctionReal    int1e_p4_sph;
extern CINTIntegralFunctionComplex int1e_p4_spinor;

extern CINTOptimizerFunction       int1e_pnucp_optimizer;
extern CINTIntegralFunctionReal    int1e_pnucp_cart;
extern CINTIntegralFunctionReal    int1e_pnucp_sph;
extern CINTIntegralFunctionComplex int1e_pnucp_spinor;

extern CINTOptimizerFunction       int1e_pnucxp_optimizer;
extern CINTIntegralFunctionReal    int1e_pnucxp_cart;
extern CINTIntegralFunctionReal    int1e_pnucxp_sph;
extern CINTIntegralFunctionComplex int1e_pnucxp_spinor;

extern CINTOptimizerFunction       int1e_prinvp_optimizer;
extern CINTIntegralFunctionReal    int1e_prinvp_cart;
extern CINTIntegralFunctionReal    int1e_prinvp_sph;
extern CINTIntegralFunctionComplex int1e_prinvp_spinor;

extern CINTOptimizerFunction       int1e_prinvxp_optimizer;
extern CINTIntegralFunctionReal    int1e_prinvxp_cart;
extern CINTIntegralFunctionReal    int1e_prinvxp_sph;
extern CINTIntegralFunctionComplex int1e_prinvxp_spinor;

extern CINTOptimizerFunction       int1e_r_optimizer;
extern CINTIntegralFunctionReal    int1e_r_cart;
extern CINTIntegralFunctionReal    int1e_r_sph;
extern CINTIntegralFunctionComplex int1e_r_spinor;

extern CINTOptimizerFunction       int1e_r2_optimizer;
extern CINTIntegralFunctionReal    int1e_r2_cart;
extern CINTIntegralFunctionReal    int1e_r2_sph;
extern CINTIntegralFunctionComplex int1e_r2_spinor;

extern CINTOptimizerFunction       int1e_r2_origi_optimizer;
extern CINTIntegralFunctionReal    int1e_r2_origi_cart;
extern CINTIntegralFunctionReal    int1e_r2_origi_sph;
extern CINTIntegralFunctionComplex int1e_r2_origi_spinor;

extern CINTOptimizerFunction       int1e_r2_origi_ip2_optimizer;
extern CINTIntegralFunctionReal    int1e_r2_origi_ip2_cart;
extern CINTIntegralFunctionReal    int1e_r2_origi_ip2_sph;
extern CINTIntegralFunctionComplex int1e_r2_origi_ip2_spinor;

extern CINTOptimizerFunction       int1e_r2_origj_optimizer;
extern CINTIntegralFunctionReal    int1e_r2_origj_cart;
extern CINTIntegralFunctionReal    int1e_r2_origj_sph;
extern CINTIntegralFunctionComplex int1e_r2_origj_spinor;

extern CINTOptimizerFunction       int1e_r4_optimizer;
extern CINTIntegralFunctionReal    int1e_r4_cart;
extern CINTIntegralFunctionReal    int1e_r4_sph;
extern CINTIntegralFunctionComplex int1e_r4_spinor;

extern CINTOptimizerFunction       int1e_r4_origi_optimizer;
extern CINTIntegralFunctionReal    int1e_r4_origi_cart;
extern CINTIntegralFunctionReal    int1e_r4_origi_sph;
extern CINTIntegralFunctionComplex int1e_r4_origi_spinor;

extern CINTOptimizerFunction       int1e_r4_origi_ip2_optimizer;
extern CINTIntegralFunctionReal    int1e_r4_origi_ip2_cart;
extern CINTIntegralFunctionReal    int1e_r4_origi_ip2_sph;
extern CINTIntegralFunctionComplex int1e_r4_origi_ip2_spinor;

extern CINTOptimizerFunction       int1e_r4_origj_optimizer;
extern CINTIntegralFunctionReal    int1e_r4_origj_cart;
extern CINTIntegralFunctionReal    int1e_r4_origj_sph;
extern CINTIntegralFunctionComplex int1e_r4_origj_spinor;

extern CINTOptimizerFunction       int1e_r_origj_optimizer;
extern CINTIntegralFunctionReal    int1e_r_origj_cart;
extern CINTIntegralFunctionReal    int1e_r_origj_sph;
extern CINTIntegralFunctionComplex int1e_r_origj_spinor;

extern CINTOptimizerFunction       int1e_rinv_optimizer;
extern CINTIntegralFunctionReal    int1e_rinv_cart;
extern CINTIntegralFunctionReal    int1e_rinv_sph;
extern CINTIntegralFunctionComplex int1e_rinv_spinor;

extern CINTOptimizerFunction       int1e_rinvipiprip_optimizer;
extern CINTIntegralFunctionReal    int1e_rinvipiprip_cart;
extern CINTIntegralFunctionReal    int1e_rinvipiprip_sph;
extern CINTIntegralFunctionComplex int1e_rinvipiprip_spinor;

extern CINTOptimizerFunction       int1e_rr_optimizer;
extern CINTIntegralFunctionReal    int1e_rr_cart;
extern CINTIntegralFunctionReal    int1e_rr_sph;
extern CINTIntegralFunctionComplex int1e_rr_spinor;

extern CINTOptimizerFunction       int1e_rr_origj_optimizer;
extern CINTIntegralFunctionReal    int1e_rr_origj_cart;
extern CINTIntegralFunctionReal    int1e_rr_origj_sph;
extern CINTIntegralFunctionComplex int1e_rr_origj_spinor;

extern CINTOptimizerFunction       int1e_rrr_optimizer;
extern CINTIntegralFunctionReal    int1e_rrr_cart;
extern CINTIntegralFunctionReal    int1e_rrr_sph;
extern CINTIntegralFunctionComplex int1e_rrr_spinor;

extern CINTOptimizerFunction       int1e_rrrr_optimizer;
extern CINTIntegralFunctionReal    int1e_rrrr_cart;
extern CINTIntegralFunctionReal    int1e_rrrr_sph;
extern CINTIntegralFunctionComplex int1e_rrrr_spinor;

extern CINTOptimizerFunction       int1e_sa01sp_optimizer;
extern CINTIntegralFunctionReal    int1e_sa01sp_cart;
extern CINTIntegralFunctionReal    int1e_sa01sp_sph;
extern CINTIntegralFunctionComplex int1e_sa01sp_spinor;

extern CINTOptimizerFunction       int1e_sigma_optimizer;
extern CINTIntegralFunctionReal    int1e_sigma_cart;
extern CINTIntegralFunctionReal    int1e_sigma_sph;
extern CINTIntegralFunctionComplex int1e_sigma_spinor;

extern CINTOptimizerFunction       int1e_sp_optimizer;
extern CINTIntegralFunctionReal    int1e_sp_cart;
extern CINTIntegralFunctionReal    int1e_sp_sph;
extern CINTIntegralFunctionComplex int1e_sp_spinor;

extern CINTOptimizerFunction       int1e_spgnucsp_optimizer;
extern CINTIntegralFunctionReal    int1e_spgnucsp_cart;
extern CINTIntegralFunctionReal    int1e_spgnucsp_sph;
extern CINTIntegralFunctionComplex int1e_spgnucsp_spinor;

extern CINTOptimizerFunction       int1e_spgsa01_optimizer;
extern CINTIntegralFunctionReal    int1e_spgsa01_cart;
extern CINTIntegralFunctionReal    int1e_spgsa01_sph;
extern CINTIntegralFunctionComplex int1e_spgsa01_spinor;

extern CINTOptimizerFunction       int1e_spgsp_optimizer;
extern CINTIntegralFunctionReal    int1e_spgsp_cart;
extern CINTIntegralFunctionReal    int1e_spgsp_sph;
extern CINTIntegralFunctionComplex int1e_spgsp_spinor;

extern CINTOptimizerFunction       int1e_spnuc_optimizer;
extern CINTIntegralFunctionReal    int1e_spnuc_cart;
extern CINTIntegralFunctionReal    int1e_spnuc_sph;
extern CINTIntegralFunctionComplex int1e_spnuc_spinor;

extern CINTOptimizerFunction       int1e_spnucsp_optimizer;
extern CINTIntegralFunctionReal    int1e_spnucsp_cart;
extern CINTIntegralFunctionReal    int1e_spnucsp_sph;
extern CINTIntegralFunctionComplex int1e_spnucsp_spinor;

extern CINTOptimizerFunction       int1e_sprinvsp_optimizer;
extern CINTIntegralFunctionReal    int1e_sprinvsp_cart;
extern CINTIntegralFunctionReal    int1e_sprinvsp_sph;
extern CINTIntegralFunctionComplex int1e_sprinvsp_spinor;

extern CINTOptimizerFunction       int1e_sprsp_optimizer;
extern CINTIntegralFunctionReal    int1e_sprsp_cart;
extern CINTIntegralFunctionReal    int1e_sprsp_sph;
extern CINTIntegralFunctionComplex int1e_sprsp_spinor;

extern CINTOptimizerFunction       int1e_spsigmasp_optimizer;
extern CINTIntegralFunctionReal    int1e_spsigmasp_cart;
extern CINTIntegralFunctionReal    int1e_spsigmasp_sph;
extern CINTIntegralFunctionComplex int1e_spsigmasp_spinor;

extern CINTOptimizerFunction       int1e_spsp_optimizer;
extern CINTIntegralFunctionReal    int1e_spsp_cart;
extern CINTIntegralFunctionReal    int1e_spsp_sph;
extern CINTIntegralFunctionComplex int1e_spsp_spinor;

extern CINTOptimizerFunction       int1e_spspsp_optimizer;
extern CINTIntegralFunctionReal    int1e_spspsp_cart;
extern CINTIntegralFunctionReal    int1e_spspsp_sph;
extern CINTIntegralFunctionComplex int1e_spspsp_spinor;

extern CINTOptimizerFunction       int1e_sr_optimizer;
extern CINTIntegralFunctionReal    int1e_sr_cart;
extern CINTIntegralFunctionReal    int1e_sr_sph;
extern CINTIntegralFunctionComplex int1e_sr_spinor;

extern CINTOptimizerFunction       int1e_srnucsr_optimizer;
extern CINTIntegralFunctionReal    int1e_srnucsr_cart;
extern CINTIntegralFunctionReal    int1e_srnucsr_sph;
extern CINTIntegralFunctionComplex int1e_srnucsr_spinor;

extern CINTOptimizerFunction       int1e_srsp_optimizer;
extern CINTIntegralFunctionReal    int1e_srsp_cart;
extern CINTIntegralFunctionReal    int1e_srsp_sph;
extern CINTIntegralFunctionComplex int1e_srsp_spinor;

extern CINTOptimizerFunction       int1e_srsr_optimizer;
extern CINTIntegralFunctionReal    int1e_srsr_cart;
extern CINTIntegralFunctionReal    int1e_srsr_sph;
extern CINTIntegralFunctionComplex int1e_srsr_spinor;

extern CINTOptimizerFunction       int1e_z_optimizer;
extern CINTIntegralFunctionReal    int1e_z_cart;
extern CINTIntegralFunctionReal    int1e_z_sph;
extern CINTIntegralFunctionComplex int1e_z_spinor;

extern CINTOptimizerFunction       int1e_z_origj_optimizer;
extern CINTIntegralFunctionReal    int1e_z_origj_cart;
extern CINTIntegralFunctionReal    int1e_z_origj_sph;
extern CINTIntegralFunctionComplex int1e_z_origj_spinor;

extern CINTOptimizerFunction       int1e_zz_optimizer;
extern CINTIntegralFunctionReal    int1e_zz_cart;
extern CINTIntegralFunctionReal    int1e_zz_sph;
extern CINTIntegralFunctionComplex int1e_zz_spinor;

extern CINTOptimizerFunction       int1e_zz_origj_optimizer;
extern CINTIntegralFunctionReal    int1e_zz_origj_cart;
extern CINTIntegralFunctionReal    int1e_zz_origj_sph;
extern CINTIntegralFunctionComplex int1e_zz_origj_spinor;

extern CINTOptimizerFunction       int2c2e_optimizer;
extern CINTIntegralFunctionReal    int2c2e_cart;
extern CINTIntegralFunctionReal    int2c2e_sph;
extern CINTIntegralFunctionComplex int2c2e_spinor;

extern CINTOptimizerFunction       int2c2e_ip1_optimizer;
extern CINTIntegralFunctionReal    int2c2e_ip1_cart;
extern CINTIntegralFunctionReal    int2c2e_ip1_sph;
extern CINTIntegralFunctionComplex int2c2e_ip1_spinor;

extern CINTOptimizerFunction       int2c2e_ip1ip2_optimizer;
extern CINTIntegralFunctionReal    int2c2e_ip1ip2_cart;
extern CINTIntegralFunctionReal    int2c2e_ip1ip2_sph;
extern CINTIntegralFunctionComplex int2c2e_ip1ip2_spinor;

extern CINTOptimizerFunction       int2c2e_ip2_optimizer;
extern CINTIntegralFunctionReal    int2c2e_ip2_cart;
extern CINTIntegralFunctionReal    int2c2e_ip2_sph;
extern CINTIntegralFunctionComplex int2c2e_ip2_spinor;

extern CINTOptimizerFunction       int2c2e_ipip1_optimizer;
extern CINTIntegralFunctionReal    int2c2e_ipip1_cart;
extern CINTIntegralFunctionReal    int2c2e_ipip1_sph;
extern CINTIntegralFunctionComplex int2c2e_ipip1_spinor;

extern CINTOptimizerFunction       int2e_optimizer;
extern CINTIntegralFunctionReal    int2e_cart;
extern CINTIntegralFunctionReal    int2e_sph;
extern CINTIntegralFunctionComplex int2e_spinor;

extern CINTOptimizerFunction       int2e_breit_r1p2_optimizer;
extern CINTIntegralFunctionReal    int2e_breit_r1p2_cart;
extern CINTIntegralFunctionReal    int2e_breit_r1p2_sph;
extern CINTIntegralFunctionComplex int2e_breit_r1p2_spinor;

extern CINTOptimizerFunction       int2e_breit_r2p2_optimizer;
extern CINTIntegralFunctionReal    int2e_breit_r2p2_cart;
extern CINTIntegralFunctionReal    int2e_breit_r2p2_sph;
extern CINTIntegralFunctionComplex int2e_breit_r2p2_spinor;

extern CINTOptimizerFunction       int2e_cg_sa10sp1_optimizer;
extern CINTIntegralFunctionReal    int2e_cg_sa10sp1_cart;
extern CINTIntegralFunctionReal    int2e_cg_sa10sp1_sph;
extern CINTIntegralFunctionComplex int2e_cg_sa10sp1_spinor;

extern CINTOptimizerFunction       int2e_cg_sa10sp1spsp2_optimizer;
extern CINTIntegralFunctionReal    int2e_cg_sa10sp1spsp2_cart;
extern CINTIntegralFunctionReal    int2e_cg_sa10sp1spsp2_sph;
extern CINTIntegralFunctionComplex int2e_cg_sa10sp1spsp2_spinor;

extern CINTOptimizerFunction       int2e_cg_ssa10ssp2_optimizer;
extern CINTIntegralFunctionReal    int2e_cg_ssa10ssp2_cart;
extern CINTIntegralFunctionReal    int2e_cg_ssa10ssp2_sph;
extern CINTIntegralFunctionComplex int2e_cg_ssa10ssp2_spinor;

extern CINTOptimizerFunction       int2e_g1_optimizer;
extern CINTIntegralFunctionReal    int2e_g1_cart;
extern CINTIntegralFunctionReal    int2e_g1_sph;
extern CINTIntegralFunctionComplex int2e_g1_spinor;

extern CINTOptimizerFunction       int2e_g1g2_optimizer;
extern CINTIntegralFunctionReal    int2e_g1g2_cart;
extern CINTIntegralFunctionReal    int2e_g1g2_sph;
extern CINTIntegralFunctionComplex int2e_g1g2_spinor;

extern CINTOptimizerFunction       int2e_g1spsp2_optimizer;
extern CINTIntegralFunctionReal    int2e_g1spsp2_cart;
extern CINTIntegralFunctionReal    int2e_g1spsp2_sph;
extern CINTIntegralFunctionComplex int2e_g1spsp2_spinor;

extern CINTOptimizerFunction       int2e_gauge_r1_sps1sps2_optimizer;
extern CINTIntegralFunctionReal    int2e_gauge_r1_sps1sps2_cart;
extern CINTIntegralFunctionReal    int2e_gauge_r1_sps1sps2_sph;
extern CINTIntegralFunctionComplex int2e_gauge_r1_sps1sps2_spinor;

extern CINTOptimizerFunction       int2e_gauge_r1_sps1ssp2_optimizer;
extern CINTIntegralFunctionReal    int2e_gauge_r1_sps1ssp2_cart;
extern CINTIntegralFunctionReal    int2e_gauge_r1_sps1ssp2_sph;
extern CINTIntegralFunctionComplex int2e_gauge_r1_sps1ssp2_spinor;

extern CINTOptimizerFunction       int2e_gauge_r1_ssp1sps2_optimizer;
extern CINTIntegralFunctionReal    int2e_gauge_r1_ssp1sps2_cart;
extern CINTIntegralFunctionReal    int2e_gauge_r1_ssp1sps2_sph;
extern CINTIntegralFunctionComplex int2e_gauge_r1_ssp1sps2_spinor;

extern CINTOptimizerFunction       int2e_gauge_r1_ssp1ssp2_optimizer;
extern CINTIntegralFunctionReal    int2e_gauge_r1_ssp1ssp2_cart;
extern CINTIntegralFunctionReal    int2e_gauge_r1_ssp1ssp2_sph;
extern CINTIntegralFunctionComplex int2e_gauge_r1_ssp1ssp2_spinor;

extern CINTOptimizerFunction       int2e_gauge_r2_sps1sps2_optimizer;
extern CINTIntegralFunctionReal    int2e_gauge_r2_sps1sps2_cart;
extern CINTIntegralFunctionReal    int2e_gauge_r2_sps1sps2_sph;
extern CINTIntegralFunctionComplex int2e_gauge_r2_sps1sps2_spinor;

extern CINTOptimizerFunction       int2e_gauge_r2_sps1ssp2_optimizer;
extern CINTIntegralFunctionReal    int2e_gauge_r2_sps1ssp2_cart;
extern CINTIntegralFunctionReal    int2e_gauge_r2_sps1ssp2_sph;
extern CINTIntegralFunctionComplex int2e_gauge_r2_sps1ssp2_spinor;

extern CINTOptimizerFunction       int2e_gauge_r2_ssp1sps2_optimizer;
extern CINTIntegralFunctionReal    int2e_gauge_r2_ssp1sps2_cart;
extern CINTIntegralFunctionReal    int2e_gauge_r2_ssp1sps2_sph;
extern CINTIntegralFunctionComplex int2e_gauge_r2_ssp1sps2_spinor;

extern CINTOptimizerFunction       int2e_gauge_r2_ssp1ssp2_optimizer;
extern CINTIntegralFunctionReal    int2e_gauge_r2_ssp1ssp2_cart;
extern CINTIntegralFunctionReal    int2e_gauge_r2_ssp1ssp2_sph;
extern CINTIntegralFunctionComplex int2e_gauge_r2_ssp1ssp2_spinor;

extern CINTOptimizerFunction       int2e_gg1_optimizer;
extern CINTIntegralFunctionReal    int2e_gg1_cart;
extern CINTIntegralFunctionReal    int2e_gg1_sph;
extern CINTIntegralFunctionComplex int2e_gg1_spinor;

extern CINTOptimizerFunction       int2e_giao_sa10sp1_optimizer;
extern CINTIntegralFunctionReal    int2e_giao_sa10sp1_cart;
extern CINTIntegralFunctionReal    int2e_giao_sa10sp1_sph;
extern CINTIntegralFunctionComplex int2e_giao_sa10sp1_spinor;

extern CINTOptimizerFunction       int2e_giao_sa10sp1spsp2_optimizer;
extern CINTIntegralFunctionReal    int2e_giao_sa10sp1spsp2_cart;
extern CINTIntegralFunctionReal    int2e_giao_sa10sp1spsp2_sph;
extern CINTIntegralFunctionComplex int2e_giao_sa10sp1spsp2_spinor;

extern CINTOptimizerFunction       int2e_giao_ssa10ssp2_optimizer;
extern CINTIntegralFunctionReal    int2e_giao_ssa10ssp2_cart;
extern CINTIntegralFunctionReal    int2e_giao_ssa10ssp2_sph;
extern CINTIntegralFunctionComplex int2e_giao_ssa10ssp2_spinor;

extern CINTOptimizerFunction       int2e_gssp1ssp2_optimizer;
extern CINTIntegralFunctionReal    int2e_gssp1ssp2_cart;
extern CINTIntegralFunctionReal    int2e_gssp1ssp2_sph;
extern CINTIntegralFunctionComplex int2e_gssp1ssp2_spinor;

extern CINTOptimizerFunction       int2e_ig1_optimizer;
extern CINTIntegralFunctionReal    int2e_ig1_cart;
extern CINTIntegralFunctionReal    int2e_ig1_sph;
extern CINTIntegralFunctionComplex int2e_ig1_spinor;

extern CINTOptimizerFunction       int2e_ip1_optimizer;
extern CINTIntegralFunctionReal    int2e_ip1_cart;
extern CINTIntegralFunctionReal    int2e_ip1_sph;
extern CINTIntegralFunctionComplex int2e_ip1_spinor;

extern CINTOptimizerFunction       int2e_ip1ip2_optimizer;
extern CINTIntegralFunctionReal    int2e_ip1ip2_cart;
extern CINTIntegralFunctionReal    int2e_ip1ip2_sph;
extern CINTIntegralFunctionComplex int2e_ip1ip2_spinor;

extern CINTOptimizerFunction       int2e_ip1spsp2_optimizer;
extern CINTIntegralFunctionReal    int2e_ip1spsp2_cart;
extern CINTIntegralFunctionReal    int2e_ip1spsp2_sph;
extern CINTIntegralFunctionComplex int2e_ip1spsp2_spinor;

extern CINTOptimizerFunction       int2e_ip1srsr2_optimizer;
extern CINTIntegralFunctionReal    int2e_ip1srsr2_cart;
extern CINTIntegralFunctionReal    int2e_ip1srsr2_sph;
extern CINTIntegralFunctionComplex int2e_ip1srsr2_spinor;

extern CINTOptimizerFunction       int2e_ip1v_r1_optimizer;
extern CINTIntegralFunctionReal    int2e_ip1v_r1_cart;
extern CINTIntegralFunctionReal    int2e_ip1v_r1_sph;
extern CINTIntegralFunctionComplex int2e_ip1v_r1_spinor;

extern CINTOptimizerFunction       int2e_ip1v_rc1_optimizer;
extern CINTIntegralFunctionReal    int2e_ip1v_rc1_cart;
extern CINTIntegralFunctionReal    int2e_ip1v_rc1_sph;
extern CINTIntegralFunctionComplex int2e_ip1v_rc1_spinor;

extern CINTOptimizerFunction       int2e_ip2_optimizer;
extern CINTIntegralFunctionReal    int2e_ip2_cart;
extern CINTIntegralFunctionReal    int2e_ip2_sph;
extern CINTIntegralFunctionComplex int2e_ip2_spinor;

extern CINTOptimizerFunction       int2e_ipip1_optimizer;
extern CINTIntegralFunctionReal    int2e_ipip1_cart;
extern CINTIntegralFunctionReal    int2e_ipip1_sph;
extern CINTIntegralFunctionComplex int2e_ipip1_spinor;

extern CINTOptimizerFunction       int2e_ipip1ipip2_optimizer;
extern CINTIntegralFunctionReal    int2e_ipip1ipip2_cart;
extern CINTIntegralFunctionReal    int2e_ipip1ipip2_sph;
extern CINTIntegralFunctionComplex int2e_ipip1ipip2_spinor;

extern CINTOptimizerFunction       int2e_ipspsp1_optimizer;
extern CINTIntegralFunctionReal    int2e_ipspsp1_cart;
extern CINTIntegralFunctionReal    int2e_ipspsp1_sph;
extern CINTIntegralFunctionComplex int2e_ipspsp1_spinor;

extern CINTOptimizerFunction       int2e_ipspsp1spsp2_optimizer;
extern CINTIntegralFunctionReal    int2e_ipspsp1spsp2_cart;
extern CINTIntegralFunctionReal    int2e_ipspsp1spsp2_sph;
extern CINTIntegralFunctionComplex int2e_ipspsp1spsp2_spinor;

extern CINTOptimizerFunction       int2e_ipsrsr1_optimizer;
extern CINTIntegralFunctionReal    int2e_ipsrsr1_cart;
extern CINTIntegralFunctionReal    int2e_ipsrsr1_sph;
extern CINTIntegralFunctionComplex int2e_ipsrsr1_spinor;

extern CINTOptimizerFunction       int2e_ipsrsr1srsr2_optimizer;
extern CINTIntegralFunctionReal    int2e_ipsrsr1srsr2_cart;
extern CINTIntegralFunctionReal    int2e_ipsrsr1srsr2_sph;
extern CINTIntegralFunctionComplex int2e_ipsrsr1srsr2_spinor;

extern CINTOptimizerFunction       int2e_ipvg1_xp1_optimizer;
extern CINTIntegralFunctionReal    int2e_ipvg1_xp1_cart;
extern CINTIntegralFunctionReal    int2e_ipvg1_xp1_sph;
extern CINTIntegralFunctionComplex int2e_ipvg1_xp1_spinor;

extern CINTOptimizerFunction       int2e_ipvg2_xp1_optimizer;
extern CINTIntegralFunctionReal    int2e_ipvg2_xp1_cart;
extern CINTIntegralFunctionReal    int2e_ipvg2_xp1_sph;
extern CINTIntegralFunctionComplex int2e_ipvg2_xp1_spinor;

extern CINTOptimizerFunction       int2e_ipvip1_optimizer;
extern CINTIntegralFunctionReal    int2e_ipvip1_cart;
extern CINTIntegralFunctionReal    int2e_ipvip1_sph;
extern CINTIntegralFunctionComplex int2e_ipvip1_spinor;

extern CINTOptimizerFunction       int2e_ipvip1ipvip2_optimizer;
extern CINTIntegralFunctionReal    int2e_ipvip1ipvip2_cart;
extern CINTIntegralFunctionReal    int2e_ipvip1ipvip2_sph;
extern CINTIntegralFunctionComplex int2e_ipvip1ipvip2_spinor;

extern CINTOptimizerFunction       int2e_p1vxp1_optimizer;
extern CINTIntegralFunctionReal    int2e_p1vxp1_cart;
extern CINTIntegralFunctionReal    int2e_p1vxp1_sph;
extern CINTIntegralFunctionComplex int2e_p1vxp1_spinor;

extern CINTOptimizerFunction       int2e_pp1_optimizer;
extern CINTIntegralFunctionReal    int2e_pp1_cart;
extern CINTIntegralFunctionReal    int2e_pp1_sph;
extern CINTIntegralFunctionComplex int2e_pp1_spinor;

extern CINTOptimizerFunction       int2e_pp1pp2_optimizer;
extern CINTIntegralFunctionReal    int2e_pp1pp2_cart;
extern CINTIntegralFunctionReal    int2e_pp1pp2_sph;
extern CINTIntegralFunctionComplex int2e_pp1pp2_spinor;

extern CINTOptimizerFunction       int2e_pp2_optimizer;
extern CINTIntegralFunctionReal    int2e_pp2_cart;
extern CINTIntegralFunctionReal    int2e_pp2_sph;
extern CINTIntegralFunctionComplex int2e_pp2_spinor;

extern CINTOptimizerFunction       int2e_spgsp1_optimizer;
extern CINTIntegralFunctionReal    int2e_spgsp1_cart;
extern CINTIntegralFunctionReal    int2e_spgsp1_sph;
extern CINTIntegralFunctionComplex int2e_spgsp1_spinor;

extern CINTOptimizerFunction       int2e_spgsp1spsp2_optimizer;
extern CINTIntegralFunctionReal    int2e_spgsp1spsp2_cart;
extern CINTIntegralFunctionReal    int2e_spgsp1spsp2_sph;
extern CINTIntegralFunctionComplex int2e_spgsp1spsp2_spinor;

extern CINTOptimizerFunction       int2e_sps1sps2_optimizer;
extern CINTIntegralFunctionReal    int2e_sps1sps2_cart;
extern CINTIntegralFunctionReal    int2e_sps1sps2_sph;
extern CINTIntegralFunctionComplex int2e_sps1sps2_spinor;

extern CINTOptimizerFunction       int2e_sps1ssp2_optimizer;
extern CINTIntegralFunctionReal    int2e_sps1ssp2_cart;
extern CINTIntegralFunctionReal    int2e_sps1ssp2_sph;
extern CINTIntegralFunctionComplex int2e_sps1ssp2_spinor;

extern CINTOptimizerFunction       int2e_spsp1_optimizer;
extern CINTIntegralFunctionReal    int2e_spsp1_cart;
extern CINTIntegralFunctionReal    int2e_spsp1_sph;
extern CINTIntegralFunctionComplex int2e_spsp1_spinor;

extern CINTOptimizerFunction       int2e_spsp1spsp2_optimizer;
extern CINTIntegralFunctionReal    int2e_spsp1spsp2_cart;
extern CINTIntegralFunctionReal    int2e_spsp1spsp2_sph;
extern CINTIntegralFunctionComplex int2e_spsp1spsp2_spinor;

extern CINTOptimizerFunction       int2e_spsp2_optimizer;
extern CINTIntegralFunctionReal    int2e_spsp2_cart;
extern CINTIntegralFunctionReal    int2e_spsp2_sph;
extern CINTIntegralFunctionComplex int2e_spsp2_spinor;

extern CINTOptimizerFunction       int2e_spv1_optimizer;
extern CINTIntegralFunctionReal    int2e_spv1_cart;
extern CINTIntegralFunctionReal    int2e_spv1_sph;
extern CINTIntegralFunctionComplex int2e_spv1_spinor;

extern CINTOptimizerFunction       int2e_spv1spsp2_optimizer;
extern CINTIntegralFunctionReal    int2e_spv1spsp2_cart;
extern CINTIntegralFunctionReal    int2e_spv1spsp2_sph;
extern CINTIntegralFunctionComplex int2e_spv1spsp2_spinor;

extern CINTOptimizerFunction       int2e_spv1spv2_optimizer;
extern CINTIntegralFunctionReal    int2e_spv1spv2_cart;
extern CINTIntegralFunctionReal    int2e_spv1spv2_sph;
extern CINTIntegralFunctionComplex int2e_spv1spv2_spinor;

extern CINTOptimizerFunction       int2e_spv1vsp2_optimizer;
extern CINTIntegralFunctionReal    int2e_spv1vsp2_cart;
extern CINTIntegralFunctionReal    int2e_spv1vsp2_sph;
extern CINTIntegralFunctionComplex int2e_spv1vsp2_spinor;

extern CINTOptimizerFunction       int2e_srsr1_optimizer;
extern CINTIntegralFunctionReal    int2e_srsr1_cart;
extern CINTIntegralFunctionReal    int2e_srsr1_sph;
extern CINTIntegralFunctionComplex int2e_srsr1_spinor;

extern CINTOptimizerFunction       int2e_srsr1srsr2_optimizer;
extern CINTIntegralFunctionReal    int2e_srsr1srsr2_cart;
extern CINTIntegralFunctionReal    int2e_srsr1srsr2_sph;
extern CINTIntegralFunctionComplex int2e_srsr1srsr2_spinor;

extern CINTOptimizerFunction       int2e_ssp1sps2_optimizer;
extern CINTIntegralFunctionReal    int2e_ssp1sps2_cart;
extern CINTIntegralFunctionReal    int2e_ssp1sps2_sph;
extern CINTIntegralFunctionComplex int2e_ssp1sps2_spinor;

extern CINTOptimizerFunction       int2e_ssp1ssp2_optimizer;
extern CINTIntegralFunctionReal    int2e_ssp1ssp2_cart;
extern CINTIntegralFunctionReal    int2e_ssp1ssp2_sph;
extern CINTIntegralFunctionComplex int2e_ssp1ssp2_spinor;

extern CINTOptimizerFunction       int2e_stg_optimizer;
extern CINTIntegralFunctionReal    int2e_stg_cart;
extern CINTIntegralFunctionReal    int2e_stg_sph;
extern CINTIntegralFunctionComplex int2e_stg_spinor;

extern CINTOptimizerFunction       int2e_stg_ip1_optimizer;
extern CINTIntegralFunctionReal    int2e_stg_ip1_cart;
extern CINTIntegralFunctionReal    int2e_stg_ip1_sph;
extern CINTIntegralFunctionComplex int2e_stg_ip1_spinor;

extern CINTOptimizerFunction       int2e_stg_ip1ip2_optimizer;
extern CINTIntegralFunctionReal    int2e_stg_ip1ip2_cart;
extern CINTIntegralFunctionReal    int2e_stg_ip1ip2_sph;
extern CINTIntegralFunctionComplex int2e_stg_ip1ip2_spinor;

extern CINTOptimizerFunction       int2e_stg_ipip1_optimizer;
extern CINTIntegralFunctionReal    int2e_stg_ipip1_cart;
extern CINTIntegralFunctionReal    int2e_stg_ipip1_sph;
extern CINTIntegralFunctionComplex int2e_stg_ipip1_spinor;

extern CINTOptimizerFunction       int2e_stg_ipvip1_optimizer;
extern CINTIntegralFunctionReal    int2e_stg_ipvip1_cart;
extern CINTIntegralFunctionReal    int2e_stg_ipvip1_sph;
extern CINTIntegralFunctionComplex int2e_stg_ipvip1_spinor;

extern CINTOptimizerFunction       int2e_vsp1_optimizer;
extern CINTIntegralFunctionReal    int2e_vsp1_cart;
extern CINTIntegralFunctionReal    int2e_vsp1_sph;
extern CINTIntegralFunctionComplex int2e_vsp1_spinor;

extern CINTOptimizerFunction       int2e_vsp1spsp2_optimizer;
extern CINTIntegralFunctionReal    int2e_vsp1spsp2_cart;
extern CINTIntegralFunctionReal    int2e_vsp1spsp2_sph;
extern CINTIntegralFunctionComplex int2e_vsp1spsp2_spinor;

extern CINTOptimizerFunction       int2e_vsp1spv2_optimizer;
extern CINTIntegralFunctionReal    int2e_vsp1spv2_cart;
extern CINTIntegralFunctionReal    int2e_vsp1spv2_sph;
extern CINTIntegralFunctionComplex int2e_vsp1spv2_spinor;

extern CINTOptimizerFunction       int2e_vsp1vsp2_optimizer;
extern CINTIntegralFunctionReal    int2e_vsp1vsp2_cart;
extern CINTIntegralFunctionReal    int2e_vsp1vsp2_sph;
extern CINTIntegralFunctionComplex int2e_vsp1vsp2_spinor;

extern CINTOptimizerFunction       int2e_yp_optimizer;
extern CINTIntegralFunctionReal    int2e_yp_cart;
extern CINTIntegralFunctionReal    int2e_yp_sph;
extern CINTIntegralFunctionComplex int2e_yp_spinor;

extern CINTOptimizerFunction       int2e_yp_ip1_optimizer;
extern CINTIntegralFunctionReal    int2e_yp_ip1_cart;
extern CINTIntegralFunctionReal    int2e_yp_ip1_sph;
extern CINTIntegralFunctionComplex int2e_yp_ip1_spinor;

extern CINTOptimizerFunction       int2e_yp_ip1ip2_optimizer;
extern CINTIntegralFunctionReal    int2e_yp_ip1ip2_cart;
extern CINTIntegralFunctionReal    int2e_yp_ip1ip2_sph;
extern CINTIntegralFunctionComplex int2e_yp_ip1ip2_spinor;

extern CINTOptimizerFunction       int2e_yp_ipip1_optimizer;
extern CINTIntegralFunctionReal    int2e_yp_ipip1_cart;
extern CINTIntegralFunctionReal    int2e_yp_ipip1_sph;
extern CINTIntegralFunctionComplex int2e_yp_ipip1_spinor;

extern CINTOptimizerFunction       int2e_yp_ipvip1_optimizer;
extern CINTIntegralFunctionReal    int2e_yp_ipvip1_cart;
extern CINTIntegralFunctionReal    int2e_yp_ipvip1_sph;
extern CINTIntegralFunctionComplex int2e_yp_ipvip1_spinor;

extern CINTOptimizerFunction       int3c1e_optimizer;
extern CINTIntegralFunctionReal    int3c1e_cart;
extern CINTIntegralFunctionReal    int3c1e_sph;
extern CINTIntegralFunctionComplex int3c1e_spinor;

extern CINTOptimizerFunction       int3c1e_ip1_optimizer;
extern CINTIntegralFunctionReal    int3c1e_ip1_cart;
extern CINTIntegralFunctionReal    int3c1e_ip1_sph;
extern CINTIntegralFunctionComplex int3c1e_ip1_spinor;

extern CINTOptimizerFunction       int3c1e_ip1_r2_origk_optimizer;
extern CINTIntegralFunctionReal    int3c1e_ip1_r2_origk_cart;
extern CINTIntegralFunctionReal    int3c1e_ip1_r2_origk_sph;
extern CINTIntegralFunctionComplex int3c1e_ip1_r2_origk_spinor;

extern CINTOptimizerFunction       int3c1e_ip1_r4_origk_optimizer;
extern CINTIntegralFunctionReal    int3c1e_ip1_r4_origk_cart;
extern CINTIntegralFunctionReal    int3c1e_ip1_r4_origk_sph;
extern CINTIntegralFunctionComplex int3c1e_ip1_r4_origk_spinor;

extern CINTOptimizerFunction       int3c1e_ip1_r6_origk_optimizer;
extern CINTIntegralFunctionReal    int3c1e_ip1_r6_origk_cart;
extern CINTIntegralFunctionReal    int3c1e_ip1_r6_origk_sph;
extern CINTIntegralFunctionComplex int3c1e_ip1_r6_origk_spinor;

extern CINTOptimizerFunction       int3c1e_iprinv_optimizer;
extern CINTIntegralFunctionReal    int3c1e_iprinv_cart;
extern CINTIntegralFunctionReal    int3c1e_iprinv_sph;
extern CINTIntegralFunctionComplex int3c1e_iprinv_spinor;

extern CINTOptimizerFunction       int3c1e_p2_optimizer;
extern CINTIntegralFunctionReal    int3c1e_p2_cart;
extern CINTIntegralFunctionReal    int3c1e_p2_sph;
extern CINTIntegralFunctionComplex int3c1e_p2_spinor;

extern CINTOptimizerFunction       int3c1e_r2_origk_optimizer;
extern CINTIntegralFunctionReal    int3c1e_r2_origk_cart;
extern CINTIntegralFunctionReal    int3c1e_r2_origk_sph;
extern CINTIntegralFunctionComplex int3c1e_r2_origk_spinor;

extern CINTOptimizerFunction       int3c1e_r4_origk_optimizer;
extern CINTIntegralFunctionReal    int3c1e_r4_origk_cart;
extern CINTIntegralFunctionReal    int3c1e_r4_origk_sph;
extern CINTIntegralFunctionComplex int3c1e_r4_origk_spinor;

extern CINTOptimizerFunction       int3c1e_r6_origk_optimizer;
extern CINTIntegralFunctionReal    int3c1e_r6_origk_cart;
extern CINTIntegralFunctionReal    int3c1e_r6_origk_sph;
extern CINTIntegralFunctionComplex int3c1e_r6_origk_spinor;

extern CINTOptimizerFunction       int3c1e_rinv_optimizer;
extern CINTIntegralFunctionReal    int3c1e_rinv_cart;
extern CINTIntegralFunctionReal    int3c1e_rinv_sph;
extern CINTIntegralFunctionComplex int3c1e_rinv_spinor;

extern CINTOptimizerFunction       int3c2e_optimizer;
extern CINTIntegralFunctionReal    int3c2e_cart;
extern CINTIntegralFunctionReal    int3c2e_sph;
extern CINTIntegralFunctionComplex int3c2e_spinor;

extern CINTOptimizerFunction       int3c2e_ig1_optimizer;
extern CINTIntegralFunctionReal    int3c2e_ig1_cart;
extern CINTIntegralFunctionReal    int3c2e_ig1_sph;
extern CINTIntegralFunctionComplex int3c2e_ig1_spinor;

extern CINTOptimizerFunction       int3c2e_ip1_optimizer;
extern CINTIntegralFunctionReal    int3c2e_ip1_cart;
extern CINTIntegralFunctionReal    int3c2e_ip1_sph;
extern CINTIntegralFunctionComplex int3c2e_ip1_spinor;

extern CINTOptimizerFunction       int3c2e_ip1ip2_optimizer;
extern CINTIntegralFunctionReal    int3c2e_ip1ip2_cart;
extern CINTIntegralFunctionReal    int3c2e_ip1ip2_sph;
extern CINTIntegralFunctionComplex int3c2e_ip1ip2_spinor;

extern CINTOptimizerFunction       int3c2e_ip2_optimizer;
extern CINTIntegralFunctionReal    int3c2e_ip2_cart;
extern CINTIntegralFunctionReal    int3c2e_ip2_sph;
extern CINTIntegralFunctionComplex int3c2e_ip2_spinor;

extern CINTOptimizerFunction       int3c2e_ipip1_optimizer;
extern CINTIntegralFunctionReal    int3c2e_ipip1_cart;
extern CINTIntegralFunctionReal    int3c2e_ipip1_sph;
extern CINTIntegralFunctionComplex int3c2e_ipip1_spinor;

extern CINTOptimizerFunction       int3c2e_ipip2_optimizer;
extern CINTIntegralFunctionReal    int3c2e_ipip2_cart;
extern CINTIntegralFunctionReal    int3c2e_ipip2_sph;
extern CINTIntegralFunctionComplex int3c2e_ipip2_spinor;

extern CINTOptimizerFunction       int3c2e_ipspsp1_optimizer;
extern CINTIntegralFunctionReal    int3c2e_ipspsp1_cart;
extern CINTIntegralFunctionReal    int3c2e_ipspsp1_sph;
extern CINTIntegralFunctionComplex int3c2e_ipspsp1_spinor;

extern CINTOptimizerFunction       int3c2e_ipvip1_optimizer;
extern CINTIntegralFunctionReal    int3c2e_ipvip1_cart;
extern CINTIntegralFunctionReal    int3c2e_ipvip1_sph;
extern CINTIntegralFunctionComplex int3c2e_ipvip1_spinor;

extern CINTOptimizerFunction       int3c2e_pvp1_optimizer;
extern CINTIntegralFunctionReal    int3c2e_pvp1_cart;
extern CINTIntegralFunctionReal    int3c2e_pvp1_sph;
extern CINTIntegralFunctionComplex int3c2e_pvp1_spinor;

extern CINTOptimizerFunction       int3c2e_pvxp1_optimizer;
extern CINTIntegralFunctionReal    int3c2e_pvxp1_cart;
extern CINTIntegralFunctionReal    int3c2e_pvxp1_sph;
extern CINTIntegralFunctionComplex int3c2e_pvxp1_spinor;

extern CINTOptimizerFunction       int3c2e_spsp1_optimizer;
extern CINTIntegralFunctionReal    int3c2e_spsp1_cart;
extern CINTIntegralFunctionReal    int3c2e_spsp1_sph;
extern CINTIntegralFunctionComplex int3c2e_spsp1_spinor;

extern CINTOptimizerFunction       int3c2e_spsp1ip2_optimizer;
extern CINTIntegralFunctionReal    int3c2e_spsp1ip2_cart;
extern CINTIntegralFunctionReal    int3c2e_spsp1ip2_sph;
extern CINTIntegralFunctionComplex int3c2e_spsp1ip2_spinor;

extern CINTOptimizerFunction       int4c1e_optimizer;
extern CINTIntegralFunctionReal    int4c1e_cart;
extern CINTIntegralFunctionReal    int4c1e_sph;
extern CINTIntegralFunctionComplex int4c1e_spinor;

extern CINTOptimizerFunction       cint1e_a01gp_optimizer;
extern CINTIntegralFunctionReal    cint1e_a01gp_cart;
extern CINTIntegralFunctionReal    cint1e_a01gp_sph;
extern CINTIntegralFunctionComplex cint1e_a01gp_spinor;

extern CINTOptimizerFunction       cint1e_cg_a11part_optimizer;
extern CINTIntegralFunctionReal    cint1e_cg_a11part_cart;
extern CINTIntegralFunctionReal    cint1e_cg_a11part_sph;
extern CINTIntegralFunctionComplex cint1e_cg_a11part_spinor;

extern CINTOptimizerFunction       cint1e_cg_irxp_optimizer;
extern CINTIntegralFunctionReal    cint1e_cg_irxp_cart;
extern CINTIntegralFunctionReal    cint1e_cg_irxp_sph;
extern CINTIntegralFunctionComplex cint1e_cg_irxp_spinor;

extern CINTOptimizerFunction       cint1e_cg_sa10nucsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_cg_sa10nucsp_cart;
extern CINTIntegralFunctionReal    cint1e_cg_sa10nucsp_sph;
extern CINTIntegralFunctionComplex cint1e_cg_sa10nucsp_spinor;

extern CINTOptimizerFunction       cint1e_cg_sa10sa01_optimizer;
extern CINTIntegralFunctionReal    cint1e_cg_sa10sa01_cart;
extern CINTIntegralFunctionReal    cint1e_cg_sa10sa01_sph;
extern CINTIntegralFunctionComplex cint1e_cg_sa10sa01_spinor;

extern CINTOptimizerFunction       cint1e_cg_sa10sp_optimizer;
extern CINTIntegralFunctionReal    cint1e_cg_sa10sp_cart;
extern CINTIntegralFunctionReal    cint1e_cg_sa10sp_sph;
extern CINTIntegralFunctionComplex cint1e_cg_sa10sp_spinor;

extern CINTOptimizerFunction       cint1e_drinv_optimizer;
extern CINTIntegralFunctionReal    cint1e_drinv_cart;
extern CINTIntegralFunctionReal    cint1e_drinv_sph;
extern CINTIntegralFunctionComplex cint1e_drinv_spinor;

extern CINTOptimizerFunction       cint1e_ggkin_optimizer;
extern CINTIntegralFunctionReal    cint1e_ggkin_cart;
extern CINTIntegralFunctionReal    cint1e_ggkin_sph;
extern CINTIntegralFunctionComplex cint1e_ggkin_spinor;

extern CINTOptimizerFunction       cint1e_ggnuc_optimizer;
extern CINTIntegralFunctionReal    cint1e_ggnuc_cart;
extern CINTIntegralFunctionReal    cint1e_ggnuc_sph;
extern CINTIntegralFunctionComplex cint1e_ggnuc_spinor;

extern CINTOptimizerFunction       cint1e_ggovlp_optimizer;
extern CINTIntegralFunctionReal    cint1e_ggovlp_cart;
extern CINTIntegralFunctionReal    cint1e_ggovlp_sph;
extern CINTIntegralFunctionComplex cint1e_ggovlp_spinor;

extern CINTOptimizerFunction       cint1e_giao_a11part_optimizer;
extern CINTIntegralFunctionReal    cint1e_giao_a11part_cart;
extern CINTIntegralFunctionReal    cint1e_giao_a11part_sph;
extern CINTIntegralFunctionComplex cint1e_giao_a11part_spinor;

extern CINTOptimizerFunction       cint1e_giao_irjxp_optimizer;
extern CINTIntegralFunctionReal    cint1e_giao_irjxp_cart;
extern CINTIntegralFunctionReal    cint1e_giao_irjxp_sph;
extern CINTIntegralFunctionComplex cint1e_giao_irjxp_spinor;

extern CINTOptimizerFunction       cint1e_giao_sa10nucsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_giao_sa10nucsp_cart;
extern CINTIntegralFunctionReal    cint1e_giao_sa10nucsp_sph;
extern CINTIntegralFunctionComplex cint1e_giao_sa10nucsp_spinor;

extern CINTOptimizerFunction       cint1e_giao_sa10sa01_optimizer;
extern CINTIntegralFunctionReal    cint1e_giao_sa10sa01_cart;
extern CINTIntegralFunctionReal    cint1e_giao_sa10sa01_sph;
extern CINTIntegralFunctionComplex cint1e_giao_sa10sa01_spinor;

extern CINTOptimizerFunction       cint1e_giao_sa10sp_optimizer;
extern CINTIntegralFunctionReal    cint1e_giao_sa10sp_cart;
extern CINTIntegralFunctionReal    cint1e_giao_sa10sp_sph;
extern CINTIntegralFunctionComplex cint1e_giao_sa10sp_spinor;

extern CINTOptimizerFunction       cint1e_gnuc_optimizer;
extern CINTIntegralFunctionReal    cint1e_gnuc_cart;
extern CINTIntegralFunctionReal    cint1e_gnuc_sph;
extern CINTIntegralFunctionComplex cint1e_gnuc_spinor;

extern CINTOptimizerFunction       cint1e_govlp_optimizer;
extern CINTIntegralFunctionReal    cint1e_govlp_cart;
extern CINTIntegralFunctionReal    cint1e_govlp_sph;
extern CINTIntegralFunctionComplex cint1e_govlp_spinor;

extern CINTOptimizerFunction       cint1e_grids_optimizer;
extern CINTIntegralFunctionReal    cint1e_grids_cart;
extern CINTIntegralFunctionReal    cint1e_grids_sph;
extern CINTIntegralFunctionComplex cint1e_grids_spinor;

extern CINTOptimizerFunction       cint1e_grids_ip_optimizer;
extern CINTIntegralFunctionReal    cint1e_grids_ip_cart;
extern CINTIntegralFunctionReal    cint1e_grids_ip_sph;
extern CINTIntegralFunctionComplex cint1e_grids_ip_spinor;

extern CINTOptimizerFunction       cint1e_grids_ipip_optimizer;
extern CINTIntegralFunctionReal    cint1e_grids_ipip_cart;
extern CINTIntegralFunctionReal    cint1e_grids_ipip_sph;
extern CINTIntegralFunctionComplex cint1e_grids_ipip_spinor;

extern CINTOptimizerFunction       cint1e_grids_ipvip_optimizer;
extern CINTIntegralFunctionReal    cint1e_grids_ipvip_cart;
extern CINTIntegralFunctionReal    cint1e_grids_ipvip_sph;
extern CINTIntegralFunctionComplex cint1e_grids_ipvip_spinor;

extern CINTOptimizerFunction       cint1e_grids_spvsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_grids_spvsp_cart;
extern CINTIntegralFunctionReal    cint1e_grids_spvsp_sph;
extern CINTIntegralFunctionComplex cint1e_grids_spvsp_spinor;

extern CINTOptimizerFunction       cint1e_grjxp_optimizer;
extern CINTIntegralFunctionReal    cint1e_grjxp_cart;
extern CINTIntegralFunctionReal    cint1e_grjxp_sph;
extern CINTIntegralFunctionComplex cint1e_grjxp_spinor;

extern CINTOptimizerFunction       cint1e_ia01p_optimizer;
extern CINTIntegralFunctionReal    cint1e_ia01p_cart;
extern CINTIntegralFunctionReal    cint1e_ia01p_sph;
extern CINTIntegralFunctionComplex cint1e_ia01p_spinor;

extern CINTOptimizerFunction       cint1e_igkin_optimizer;
extern CINTIntegralFunctionReal    cint1e_igkin_cart;
extern CINTIntegralFunctionReal    cint1e_igkin_sph;
extern CINTIntegralFunctionComplex cint1e_igkin_spinor;

extern CINTOptimizerFunction       cint1e_ignuc_optimizer;
extern CINTIntegralFunctionReal    cint1e_ignuc_cart;
extern CINTIntegralFunctionReal    cint1e_ignuc_sph;
extern CINTIntegralFunctionComplex cint1e_ignuc_spinor;

extern CINTOptimizerFunction       cint1e_igovlp_optimizer;
extern CINTIntegralFunctionReal    cint1e_igovlp_cart;
extern CINTIntegralFunctionReal    cint1e_igovlp_sph;
extern CINTIntegralFunctionComplex cint1e_igovlp_spinor;

extern CINTOptimizerFunction       cint1e_inuc_rcxp_optimizer;
extern CINTIntegralFunctionReal    cint1e_inuc_rcxp_cart;
extern CINTIntegralFunctionReal    cint1e_inuc_rcxp_sph;
extern CINTIntegralFunctionComplex cint1e_inuc_rcxp_spinor;

extern CINTOptimizerFunction       cint1e_inuc_rxp_optimizer;
extern CINTIntegralFunctionReal    cint1e_inuc_rxp_cart;
extern CINTIntegralFunctionReal    cint1e_inuc_rxp_sph;
extern CINTIntegralFunctionComplex cint1e_inuc_rxp_spinor;

extern CINTOptimizerFunction       cint1e_ipipipiprinv_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipipipiprinv_cart;
extern CINTIntegralFunctionReal    cint1e_ipipipiprinv_sph;
extern CINTIntegralFunctionComplex cint1e_ipipipiprinv_spinor;

extern CINTOptimizerFunction       cint1e_ipipipnuc_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipipipnuc_cart;
extern CINTIntegralFunctionReal    cint1e_ipipipnuc_sph;
extern CINTIntegralFunctionComplex cint1e_ipipipnuc_spinor;

extern CINTOptimizerFunction       cint1e_ipipiprinv_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipipiprinv_cart;
extern CINTIntegralFunctionReal    cint1e_ipipiprinv_sph;
extern CINTIntegralFunctionComplex cint1e_ipipiprinv_spinor;

extern CINTOptimizerFunction       cint1e_ipipiprinvip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipipiprinvip_cart;
extern CINTIntegralFunctionReal    cint1e_ipipiprinvip_sph;
extern CINTIntegralFunctionComplex cint1e_ipipiprinvip_spinor;

extern CINTOptimizerFunction       cint1e_ipipkin_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipipkin_cart;
extern CINTIntegralFunctionReal    cint1e_ipipkin_sph;
extern CINTIntegralFunctionComplex cint1e_ipipkin_spinor;

extern CINTOptimizerFunction       cint1e_ipipnuc_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipipnuc_cart;
extern CINTIntegralFunctionReal    cint1e_ipipnuc_sph;
extern CINTIntegralFunctionComplex cint1e_ipipnuc_spinor;

extern CINTOptimizerFunction       cint1e_ipipnucip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipipnucip_cart;
extern CINTIntegralFunctionReal    cint1e_ipipnucip_sph;
extern CINTIntegralFunctionComplex cint1e_ipipnucip_spinor;

extern CINTOptimizerFunction       cint1e_ipipovlp_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipipovlp_cart;
extern CINTIntegralFunctionReal    cint1e_ipipovlp_sph;
extern CINTIntegralFunctionComplex cint1e_ipipovlp_spinor;

extern CINTOptimizerFunction       cint1e_ipippnucp_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipippnucp_cart;
extern CINTIntegralFunctionReal    cint1e_ipippnucp_sph;
extern CINTIntegralFunctionComplex cint1e_ipippnucp_spinor;

extern CINTOptimizerFunction       cint1e_ipipprinvp_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipipprinvp_cart;
extern CINTIntegralFunctionReal    cint1e_ipipprinvp_sph;
extern CINTIntegralFunctionComplex cint1e_ipipprinvp_spinor;

extern CINTOptimizerFunction       cint1e_ipipr_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipipr_cart;
extern CINTIntegralFunctionReal    cint1e_ipipr_sph;
extern CINTIntegralFunctionComplex cint1e_ipipr_spinor;

extern CINTOptimizerFunction       cint1e_ipiprinv_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipiprinv_cart;
extern CINTIntegralFunctionReal    cint1e_ipiprinv_sph;
extern CINTIntegralFunctionComplex cint1e_ipiprinv_spinor;

extern CINTOptimizerFunction       cint1e_ipiprinvip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipiprinvip_cart;
extern CINTIntegralFunctionReal    cint1e_ipiprinvip_sph;
extern CINTIntegralFunctionComplex cint1e_ipiprinvip_spinor;

extern CINTOptimizerFunction       cint1e_ipiprinvipip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipiprinvipip_cart;
extern CINTIntegralFunctionReal    cint1e_ipiprinvipip_sph;
extern CINTIntegralFunctionComplex cint1e_ipiprinvipip_spinor;

extern CINTOptimizerFunction       cint1e_ipiprinvrip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipiprinvrip_cart;
extern CINTIntegralFunctionReal    cint1e_ipiprinvrip_sph;
extern CINTIntegralFunctionComplex cint1e_ipiprinvrip_spinor;

extern CINTOptimizerFunction       cint1e_ipipspnucsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipipspnucsp_cart;
extern CINTIntegralFunctionReal    cint1e_ipipspnucsp_sph;
extern CINTIntegralFunctionComplex cint1e_ipipspnucsp_spinor;

extern CINTOptimizerFunction       cint1e_ipipsprinvsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipipsprinvsp_cart;
extern CINTIntegralFunctionReal    cint1e_ipipsprinvsp_sph;
extern CINTIntegralFunctionComplex cint1e_ipipsprinvsp_spinor;

extern CINTOptimizerFunction       cint1e_ipkin_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipkin_cart;
extern CINTIntegralFunctionReal    cint1e_ipkin_sph;
extern CINTIntegralFunctionComplex cint1e_ipkin_spinor;

extern CINTOptimizerFunction       cint1e_ipkinip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipkinip_cart;
extern CINTIntegralFunctionReal    cint1e_ipkinip_sph;
extern CINTIntegralFunctionComplex cint1e_ipkinip_spinor;

extern CINTOptimizerFunction       cint1e_ipnuc_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipnuc_cart;
extern CINTIntegralFunctionReal    cint1e_ipnuc_sph;
extern CINTIntegralFunctionComplex cint1e_ipnuc_spinor;

extern CINTOptimizerFunction       cint1e_ipnucip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipnucip_cart;
extern CINTIntegralFunctionReal    cint1e_ipnucip_sph;
extern CINTIntegralFunctionComplex cint1e_ipnucip_spinor;

extern CINTOptimizerFunction       cint1e_ipovlp_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipovlp_cart;
extern CINTIntegralFunctionReal    cint1e_ipovlp_sph;
extern CINTIntegralFunctionComplex cint1e_ipovlp_spinor;

extern CINTOptimizerFunction       cint1e_ipovlpip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipovlpip_cart;
extern CINTIntegralFunctionReal    cint1e_ipovlpip_sph;
extern CINTIntegralFunctionComplex cint1e_ipovlpip_spinor;

extern CINTOptimizerFunction       cint1e_ippnucp_optimizer;
extern CINTIntegralFunctionReal    cint1e_ippnucp_cart;
extern CINTIntegralFunctionReal    cint1e_ippnucp_sph;
extern CINTIntegralFunctionComplex cint1e_ippnucp_spinor;

extern CINTOptimizerFunction       cint1e_ippnucpip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ippnucpip_cart;
extern CINTIntegralFunctionReal    cint1e_ippnucpip_sph;
extern CINTIntegralFunctionComplex cint1e_ippnucpip_spinor;

extern CINTOptimizerFunction       cint1e_ipprinvp_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipprinvp_cart;
extern CINTIntegralFunctionReal    cint1e_ipprinvp_sph;
extern CINTIntegralFunctionComplex cint1e_ipprinvp_spinor;

extern CINTOptimizerFunction       cint1e_ipprinvpip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipprinvpip_cart;
extern CINTIntegralFunctionReal    cint1e_ipprinvpip_sph;
extern CINTIntegralFunctionComplex cint1e_ipprinvpip_spinor;

extern CINTOptimizerFunction       cint1e_iprinv_optimizer;
extern CINTIntegralFunctionReal    cint1e_iprinv_cart;
extern CINTIntegralFunctionReal    cint1e_iprinv_sph;
extern CINTIntegralFunctionComplex cint1e_iprinv_spinor;

extern CINTOptimizerFunction       cint1e_iprinvip_optimizer;
extern CINTIntegralFunctionReal    cint1e_iprinvip_cart;
extern CINTIntegralFunctionReal    cint1e_iprinvip_sph;
extern CINTIntegralFunctionComplex cint1e_iprinvip_spinor;

extern CINTOptimizerFunction       cint1e_iprinviprip_optimizer;
extern CINTIntegralFunctionReal    cint1e_iprinviprip_cart;
extern CINTIntegralFunctionReal    cint1e_iprinviprip_sph;
extern CINTIntegralFunctionComplex cint1e_iprinviprip_spinor;

extern CINTOptimizerFunction       cint1e_iprinvr_optimizer;
extern CINTIntegralFunctionReal    cint1e_iprinvr_cart;
extern CINTIntegralFunctionReal    cint1e_iprinvr_sph;
extern CINTIntegralFunctionComplex cint1e_iprinvr_spinor;

extern CINTOptimizerFunction       cint1e_iprip_optimizer;
extern CINTIntegralFunctionReal    cint1e_iprip_cart;
extern CINTIntegralFunctionReal    cint1e_iprip_sph;
extern CINTIntegralFunctionComplex cint1e_iprip_spinor;

extern CINTOptimizerFunction       cint1e_ipspnucsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipspnucsp_cart;
extern CINTIntegralFunctionReal    cint1e_ipspnucsp_sph;
extern CINTIntegralFunctionComplex cint1e_ipspnucsp_spinor;

extern CINTOptimizerFunction       cint1e_ipspnucspip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipspnucspip_cart;
extern CINTIntegralFunctionReal    cint1e_ipspnucspip_sph;
extern CINTIntegralFunctionComplex cint1e_ipspnucspip_spinor;

extern CINTOptimizerFunction       cint1e_ipsprinvsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipsprinvsp_cart;
extern CINTIntegralFunctionReal    cint1e_ipsprinvsp_sph;
extern CINTIntegralFunctionComplex cint1e_ipsprinvsp_spinor;

extern CINTOptimizerFunction       cint1e_ipsprinvspip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ipsprinvspip_cart;
extern CINTIntegralFunctionReal    cint1e_ipsprinvspip_sph;
extern CINTIntegralFunctionComplex cint1e_ipsprinvspip_spinor;

extern CINTOptimizerFunction       cint1e_irp_optimizer;
extern CINTIntegralFunctionReal    cint1e_irp_cart;
extern CINTIntegralFunctionReal    cint1e_irp_sph;
extern CINTIntegralFunctionComplex cint1e_irp_spinor;

extern CINTOptimizerFunction       cint1e_irpr_optimizer;
extern CINTIntegralFunctionReal    cint1e_irpr_cart;
extern CINTIntegralFunctionReal    cint1e_irpr_sph;
extern CINTIntegralFunctionComplex cint1e_irpr_spinor;

extern CINTOptimizerFunction       cint1e_irrp_optimizer;
extern CINTIntegralFunctionReal    cint1e_irrp_cart;
extern CINTIntegralFunctionReal    cint1e_irrp_sph;
extern CINTIntegralFunctionComplex cint1e_irrp_spinor;

extern CINTOptimizerFunction       cint1e_kin_optimizer;
extern CINTIntegralFunctionReal    cint1e_kin_cart;
extern CINTIntegralFunctionReal    cint1e_kin_sph;
extern CINTIntegralFunctionComplex cint1e_kin_spinor;

extern CINTOptimizerFunction       cint1e_kinip_optimizer;
extern CINTIntegralFunctionReal    cint1e_kinip_cart;
extern CINTIntegralFunctionReal    cint1e_kinip_sph;
extern CINTIntegralFunctionComplex cint1e_kinip_spinor;

extern CINTOptimizerFunction       cint1e_nuc_optimizer;
extern CINTIntegralFunctionReal    cint1e_nuc_cart;
extern CINTIntegralFunctionReal    cint1e_nuc_sph;
extern CINTIntegralFunctionComplex cint1e_nuc_spinor;

extern CINTOptimizerFunction       cint1e_ovlp_optimizer;
extern CINTIntegralFunctionReal    cint1e_ovlp_cart;
extern CINTIntegralFunctionReal    cint1e_ovlp_sph;
extern CINTIntegralFunctionComplex cint1e_ovlp_spinor;

extern CINTOptimizerFunction       cint1e_ovlpip_optimizer;
extern CINTIntegralFunctionReal    cint1e_ovlpip_cart;
extern CINTIntegralFunctionReal    cint1e_ovlpip_sph;
extern CINTIntegralFunctionComplex cint1e_ovlpip_spinor;

extern CINTOptimizerFunction       cint1e_p4_optimizer;
extern CINTIntegralFunctionReal    cint1e_p4_cart;
extern CINTIntegralFunctionReal    cint1e_p4_sph;
extern CINTIntegralFunctionComplex cint1e_p4_spinor;

extern CINTOptimizerFunction       cint1e_pnucp_optimizer;
extern CINTIntegralFunctionReal    cint1e_pnucp_cart;
extern CINTIntegralFunctionReal    cint1e_pnucp_sph;
extern CINTIntegralFunctionComplex cint1e_pnucp_spinor;

extern CINTOptimizerFunction       cint1e_pnucxp_optimizer;
extern CINTIntegralFunctionReal    cint1e_pnucxp_cart;
extern CINTIntegralFunctionReal    cint1e_pnucxp_sph;
extern CINTIntegralFunctionComplex cint1e_pnucxp_spinor;

extern CINTOptimizerFunction       cint1e_prinvp_optimizer;
extern CINTIntegralFunctionReal    cint1e_prinvp_cart;
extern CINTIntegralFunctionReal    cint1e_prinvp_sph;
extern CINTIntegralFunctionComplex cint1e_prinvp_spinor;

extern CINTOptimizerFunction       cint1e_prinvxp_optimizer;
extern CINTIntegralFunctionReal    cint1e_prinvxp_cart;
extern CINTIntegralFunctionReal    cint1e_prinvxp_sph;
extern CINTIntegralFunctionComplex cint1e_prinvxp_spinor;

extern CINTOptimizerFunction       cint1e_r_optimizer;
extern CINTIntegralFunctionReal    cint1e_r_cart;
extern CINTIntegralFunctionReal    cint1e_r_sph;
extern CINTIntegralFunctionComplex cint1e_r_spinor;

extern CINTOptimizerFunction       cint1e_r2_optimizer;
extern CINTIntegralFunctionReal    cint1e_r2_cart;
extern CINTIntegralFunctionReal    cint1e_r2_sph;
extern CINTIntegralFunctionComplex cint1e_r2_spinor;

extern CINTOptimizerFunction       cint1e_r2_origi_optimizer;
extern CINTIntegralFunctionReal    cint1e_r2_origi_cart;
extern CINTIntegralFunctionReal    cint1e_r2_origi_sph;
extern CINTIntegralFunctionComplex cint1e_r2_origi_spinor;

extern CINTOptimizerFunction       cint1e_r2_origi_ip2_optimizer;
extern CINTIntegralFunctionReal    cint1e_r2_origi_ip2_cart;
extern CINTIntegralFunctionReal    cint1e_r2_origi_ip2_sph;
extern CINTIntegralFunctionComplex cint1e_r2_origi_ip2_spinor;

extern CINTOptimizerFunction       cint1e_r2_origj_optimizer;
extern CINTIntegralFunctionReal    cint1e_r2_origj_cart;
extern CINTIntegralFunctionReal    cint1e_r2_origj_sph;
extern CINTIntegralFunctionComplex cint1e_r2_origj_spinor;

extern CINTOptimizerFunction       cint1e_r4_optimizer;
extern CINTIntegralFunctionReal    cint1e_r4_cart;
extern CINTIntegralFunctionReal    cint1e_r4_sph;
extern CINTIntegralFunctionComplex cint1e_r4_spinor;

extern CINTOptimizerFunction       cint1e_r4_origi_optimizer;
extern CINTIntegralFunctionReal    cint1e_r4_origi_cart;
extern CINTIntegralFunctionReal    cint1e_r4_origi_sph;
extern CINTIntegralFunctionComplex cint1e_r4_origi_spinor;

extern CINTOptimizerFunction       cint1e_r4_origi_ip2_optimizer;
extern CINTIntegralFunctionReal    cint1e_r4_origi_ip2_cart;
extern CINTIntegralFunctionReal    cint1e_r4_origi_ip2_sph;
extern CINTIntegralFunctionComplex cint1e_r4_origi_ip2_spinor;

extern CINTOptimizerFunction       cint1e_r4_origj_optimizer;
extern CINTIntegralFunctionReal    cint1e_r4_origj_cart;
extern CINTIntegralFunctionReal    cint1e_r4_origj_sph;
extern CINTIntegralFunctionComplex cint1e_r4_origj_spinor;

extern CINTOptimizerFunction       cint1e_r_origj_optimizer;
extern CINTIntegralFunctionReal    cint1e_r_origj_cart;
extern CINTIntegralFunctionReal    cint1e_r_origj_sph;
extern CINTIntegralFunctionComplex cint1e_r_origj_spinor;

extern CINTOptimizerFunction       cint1e_rinv_optimizer;
extern CINTIntegralFunctionReal    cint1e_rinv_cart;
extern CINTIntegralFunctionReal    cint1e_rinv_sph;
extern CINTIntegralFunctionComplex cint1e_rinv_spinor;

extern CINTOptimizerFunction       cint1e_rinvipiprip_optimizer;
extern CINTIntegralFunctionReal    cint1e_rinvipiprip_cart;
extern CINTIntegralFunctionReal    cint1e_rinvipiprip_sph;
extern CINTIntegralFunctionComplex cint1e_rinvipiprip_spinor;

extern CINTOptimizerFunction       cint1e_rr_optimizer;
extern CINTIntegralFunctionReal    cint1e_rr_cart;
extern CINTIntegralFunctionReal    cint1e_rr_sph;
extern CINTIntegralFunctionComplex cint1e_rr_spinor;

extern CINTOptimizerFunction       cint1e_rr_origj_optimizer;
extern CINTIntegralFunctionReal    cint1e_rr_origj_cart;
extern CINTIntegralFunctionReal    cint1e_rr_origj_sph;
extern CINTIntegralFunctionComplex cint1e_rr_origj_spinor;

extern CINTOptimizerFunction       cint1e_rrr_optimizer;
extern CINTIntegralFunctionReal    cint1e_rrr_cart;
extern CINTIntegralFunctionReal    cint1e_rrr_sph;
extern CINTIntegralFunctionComplex cint1e_rrr_spinor;

extern CINTOptimizerFunction       cint1e_rrrr_optimizer;
extern CINTIntegralFunctionReal    cint1e_rrrr_cart;
extern CINTIntegralFunctionReal    cint1e_rrrr_sph;
extern CINTIntegralFunctionComplex cint1e_rrrr_spinor;

extern CINTOptimizerFunction       cint1e_sa01sp_optimizer;
extern CINTIntegralFunctionReal    cint1e_sa01sp_cart;
extern CINTIntegralFunctionReal    cint1e_sa01sp_sph;
extern CINTIntegralFunctionComplex cint1e_sa01sp_spinor;

extern CINTOptimizerFunction       cint1e_sigma_optimizer;
extern CINTIntegralFunctionReal    cint1e_sigma_cart;
extern CINTIntegralFunctionReal    cint1e_sigma_sph;
extern CINTIntegralFunctionComplex cint1e_sigma_spinor;

extern CINTOptimizerFunction       cint1e_sp_optimizer;
extern CINTIntegralFunctionReal    cint1e_sp_cart;
extern CINTIntegralFunctionReal    cint1e_sp_sph;
extern CINTIntegralFunctionComplex cint1e_sp_spinor;

extern CINTOptimizerFunction       cint1e_spgnucsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_spgnucsp_cart;
extern CINTIntegralFunctionReal    cint1e_spgnucsp_sph;
extern CINTIntegralFunctionComplex cint1e_spgnucsp_spinor;

extern CINTOptimizerFunction       cint1e_spgsa01_optimizer;
extern CINTIntegralFunctionReal    cint1e_spgsa01_cart;
extern CINTIntegralFunctionReal    cint1e_spgsa01_sph;
extern CINTIntegralFunctionComplex cint1e_spgsa01_spinor;

extern CINTOptimizerFunction       cint1e_spgsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_spgsp_cart;
extern CINTIntegralFunctionReal    cint1e_spgsp_sph;
extern CINTIntegralFunctionComplex cint1e_spgsp_spinor;

extern CINTOptimizerFunction       cint1e_spnuc_optimizer;
extern CINTIntegralFunctionReal    cint1e_spnuc_cart;
extern CINTIntegralFunctionReal    cint1e_spnuc_sph;
extern CINTIntegralFunctionComplex cint1e_spnuc_spinor;

extern CINTOptimizerFunction       cint1e_spnucsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_spnucsp_cart;
extern CINTIntegralFunctionReal    cint1e_spnucsp_sph;
extern CINTIntegralFunctionComplex cint1e_spnucsp_spinor;

extern CINTOptimizerFunction       cint1e_sprinvsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_sprinvsp_cart;
extern CINTIntegralFunctionReal    cint1e_sprinvsp_sph;
extern CINTIntegralFunctionComplex cint1e_sprinvsp_spinor;

extern CINTOptimizerFunction       cint1e_sprsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_sprsp_cart;
extern CINTIntegralFunctionReal    cint1e_sprsp_sph;
extern CINTIntegralFunctionComplex cint1e_sprsp_spinor;

extern CINTOptimizerFunction       cint1e_spsigmasp_optimizer;
extern CINTIntegralFunctionReal    cint1e_spsigmasp_cart;
extern CINTIntegralFunctionReal    cint1e_spsigmasp_sph;
extern CINTIntegralFunctionComplex cint1e_spsigmasp_spinor;

extern CINTOptimizerFunction       cint1e_spsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_spsp_cart;
extern CINTIntegralFunctionReal    cint1e_spsp_sph;
extern CINTIntegralFunctionComplex cint1e_spsp_spinor;

extern CINTOptimizerFunction       cint1e_spspsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_spspsp_cart;
extern CINTIntegralFunctionReal    cint1e_spspsp_sph;
extern CINTIntegralFunctionComplex cint1e_spspsp_spinor;

extern CINTOptimizerFunction       cint1e_sr_optimizer;
extern CINTIntegralFunctionReal    cint1e_sr_cart;
extern CINTIntegralFunctionReal    cint1e_sr_sph;
extern CINTIntegralFunctionComplex cint1e_sr_spinor;

extern CINTOptimizerFunction       cint1e_srnucsr_optimizer;
extern CINTIntegralFunctionReal    cint1e_srnucsr_cart;
extern CINTIntegralFunctionReal    cint1e_srnucsr_sph;
extern CINTIntegralFunctionComplex cint1e_srnucsr_spinor;

extern CINTOptimizerFunction       cint1e_srsp_optimizer;
extern CINTIntegralFunctionReal    cint1e_srsp_cart;
extern CINTIntegralFunctionReal    cint1e_srsp_sph;
extern CINTIntegralFunctionComplex cint1e_srsp_spinor;

extern CINTOptimizerFunction       cint1e_srsr_optimizer;
extern CINTIntegralFunctionReal    cint1e_srsr_cart;
extern CINTIntegralFunctionReal    cint1e_srsr_sph;
extern CINTIntegralFunctionComplex cint1e_srsr_spinor;

extern CINTOptimizerFunction       cint1e_z_optimizer;
extern CINTIntegralFunctionReal    cint1e_z_cart;
extern CINTIntegralFunctionReal    cint1e_z_sph;
extern CINTIntegralFunctionComplex cint1e_z_spinor;

extern CINTOptimizerFunction       cint1e_z_origj_optimizer;
extern CINTIntegralFunctionReal    cint1e_z_origj_cart;
extern CINTIntegralFunctionReal    cint1e_z_origj_sph;
extern CINTIntegralFunctionComplex cint1e_z_origj_spinor;

extern CINTOptimizerFunction       cint1e_zz_optimizer;
extern CINTIntegralFunctionReal    cint1e_zz_cart;
extern CINTIntegralFunctionReal    cint1e_zz_sph;
extern CINTIntegralFunctionComplex cint1e_zz_spinor;

extern CINTOptimizerFunction       cint1e_zz_origj_optimizer;
extern CINTIntegralFunctionReal    cint1e_zz_origj_cart;
extern CINTIntegralFunctionReal    cint1e_zz_origj_sph;
extern CINTIntegralFunctionComplex cint1e_zz_origj_spinor;

extern CINTOptimizerFunction       cint2c2e_optimizer;
extern CINTIntegralFunctionReal    cint2c2e_cart;
extern CINTIntegralFunctionReal    cint2c2e_sph;
extern CINTIntegralFunctionComplex cint2c2e_spinor;

extern CINTOptimizerFunction       cint2c2e_ip1_optimizer;
extern CINTIntegralFunctionReal    cint2c2e_ip1_cart;
extern CINTIntegralFunctionReal    cint2c2e_ip1_sph;
extern CINTIntegralFunctionComplex cint2c2e_ip1_spinor;

extern CINTOptimizerFunction       cint2c2e_ip1ip2_optimizer;
extern CINTIntegralFunctionReal    cint2c2e_ip1ip2_cart;
extern CINTIntegralFunctionReal    cint2c2e_ip1ip2_sph;
extern CINTIntegralFunctionComplex cint2c2e_ip1ip2_spinor;

extern CINTOptimizerFunction       cint2c2e_ip2_optimizer;
extern CINTIntegralFunctionReal    cint2c2e_ip2_cart;
extern CINTIntegralFunctionReal    cint2c2e_ip2_sph;
extern CINTIntegralFunctionComplex cint2c2e_ip2_spinor;

extern CINTOptimizerFunction       cint2c2e_ipip1_optimizer;
extern CINTIntegralFunctionReal    cint2c2e_ipip1_cart;
extern CINTIntegralFunctionReal    cint2c2e_ipip1_sph;
extern CINTIntegralFunctionComplex cint2c2e_ipip1_spinor;

extern CINTOptimizerFunction       cint2e_breit_r1p2_optimizer;
extern CINTIntegralFunctionReal    cint2e_breit_r1p2_cart;
extern CINTIntegralFunctionReal    cint2e_breit_r1p2_sph;
extern CINTIntegralFunctionComplex cint2e_breit_r1p2_spinor;

extern CINTOptimizerFunction       cint2e_breit_r2p2_optimizer;
extern CINTIntegralFunctionReal    cint2e_breit_r2p2_cart;
extern CINTIntegralFunctionReal    cint2e_breit_r2p2_sph;
extern CINTIntegralFunctionComplex cint2e_breit_r2p2_spinor;

extern CINTOptimizerFunction       cint2e_cg_sa10sp1_optimizer;
extern CINTIntegralFunctionReal    cint2e_cg_sa10sp1_cart;
extern CINTIntegralFunctionReal    cint2e_cg_sa10sp1_sph;
extern CINTIntegralFunctionComplex cint2e_cg_sa10sp1_spinor;

extern CINTOptimizerFunction       cint2e_cg_sa10sp1spsp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_cg_sa10sp1spsp2_cart;
extern CINTIntegralFunctionReal    cint2e_cg_sa10sp1spsp2_sph;
extern CINTIntegralFunctionComplex cint2e_cg_sa10sp1spsp2_spinor;

extern CINTOptimizerFunction       cint2e_cg_ssa10ssp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_cg_ssa10ssp2_cart;
extern CINTIntegralFunctionReal    cint2e_cg_ssa10ssp2_sph;
extern CINTIntegralFunctionComplex cint2e_cg_ssa10ssp2_spinor;

extern CINTOptimizerFunction       cint2e_g1_optimizer;
extern CINTIntegralFunctionReal    cint2e_g1_cart;
extern CINTIntegralFunctionReal    cint2e_g1_sph;
extern CINTIntegralFunctionComplex cint2e_g1_spinor;

extern CINTOptimizerFunction       cint2e_g1g2_optimizer;
extern CINTIntegralFunctionReal    cint2e_g1g2_cart;
extern CINTIntegralFunctionReal    cint2e_g1g2_sph;
extern CINTIntegralFunctionComplex cint2e_g1g2_spinor;

extern CINTOptimizerFunction       cint2e_g1spsp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_g1spsp2_cart;
extern CINTIntegralFunctionReal    cint2e_g1spsp2_sph;
extern CINTIntegralFunctionComplex cint2e_g1spsp2_spinor;

extern CINTOptimizerFunction       cint2e_gauge_r1_sps1sps2_optimizer;
extern CINTIntegralFunctionReal    cint2e_gauge_r1_sps1sps2_cart;
extern CINTIntegralFunctionReal    cint2e_gauge_r1_sps1sps2_sph;
extern CINTIntegralFunctionComplex cint2e_gauge_r1_sps1sps2_spinor;

extern CINTOptimizerFunction       cint2e_gauge_r1_sps1ssp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_gauge_r1_sps1ssp2_cart;
extern CINTIntegralFunctionReal    cint2e_gauge_r1_sps1ssp2_sph;
extern CINTIntegralFunctionComplex cint2e_gauge_r1_sps1ssp2_spinor;

extern CINTOptimizerFunction       cint2e_gauge_r1_ssp1sps2_optimizer;
extern CINTIntegralFunctionReal    cint2e_gauge_r1_ssp1sps2_cart;
extern CINTIntegralFunctionReal    cint2e_gauge_r1_ssp1sps2_sph;
extern CINTIntegralFunctionComplex cint2e_gauge_r1_ssp1sps2_spinor;

extern CINTOptimizerFunction       cint2e_gauge_r1_ssp1ssp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_gauge_r1_ssp1ssp2_cart;
extern CINTIntegralFunctionReal    cint2e_gauge_r1_ssp1ssp2_sph;
extern CINTIntegralFunctionComplex cint2e_gauge_r1_ssp1ssp2_spinor;

extern CINTOptimizerFunction       cint2e_gauge_r2_sps1sps2_optimizer;
extern CINTIntegralFunctionReal    cint2e_gauge_r2_sps1sps2_cart;
extern CINTIntegralFunctionReal    cint2e_gauge_r2_sps1sps2_sph;
extern CINTIntegralFunctionComplex cint2e_gauge_r2_sps1sps2_spinor;

extern CINTOptimizerFunction       cint2e_gauge_r2_sps1ssp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_gauge_r2_sps1ssp2_cart;
extern CINTIntegralFunctionReal    cint2e_gauge_r2_sps1ssp2_sph;
extern CINTIntegralFunctionComplex cint2e_gauge_r2_sps1ssp2_spinor;

extern CINTOptimizerFunction       cint2e_gauge_r2_ssp1sps2_optimizer;
extern CINTIntegralFunctionReal    cint2e_gauge_r2_ssp1sps2_cart;
extern CINTIntegralFunctionReal    cint2e_gauge_r2_ssp1sps2_sph;
extern CINTIntegralFunctionComplex cint2e_gauge_r2_ssp1sps2_spinor;

extern CINTOptimizerFunction       cint2e_gauge_r2_ssp1ssp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_gauge_r2_ssp1ssp2_cart;
extern CINTIntegralFunctionReal    cint2e_gauge_r2_ssp1ssp2_sph;
extern CINTIntegralFunctionComplex cint2e_gauge_r2_ssp1ssp2_spinor;

extern CINTOptimizerFunction       cint2e_gg1_optimizer;
extern CINTIntegralFunctionReal    cint2e_gg1_cart;
extern CINTIntegralFunctionReal    cint2e_gg1_sph;
extern CINTIntegralFunctionComplex cint2e_gg1_spinor;

extern CINTOptimizerFunction       cint2e_giao_sa10sp1_optimizer;
extern CINTIntegralFunctionReal    cint2e_giao_sa10sp1_cart;
extern CINTIntegralFunctionReal    cint2e_giao_sa10sp1_sph;
extern CINTIntegralFunctionComplex cint2e_giao_sa10sp1_spinor;

extern CINTOptimizerFunction       cint2e_giao_sa10sp1spsp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_giao_sa10sp1spsp2_cart;
extern CINTIntegralFunctionReal    cint2e_giao_sa10sp1spsp2_sph;
extern CINTIntegralFunctionComplex cint2e_giao_sa10sp1spsp2_spinor;

extern CINTOptimizerFunction       cint2e_giao_ssa10ssp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_giao_ssa10ssp2_cart;
extern CINTIntegralFunctionReal    cint2e_giao_ssa10ssp2_sph;
extern CINTIntegralFunctionComplex cint2e_giao_ssa10ssp2_spinor;

extern CINTOptimizerFunction       cint2e_gssp1ssp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_gssp1ssp2_cart;
extern CINTIntegralFunctionReal    cint2e_gssp1ssp2_sph;
extern CINTIntegralFunctionComplex cint2e_gssp1ssp2_spinor;

extern CINTOptimizerFunction       cint2e_ig1_optimizer;
extern CINTIntegralFunctionReal    cint2e_ig1_cart;
extern CINTIntegralFunctionReal    cint2e_ig1_sph;
extern CINTIntegralFunctionComplex cint2e_ig1_spinor;

extern CINTOptimizerFunction       cint2e_ip1_optimizer;
extern CINTIntegralFunctionReal    cint2e_ip1_cart;
extern CINTIntegralFunctionReal    cint2e_ip1_sph;
extern CINTIntegralFunctionComplex cint2e_ip1_spinor;

extern CINTOptimizerFunction       cint2e_ip1ip2_optimizer;
extern CINTIntegralFunctionReal    cint2e_ip1ip2_cart;
extern CINTIntegralFunctionReal    cint2e_ip1ip2_sph;
extern CINTIntegralFunctionComplex cint2e_ip1ip2_spinor;

extern CINTOptimizerFunction       cint2e_ip1spsp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_ip1spsp2_cart;
extern CINTIntegralFunctionReal    cint2e_ip1spsp2_sph;
extern CINTIntegralFunctionComplex cint2e_ip1spsp2_spinor;

extern CINTOptimizerFunction       cint2e_ip1srsr2_optimizer;
extern CINTIntegralFunctionReal    cint2e_ip1srsr2_cart;
extern CINTIntegralFunctionReal    cint2e_ip1srsr2_sph;
extern CINTIntegralFunctionComplex cint2e_ip1srsr2_spinor;

extern CINTOptimizerFunction       cint2e_ip1v_r1_optimizer;
extern CINTIntegralFunctionReal    cint2e_ip1v_r1_cart;
extern CINTIntegralFunctionReal    cint2e_ip1v_r1_sph;
extern CINTIntegralFunctionComplex cint2e_ip1v_r1_spinor;

extern CINTOptimizerFunction       cint2e_ip1v_rc1_optimizer;
extern CINTIntegralFunctionReal    cint2e_ip1v_rc1_cart;
extern CINTIntegralFunctionReal    cint2e_ip1v_rc1_sph;
extern CINTIntegralFunctionComplex cint2e_ip1v_rc1_spinor;

extern CINTOptimizerFunction       cint2e_ip2_optimizer;
extern CINTIntegralFunctionReal    cint2e_ip2_cart;
extern CINTIntegralFunctionReal    cint2e_ip2_sph;
extern CINTIntegralFunctionComplex cint2e_ip2_spinor;

extern CINTOptimizerFunction       cint2e_ipip1_optimizer;
extern CINTIntegralFunctionReal    cint2e_ipip1_cart;
extern CINTIntegralFunctionReal    cint2e_ipip1_sph;
extern CINTIntegralFunctionComplex cint2e_ipip1_spinor;

extern CINTOptimizerFunction       cint2e_ipip1ipip2_optimizer;
extern CINTIntegralFunctionReal    cint2e_ipip1ipip2_cart;
extern CINTIntegralFunctionReal    cint2e_ipip1ipip2_sph;
extern CINTIntegralFunctionComplex cint2e_ipip1ipip2_spinor;

extern CINTOptimizerFunction       cint2e_ipspsp1_optimizer;
extern CINTIntegralFunctionReal    cint2e_ipspsp1_cart;
extern CINTIntegralFunctionReal    cint2e_ipspsp1_sph;
extern CINTIntegralFunctionComplex cint2e_ipspsp1_spinor;

extern CINTOptimizerFunction       cint2e_ipspsp1spsp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_ipspsp1spsp2_cart;
extern CINTIntegralFunctionReal    cint2e_ipspsp1spsp2_sph;
extern CINTIntegralFunctionComplex cint2e_ipspsp1spsp2_spinor;

extern CINTOptimizerFunction       cint2e_ipsrsr1_optimizer;
extern CINTIntegralFunctionReal    cint2e_ipsrsr1_cart;
extern CINTIntegralFunctionReal    cint2e_ipsrsr1_sph;
extern CINTIntegralFunctionComplex cint2e_ipsrsr1_spinor;

extern CINTOptimizerFunction       cint2e_ipsrsr1srsr2_optimizer;
extern CINTIntegralFunctionReal    cint2e_ipsrsr1srsr2_cart;
extern CINTIntegralFunctionReal    cint2e_ipsrsr1srsr2_sph;
extern CINTIntegralFunctionComplex cint2e_ipsrsr1srsr2_spinor;

extern CINTOptimizerFunction       cint2e_ipvg1_xp1_optimizer;
extern CINTIntegralFunctionReal    cint2e_ipvg1_xp1_cart;
extern CINTIntegralFunctionReal    cint2e_ipvg1_xp1_sph;
extern CINTIntegralFunctionComplex cint2e_ipvg1_xp1_spinor;

extern CINTOptimizerFunction       cint2e_ipvg2_xp1_optimizer;
extern CINTIntegralFunctionReal    cint2e_ipvg2_xp1_cart;
extern CINTIntegralFunctionReal    cint2e_ipvg2_xp1_sph;
extern CINTIntegralFunctionComplex cint2e_ipvg2_xp1_spinor;

extern CINTOptimizerFunction       cint2e_ipvip1_optimizer;
extern CINTIntegralFunctionReal    cint2e_ipvip1_cart;
extern CINTIntegralFunctionReal    cint2e_ipvip1_sph;
extern CINTIntegralFunctionComplex cint2e_ipvip1_spinor;

extern CINTOptimizerFunction       cint2e_ipvip1ipvip2_optimizer;
extern CINTIntegralFunctionReal    cint2e_ipvip1ipvip2_cart;
extern CINTIntegralFunctionReal    cint2e_ipvip1ipvip2_sph;
extern CINTIntegralFunctionComplex cint2e_ipvip1ipvip2_spinor;

extern CINTOptimizerFunction       cint2e_p1vxp1_optimizer;
extern CINTIntegralFunctionReal    cint2e_p1vxp1_cart;
extern CINTIntegralFunctionReal    cint2e_p1vxp1_sph;
extern CINTIntegralFunctionComplex cint2e_p1vxp1_spinor;

extern CINTOptimizerFunction       cint2e_pp1_optimizer;
extern CINTIntegralFunctionReal    cint2e_pp1_cart;
extern CINTIntegralFunctionReal    cint2e_pp1_sph;
extern CINTIntegralFunctionComplex cint2e_pp1_spinor;

extern CINTOptimizerFunction       cint2e_pp1pp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_pp1pp2_cart;
extern CINTIntegralFunctionReal    cint2e_pp1pp2_sph;
extern CINTIntegralFunctionComplex cint2e_pp1pp2_spinor;

extern CINTOptimizerFunction       cint2e_pp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_pp2_cart;
extern CINTIntegralFunctionReal    cint2e_pp2_sph;
extern CINTIntegralFunctionComplex cint2e_pp2_spinor;

extern CINTOptimizerFunction       cint2e_spgsp1_optimizer;
extern CINTIntegralFunctionReal    cint2e_spgsp1_cart;
extern CINTIntegralFunctionReal    cint2e_spgsp1_sph;
extern CINTIntegralFunctionComplex cint2e_spgsp1_spinor;

extern CINTOptimizerFunction       cint2e_spgsp1spsp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_spgsp1spsp2_cart;
extern CINTIntegralFunctionReal    cint2e_spgsp1spsp2_sph;
extern CINTIntegralFunctionComplex cint2e_spgsp1spsp2_spinor;

extern CINTOptimizerFunction       cint2e_sps1sps2_optimizer;
extern CINTIntegralFunctionReal    cint2e_sps1sps2_cart;
extern CINTIntegralFunctionReal    cint2e_sps1sps2_sph;
extern CINTIntegralFunctionComplex cint2e_sps1sps2_spinor;

extern CINTOptimizerFunction       cint2e_sps1ssp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_sps1ssp2_cart;
extern CINTIntegralFunctionReal    cint2e_sps1ssp2_sph;
extern CINTIntegralFunctionComplex cint2e_sps1ssp2_spinor;

extern CINTOptimizerFunction       cint2e_spsp1_optimizer;
extern CINTIntegralFunctionReal    cint2e_spsp1_cart;
extern CINTIntegralFunctionReal    cint2e_spsp1_sph;
extern CINTIntegralFunctionComplex cint2e_spsp1_spinor;

extern CINTOptimizerFunction       cint2e_spsp1spsp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_spsp1spsp2_cart;
extern CINTIntegralFunctionReal    cint2e_spsp1spsp2_sph;
extern CINTIntegralFunctionComplex cint2e_spsp1spsp2_spinor;

extern CINTOptimizerFunction       cint2e_spsp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_spsp2_cart;
extern CINTIntegralFunctionReal    cint2e_spsp2_sph;
extern CINTIntegralFunctionComplex cint2e_spsp2_spinor;

extern CINTOptimizerFunction       cint2e_spv1_optimizer;
extern CINTIntegralFunctionReal    cint2e_spv1_cart;
extern CINTIntegralFunctionReal    cint2e_spv1_sph;
extern CINTIntegralFunctionComplex cint2e_spv1_spinor;

extern CINTOptimizerFunction       cint2e_spv1spsp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_spv1spsp2_cart;
extern CINTIntegralFunctionReal    cint2e_spv1spsp2_sph;
extern CINTIntegralFunctionComplex cint2e_spv1spsp2_spinor;

extern CINTOptimizerFunction       cint2e_spv1spv2_optimizer;
extern CINTIntegralFunctionReal    cint2e_spv1spv2_cart;
extern CINTIntegralFunctionReal    cint2e_spv1spv2_sph;
extern CINTIntegralFunctionComplex cint2e_spv1spv2_spinor;

extern CINTOptimizerFunction       cint2e_spv1vsp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_spv1vsp2_cart;
extern CINTIntegralFunctionReal    cint2e_spv1vsp2_sph;
extern CINTIntegralFunctionComplex cint2e_spv1vsp2_spinor;

extern CINTOptimizerFunction       cint2e_srsr1_optimizer;
extern CINTIntegralFunctionReal    cint2e_srsr1_cart;
extern CINTIntegralFunctionReal    cint2e_srsr1_sph;
extern CINTIntegralFunctionComplex cint2e_srsr1_spinor;

extern CINTOptimizerFunction       cint2e_srsr1srsr2_optimizer;
extern CINTIntegralFunctionReal    cint2e_srsr1srsr2_cart;
extern CINTIntegralFunctionReal    cint2e_srsr1srsr2_sph;
extern CINTIntegralFunctionComplex cint2e_srsr1srsr2_spinor;

extern CINTOptimizerFunction       cint2e_ssp1sps2_optimizer;
extern CINTIntegralFunctionReal    cint2e_ssp1sps2_cart;
extern CINTIntegralFunctionReal    cint2e_ssp1sps2_sph;
extern CINTIntegralFunctionComplex cint2e_ssp1sps2_spinor;

extern CINTOptimizerFunction       cint2e_ssp1ssp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_ssp1ssp2_cart;
extern CINTIntegralFunctionReal    cint2e_ssp1ssp2_sph;
extern CINTIntegralFunctionComplex cint2e_ssp1ssp2_spinor;

extern CINTOptimizerFunction       cint2e_stg_optimizer;
extern CINTIntegralFunctionReal    cint2e_stg_cart;
extern CINTIntegralFunctionReal    cint2e_stg_sph;
extern CINTIntegralFunctionComplex cint2e_stg_spinor;

extern CINTOptimizerFunction       cint2e_stg_ip1_optimizer;
extern CINTIntegralFunctionReal    cint2e_stg_ip1_cart;
extern CINTIntegralFunctionReal    cint2e_stg_ip1_sph;
extern CINTIntegralFunctionComplex cint2e_stg_ip1_spinor;

extern CINTOptimizerFunction       cint2e_stg_ip1ip2_optimizer;
extern CINTIntegralFunctionReal    cint2e_stg_ip1ip2_cart;
extern CINTIntegralFunctionReal    cint2e_stg_ip1ip2_sph;
extern CINTIntegralFunctionComplex cint2e_stg_ip1ip2_spinor;

extern CINTOptimizerFunction       cint2e_stg_ipip1_optimizer;
extern CINTIntegralFunctionReal    cint2e_stg_ipip1_cart;
extern CINTIntegralFunctionReal    cint2e_stg_ipip1_sph;
extern CINTIntegralFunctionComplex cint2e_stg_ipip1_spinor;

extern CINTOptimizerFunction       cint2e_stg_ipvip1_optimizer;
extern CINTIntegralFunctionReal    cint2e_stg_ipvip1_cart;
extern CINTIntegralFunctionReal    cint2e_stg_ipvip1_sph;
extern CINTIntegralFunctionComplex cint2e_stg_ipvip1_spinor;

extern CINTOptimizerFunction       cint2e_vsp1_optimizer;
extern CINTIntegralFunctionReal    cint2e_vsp1_cart;
extern CINTIntegralFunctionReal    cint2e_vsp1_sph;
extern CINTIntegralFunctionComplex cint2e_vsp1_spinor;

extern CINTOptimizerFunction       cint2e_vsp1spsp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_vsp1spsp2_cart;
extern CINTIntegralFunctionReal    cint2e_vsp1spsp2_sph;
extern CINTIntegralFunctionComplex cint2e_vsp1spsp2_spinor;

extern CINTOptimizerFunction       cint2e_vsp1spv2_optimizer;
extern CINTIntegralFunctionReal    cint2e_vsp1spv2_cart;
extern CINTIntegralFunctionReal    cint2e_vsp1spv2_sph;
extern CINTIntegralFunctionComplex cint2e_vsp1spv2_spinor;

extern CINTOptimizerFunction       cint2e_vsp1vsp2_optimizer;
extern CINTIntegralFunctionReal    cint2e_vsp1vsp2_cart;
extern CINTIntegralFunctionReal    cint2e_vsp1vsp2_sph;
extern CINTIntegralFunctionComplex cint2e_vsp1vsp2_spinor;

extern CINTOptimizerFunction       cint2e_yp_optimizer;
extern CINTIntegralFunctionReal    cint2e_yp_cart;
extern CINTIntegralFunctionReal    cint2e_yp_sph;
extern CINTIntegralFunctionComplex cint2e_yp_spinor;

extern CINTOptimizerFunction       cint2e_yp_ip1_optimizer;
extern CINTIntegralFunctionReal    cint2e_yp_ip1_cart;
extern CINTIntegralFunctionReal    cint2e_yp_ip1_sph;
extern CINTIntegralFunctionComplex cint2e_yp_ip1_spinor;

extern CINTOptimizerFunction       cint2e_yp_ip1ip2_optimizer;
extern CINTIntegralFunctionReal    cint2e_yp_ip1ip2_cart;
extern CINTIntegralFunctionReal    cint2e_yp_ip1ip2_sph;
extern CINTIntegralFunctionComplex cint2e_yp_ip1ip2_spinor;

extern CINTOptimizerFunction       cint2e_yp_ipip1_optimizer;
extern CINTIntegralFunctionReal    cint2e_yp_ipip1_cart;
extern CINTIntegralFunctionReal    cint2e_yp_ipip1_sph;
extern CINTIntegralFunctionComplex cint2e_yp_ipip1_spinor;

extern CINTOptimizerFunction       cint2e_yp_ipvip1_optimizer;
extern CINTIntegralFunctionReal    cint2e_yp_ipvip1_cart;
extern CINTIntegralFunctionReal    cint2e_yp_ipvip1_sph;
extern CINTIntegralFunctionComplex cint2e_yp_ipvip1_spinor;

extern CINTOptimizerFunction       cint3c1e_optimizer;
extern CINTIntegralFunctionReal    cint3c1e_cart;
extern CINTIntegralFunctionReal    cint3c1e_sph;
extern CINTIntegralFunctionComplex cint3c1e_spinor;

extern CINTOptimizerFunction       cint3c1e_ip1_optimizer;
extern CINTIntegralFunctionReal    cint3c1e_ip1_cart;
extern CINTIntegralFunctionReal    cint3c1e_ip1_sph;
extern CINTIntegralFunctionComplex cint3c1e_ip1_spinor;

extern CINTOptimizerFunction       cint3c1e_ip1_r2_origk_optimizer;
extern CINTIntegralFunctionReal    cint3c1e_ip1_r2_origk_cart;
extern CINTIntegralFunctionReal    cint3c1e_ip1_r2_origk_sph;
extern CINTIntegralFunctionComplex cint3c1e_ip1_r2_origk_spinor;

extern CINTOptimizerFunction       cint3c1e_ip1_r4_origk_optimizer;
extern CINTIntegralFunctionReal    cint3c1e_ip1_r4_origk_cart;
extern CINTIntegralFunctionReal    cint3c1e_ip1_r4_origk_sph;
extern CINTIntegralFunctionComplex cint3c1e_ip1_r4_origk_spinor;

extern CINTOptimizerFunction       cint3c1e_ip1_r6_origk_optimizer;
extern CINTIntegralFunctionReal    cint3c1e_ip1_r6_origk_cart;
extern CINTIntegralFunctionReal    cint3c1e_ip1_r6_origk_sph;
extern CINTIntegralFunctionComplex cint3c1e_ip1_r6_origk_spinor;

extern CINTOptimizerFunction       cint3c1e_iprinv_optimizer;
extern CINTIntegralFunctionReal    cint3c1e_iprinv_cart;
extern CINTIntegralFunctionReal    cint3c1e_iprinv_sph;
extern CINTIntegralFunctionComplex cint3c1e_iprinv_spinor;

extern CINTOptimizerFunction       cint3c1e_p2_optimizer;
extern CINTIntegralFunctionReal    cint3c1e_p2_cart;
extern CINTIntegralFunctionReal    cint3c1e_p2_sph;
extern CINTIntegralFunctionComplex cint3c1e_p2_spinor;

extern CINTOptimizerFunction       cint3c1e_r2_origk_optimizer;
extern CINTIntegralFunctionReal    cint3c1e_r2_origk_cart;
extern CINTIntegralFunctionReal    cint3c1e_r2_origk_sph;
extern CINTIntegralFunctionComplex cint3c1e_r2_origk_spinor;

extern CINTOptimizerFunction       cint3c1e_r4_origk_optimizer;
extern CINTIntegralFunctionReal    cint3c1e_r4_origk_cart;
extern CINTIntegralFunctionReal    cint3c1e_r4_origk_sph;
extern CINTIntegralFunctionComplex cint3c1e_r4_origk_spinor;

extern CINTOptimizerFunction       cint3c1e_r6_origk_optimizer;
extern CINTIntegralFunctionReal    cint3c1e_r6_origk_cart;
extern CINTIntegralFunctionReal    cint3c1e_r6_origk_sph;
extern CINTIntegralFunctionComplex cint3c1e_r6_origk_spinor;

extern CINTOptimizerFunction       cint3c1e_rinv_optimizer;
extern CINTIntegralFunctionReal    cint3c1e_rinv_cart;
extern CINTIntegralFunctionReal    cint3c1e_rinv_sph;
extern CINTIntegralFunctionComplex cint3c1e_rinv_spinor;

extern CINTOptimizerFunction       cint3c2e_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_cart;
extern CINTIntegralFunctionReal    cint3c2e_sph;
extern CINTIntegralFunctionComplex cint3c2e_spinor;

extern CINTOptimizerFunction       cint3c2e_ig1_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_ig1_cart;
extern CINTIntegralFunctionReal    cint3c2e_ig1_sph;
extern CINTIntegralFunctionComplex cint3c2e_ig1_spinor;

extern CINTOptimizerFunction       cint3c2e_ip1_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_ip1_cart;
extern CINTIntegralFunctionReal    cint3c2e_ip1_sph;
extern CINTIntegralFunctionComplex cint3c2e_ip1_spinor;

extern CINTOptimizerFunction       cint3c2e_ip1ip2_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_ip1ip2_cart;
extern CINTIntegralFunctionReal    cint3c2e_ip1ip2_sph;
extern CINTIntegralFunctionComplex cint3c2e_ip1ip2_spinor;

extern CINTOptimizerFunction       cint3c2e_ip2_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_ip2_cart;
extern CINTIntegralFunctionReal    cint3c2e_ip2_sph;
extern CINTIntegralFunctionComplex cint3c2e_ip2_spinor;

extern CINTOptimizerFunction       cint3c2e_ipip1_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_ipip1_cart;
extern CINTIntegralFunctionReal    cint3c2e_ipip1_sph;
extern CINTIntegralFunctionComplex cint3c2e_ipip1_spinor;

extern CINTOptimizerFunction       cint3c2e_ipip2_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_ipip2_cart;
extern CINTIntegralFunctionReal    cint3c2e_ipip2_sph;
extern CINTIntegralFunctionComplex cint3c2e_ipip2_spinor;

extern CINTOptimizerFunction       cint3c2e_ipspsp1_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_ipspsp1_cart;
extern CINTIntegralFunctionReal    cint3c2e_ipspsp1_sph;
extern CINTIntegralFunctionComplex cint3c2e_ipspsp1_spinor;

extern CINTOptimizerFunction       cint3c2e_ipvip1_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_ipvip1_cart;
extern CINTIntegralFunctionReal    cint3c2e_ipvip1_sph;
extern CINTIntegralFunctionComplex cint3c2e_ipvip1_spinor;

extern CINTOptimizerFunction       cint3c2e_pvp1_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_pvp1_cart;
extern CINTIntegralFunctionReal    cint3c2e_pvp1_sph;
extern CINTIntegralFunctionComplex cint3c2e_pvp1_spinor;

extern CINTOptimizerFunction       cint3c2e_pvxp1_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_pvxp1_cart;
extern CINTIntegralFunctionReal    cint3c2e_pvxp1_sph;
extern CINTIntegralFunctionComplex cint3c2e_pvxp1_spinor;

extern CINTOptimizerFunction       cint3c2e_spsp1_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_spsp1_cart;
extern CINTIntegralFunctionReal    cint3c2e_spsp1_sph;
extern CINTIntegralFunctionComplex cint3c2e_spsp1_spinor;

extern CINTOptimizerFunction       cint3c2e_spsp1ip2_optimizer;
extern CINTIntegralFunctionReal    cint3c2e_spsp1ip2_cart;
extern CINTIntegralFunctionReal    cint3c2e_spsp1ip2_sph;
extern CINTIntegralFunctionComplex cint3c2e_spsp1ip2_spinor;

extern CINTOptimizerFunction       cint4c1e_optimizer;
extern CINTIntegralFunctionReal    cint4c1e_cart;
extern CINTIntegralFunctionReal    cint4c1e_sph;
extern CINTIntegralFunctionComplex cint4c1e_spinor;
