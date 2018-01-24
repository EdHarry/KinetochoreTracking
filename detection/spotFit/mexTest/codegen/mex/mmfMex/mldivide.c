/*
 * mldivide.c
 *
 * Code generation for function 'mldivide'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "mldivide.h"
#include "mmfMex_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "eml_warning.h"
#include "colon.h"
#include "mmfMex_mexutil.h"
#include "mmfMex_data.h"

/* Variable Definitions */
static emlrtRSInfo tb_emlrtRSI = { 19, "sqrt",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elfun/sqrt.m" };

static emlrtRSInfo ed_emlrtRSI = { 19, "abs",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elfun/abs.m" };

static emlrtRSInfo ge_emlrtRSI = { 1, "mldivide",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mldivide.p" };

static emlrtRSInfo he_emlrtRSI = { 20, "eml_lusolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_lusolve.m" };

static emlrtRSInfo ie_emlrtRSI = { 70, "eml_lusolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_lusolve.m" };

static emlrtRSInfo je_emlrtRSI = { 68, "eml_lusolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_lusolve.m" };

static emlrtRSInfo ke_emlrtRSI = { 77, "eml_lusolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_lusolve.m" };

static emlrtRSInfo le_emlrtRSI = { 88, "eml_lusolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_lusolve.m" };

static emlrtRSInfo me_emlrtRSI = { 90, "eml_lusolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_lusolve.m" };

static emlrtRSInfo ne_emlrtRSI = { 16, "min",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/min.m" };

static emlrtRSInfo oe_emlrtRSI = { 18, "eml_min_or_max",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"
};

static emlrtRSInfo pe_emlrtRSI = { 59, "eml_min_or_max",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"
};

static emlrtRSInfo qe_emlrtRSI = { 124, "eml_min_or_max",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_min_or_max.m"
};

static emlrtRSInfo re_emlrtRSI = { 8, "eml_xgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/eml_xgetrf.m"
};

static emlrtRSInfo se_emlrtRSI = { 8, "eml_lapack_xgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgetrf.m"
};

static emlrtRSInfo te_emlrtRSI = { 58, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo ue_emlrtRSI = { 60, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo ve_emlrtRSI = { 59, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo we_emlrtRSI = { 51, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo xe_emlrtRSI = { 50, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo ye_emlrtRSI = { 44, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo af_emlrtRSI = { 43, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo bf_emlrtRSI = { 42, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo cf_emlrtRSI = { 41, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo df_emlrtRSI = { 37, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo ef_emlrtRSI = { 36, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo ff_emlrtRSI = { 34, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo gf_emlrtRSI = { 33, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo hf_emlrtRSI = { 32, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo if_emlrtRSI = { 31, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo jf_emlrtRSI = { 30, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo kf_emlrtRSI = { 29, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo lf_emlrtRSI = { 23, "eml_matlab_zgetrf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgetrf.m"
};

static emlrtRSInfo nf_emlrtRSI = { 75, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo of_emlrtRSI = { 106, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo pf_emlrtRSI = { 112, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo yf_emlrtRSI = { 20, "eml_ixamax",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/eml_ixamax.m"
};

static emlrtRSInfo ag_emlrtRSI = { 15, "eml_blas_ixamax",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m"
};

static emlrtRSInfo cg_emlrtRSI = { 24, "eml_blas_ixamax",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m"
};

static emlrtRSInfo hg_emlrtRSI = { 15, "eml_xcabs1",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/eml_xcabs1.m"
};

static emlrtRSInfo ig_emlrtRSI = { 69, "eml_blas_ixamax",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m"
};

static emlrtRSInfo jg_emlrtRSI = { 70, "eml_blas_ixamax",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_ixamax.m"
};

static emlrtRSInfo kg_emlrtRSI = { 26, "eml_xswap",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/eml_xswap.m"
};

static emlrtRSInfo lg_emlrtRSI = { 15, "eml_blas_xswap",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xswap.m"
};

static emlrtRSInfo mg_emlrtRSI = { 36, "eml_refblas_xswap",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"
};

static emlrtRSInfo ng_emlrtRSI = { 31, "eml_refblas_xswap",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"
};

static emlrtRSInfo og_emlrtRSI = { 19, "eml_refblas_xswap",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"
};

static emlrtRSInfo pg_emlrtRSI = { 18, "eml_refblas_xswap",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"
};

static emlrtRSInfo qg_emlrtRSI = { 17, "eml_refblas_xswap",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/refblas/eml_refblas_xswap.m"
};

static emlrtRSInfo kh_emlrtRSI = { 54, "eml_lusolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_lusolve.m" };

static emlrtRSInfo lh_emlrtRSI = { 54, "eml_xtrsm",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/eml_xtrsm.m"
};

static emlrtRSInfo mh_emlrtRSI = { 28, "eml_blas_xtrsm",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m"
};

static emlrtRSInfo oh_emlrtRSI = { 17, "eml_blas_xtrsm",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m"
};

static emlrtRSInfo ai_emlrtRSI = { 87, "eml_blas_xtrsm",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m"
};

static emlrtRSInfo bi_emlrtRSI = { 88, "eml_blas_xtrsm",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m"
};

static emlrtRSInfo ci_emlrtRSI = { 89, "eml_blas_xtrsm",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m"
};

static emlrtRSInfo di_emlrtRSI = { 90, "eml_blas_xtrsm",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xtrsm.m"
};

static emlrtRSInfo pi_emlrtRSI = { 28, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRSInfo qi_emlrtRSI = { 29, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRSInfo ri_emlrtRSI = { 32, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRSInfo si_emlrtRSI = { 34, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRSInfo ti_emlrtRSI = { 38, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRSInfo ui_emlrtRSI = { 37, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRSInfo vi_emlrtRSI = { 77, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRSInfo wi_emlrtRSI = { 79, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRSInfo xi_emlrtRSI = { 83, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRSInfo yi_emlrtRSI = { 108, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRSInfo aj_emlrtRSI = { 110, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRSInfo bj_emlrtRSI = { 8, "eml_xgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/eml_xgeqp3.m"
};

static emlrtRSInfo cj_emlrtRSI = { 8, "eml_lapack_xgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/internal/eml_lapack_xgeqp3.m"
};

static emlrtRSInfo dj_emlrtRSI = { 15, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo ej_emlrtRSI = { 19, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo fj_emlrtRSI = { 31, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo gj_emlrtRSI = { 32, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo hj_emlrtRSI = { 34, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo ij_emlrtRSI = { 37, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo jj_emlrtRSI = { 38, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo kj_emlrtRSI = { 39, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo lj_emlrtRSI = { 40, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo mj_emlrtRSI = { 41, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo nj_emlrtRSI = { 42, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo oj_emlrtRSI = { 43, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo pj_emlrtRSI = { 44, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo qj_emlrtRSI = { 47, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo rj_emlrtRSI = { 49, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo sj_emlrtRSI = { 50, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo tj_emlrtRSI = { 51, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo uj_emlrtRSI = { 64, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo vj_emlrtRSI = { 66, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo wj_emlrtRSI = { 73, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo xj_emlrtRSI = { 74, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo yj_emlrtRSI = { 79, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo ak_emlrtRSI = { 80, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo bk_emlrtRSI = { 84, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo ck_emlrtRSI = { 85, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo dk_emlrtRSI = { 90, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo ek_emlrtRSI = { 93, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo fk_emlrtRSI = { 100, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRSInfo gk_emlrtRSI = { 19, "eml_xnrm2",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/eml_xnrm2.m"
};

static emlrtRSInfo hk_emlrtRSI = { 17, "eml_blas_xnrm2",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xnrm2.m"
};

static emlrtRSInfo jk_emlrtRSI = { 24, "eml_blas_xnrm2",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xnrm2.m"
};

static emlrtRSInfo ok_emlrtRSI = { 65, "eml_blas_xnrm2",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xnrm2.m"
};

static emlrtRSInfo pk_emlrtRSI = { 66, "eml_blas_xnrm2",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/external/eml_blas_xnrm2.m"
};

static emlrtRSInfo qk_emlrtRSI = { 18, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo rk_emlrtRSI = { 22, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo sk_emlrtRSI = { 31, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo tk_emlrtRSI = { 35, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo uk_emlrtRSI = { 39, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo vk_emlrtRSI = { 41, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo wk_emlrtRSI = { 42, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo xk_emlrtRSI = { 43, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo yk_emlrtRSI = { 51, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo al_emlrtRSI = { 54, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo bl_emlrtRSI = { 61, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo cl_emlrtRSI = { 62, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo dl_emlrtRSI = { 66, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo el_emlrtRSI = { 69, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo fl_emlrtRSI = { 70, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo gl_emlrtRSI = { 74, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo hl_emlrtRSI = { 75, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo il_emlrtRSI = { 79, "eml_matlab_zlarfg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarfg.m"
};

static emlrtRSInfo rl_emlrtRSI = { 101, "eml_matlab_zlarf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"
};

static emlrtRSInfo sl_emlrtRSI = { 102, "eml_matlab_zlarf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"
};

static emlrtRSInfo tl_emlrtRSI = { 103, "eml_matlab_zlarf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"
};

static emlrtRSInfo ul_emlrtRSI = { 108, "eml_matlab_zlarf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"
};

static emlrtRSInfo vl_emlrtRSI = { 32, "eml_matlab_zlarf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"
};

static emlrtRSInfo wl_emlrtRSI = { 39, "eml_matlab_zlarf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"
};

static emlrtRSInfo xl_emlrtRSI = { 40, "eml_matlab_zlarf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"
};

static emlrtRSInfo yl_emlrtRSI = { 50, "eml_matlab_zlarf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"
};

static emlrtRSInfo am_emlrtRSI = { 68, "eml_matlab_zlarf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"
};

static emlrtRSInfo bm_emlrtRSI = { 75, "eml_matlab_zlarf",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zlarf.m"
};

static emlrtRSInfo cm_emlrtRSI = { 52, "eml_xgemv",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/eml_xgemv.m"
};

static emlrtRSInfo um_emlrtRSI = { 42, "eml_xgerc",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/blas/eml_xgerc.m"
};

static emlrtRSInfo vm_emlrtRSI = { 16, "max",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/max.m" };

static emlrtMCInfo o_emlrtMCI = { 1, 1, "mldivide",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mldivide.p" };

static emlrtMCInfo t_emlrtMCI = { 29, 23, "eml_flt2str",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_flt2str.m" };

static emlrtMCInfo u_emlrtMCI = { 29, 15, "eml_flt2str",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_flt2str.m" };

static emlrtRTEInfo q_emlrtRTEI = { 1, 2, "mldivide",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mldivide.p" };

static emlrtRTEInfo r_emlrtRTEI = { 1, 19, "eml_lusolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_lusolve.m" };

static emlrtRTEInfo u_emlrtRTEI = { 1, 24, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRTEInfo v_emlrtRTEI = { 16, 1, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRTEInfo w_emlrtRTEI = { 28, 5, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRTEInfo x_emlrtRTEI = { 29, 5, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRTEInfo y_emlrtRTEI = { 24, 1, "eml_matlab_zgeqp3",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/lapack/matlab/eml_matlab_zgeqp3.m"
};

static emlrtRTEInfo nb_emlrtRTEI = { 76, 17, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRTEInfo ob_emlrtRTEI = { 82, 21, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRTEInfo pb_emlrtRTEI = { 99, 5, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRTEInfo qb_emlrtRTEI = { 106, 5, "eml_qrsolve",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_qrsolve.m" };

static emlrtRSInfo qq_emlrtRSI = { 29, "eml_flt2str",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_flt2str.m" };

/* Function Declarations */
static real_T b_eml_matlab_zlarfg(void);
static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[14]);
static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, const mxArray *d, emlrtMCInfo *location);
static const mxArray *c_sprintf(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, emlrtMCInfo *location);
static int32_T eml_ixamax(int32_T n, const emxArray_real_T *x, int32_T ix0);
static void eml_lusolve(const emlrtStack *sp, const emxArray_real_T *A,
  emxArray_real_T *B);
static void eml_matlab_zlarf(const emlrtStack *sp, int32_T m, int32_T n, int32_T
  iv0, real_T tau, emxArray_real_T *C, int32_T ic0, int32_T ldc, emxArray_real_T
  *work);
static real_T eml_matlab_zlarfg(const emlrtStack *sp, int32_T n, real_T *alpha1,
  emxArray_real_T *x, int32_T ix0);
static void eml_qrsolve(const emlrtStack *sp, const emxArray_real_T *A,
  emxArray_real_T *B, emxArray_real_T *Y);
static void eml_xgeru(int32_T m, int32_T n, int32_T ix0, int32_T iy0, int32_T
                      incy, emxArray_real_T *A, int32_T ia0, int32_T lda);
static real_T eml_xnrm2(int32_T n, const emxArray_real_T *x, int32_T ix0);
static void eml_xscal(int32_T n, real_T a, emxArray_real_T *x, int32_T ix0);
static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *d_sprintf,
  const char_T *identifier, char_T y[14]);
static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[14]);
static void warn_singular(const emlrtStack *sp);

/* Function Definitions */
static real_T b_eml_matlab_zlarfg(void)
{
  return 0.0;
}

static void b_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, char_T y[14])
{
  i_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static const mxArray *b_sprintf(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, const mxArray *d, emlrtMCInfo *location)
{
  const mxArray *pArrays[3];
  const mxArray *m20;
  pArrays[0] = b;
  pArrays[1] = c;
  pArrays[2] = d;
  return emlrtCallMATLABR2012b(sp, 1, &m20, 3, pArrays, "sprintf", TRUE,
    location);
}

static const mxArray *c_sprintf(const emlrtStack *sp, const mxArray *b, const
  mxArray *c, emlrtMCInfo *location)
{
  const mxArray *pArrays[2];
  const mxArray *m21;
  pArrays[0] = b;
  pArrays[1] = c;
  return emlrtCallMATLABR2012b(sp, 1, &m21, 2, pArrays, "sprintf", TRUE,
    location);
}

static int32_T eml_ixamax(int32_T n, const emxArray_real_T *x, int32_T ix0)
{
  int32_T idxmax;
  ptrdiff_t n_t;
  ptrdiff_t incx_t;
  double * xix0_t;
  if (n < 1) {
    idxmax = 0;
  } else {
    n_t = (ptrdiff_t)(n);
    incx_t = (ptrdiff_t)(1);
    xix0_t = (double *)(&x->data[ix0 - 1]);
    n_t = idamax(&n_t, xix0_t, &incx_t);
    idxmax = (int32_T)(n_t);
  }

  return idxmax;
}

static void eml_lusolve(const emlrtStack *sp, const emxArray_real_T *A,
  emxArray_real_T *B)
{
  emxArray_real_T *b_A;
  int32_T i;
  int32_T iy;
  emxArray_int32_T *ipiv;
  int32_T info;
  int32_T b;
  int32_T j;
  int32_T mmj;
  int32_T b_b;
  ptrdiff_t n_t;
  ptrdiff_t m_t;
  double * Aia0_t;
  int32_T ix;
  boolean_T overflow;
  real_T temp;
  boolean_T c_b;
  char_T DIAGA;
  char_T TRANSA;
  char_T UPLO;
  char_T SIDE;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  double * Bib0_t;
  double * alpha1_t;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  emlrtStack i_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  h_st.prev = &g_st;
  h_st.tls = g_st.tls;
  i_st.prev = &h_st;
  i_st.tls = h_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &b_A, 2, &r_emlrtRTEI, TRUE);
  st.site = &he_emlrtRSI;
  b_st.site = &je_emlrtRSI;
  c_st.site = &re_emlrtRSI;
  d_st.site = &se_emlrtRSI;
  i = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity(&d_st, (emxArray__common *)b_A, i, (int32_T)sizeof(real_T),
                    &r_emlrtRTEI);
  iy = A->size[0] * A->size[1];
  for (i = 0; i < iy; i++) {
    b_A->data[i] = A->data[i];
  }

  b_emxInit_int32_T(&d_st, &ipiv, 2, &r_emlrtRTEI, TRUE);
  e_st.site = &lf_emlrtRSI;
  f_st.site = &ne_emlrtRSI;
  g_st.site = &oe_emlrtRSI;
  h_st.site = &pe_emlrtRSI;
  e_st.site = &lf_emlrtRSI;
  f_st.site = &mf_emlrtRSI;
  g_st.site = &nf_emlrtRSI;
  h_st.site = &of_emlrtRSI;
  h_st.site = &pf_emlrtRSI;
  eml_signed_integer_colon(&h_st, muIntScalarMin_sint32(A->size[1], A->size[1]),
    ipiv);
  info = 0;
  e_st.site = &kf_emlrtRSI;
  e_st.site = &jf_emlrtRSI;
  f_st.site = &ne_emlrtRSI;
  g_st.site = &oe_emlrtRSI;
  h_st.site = &pe_emlrtRSI;
  b = muIntScalarMin_sint32(A->size[1] - 1, A->size[1]);
  e_st.site = &jf_emlrtRSI;
  f_st.site = &gb_emlrtRSI;
  for (j = 0; j + 1 <= b; j++) {
    e_st.site = &if_emlrtRSI;
    e_st.site = &hf_emlrtRSI;
    mmj = A->size[1] - j;
    e_st.site = &gf_emlrtRSI;
    b_b = j * (A->size[1] + 1);
    e_st.site = &gf_emlrtRSI;
    e_st.site = &ff_emlrtRSI;
    e_st.site = &ef_emlrtRSI;
    e_st.site = &ef_emlrtRSI;
    f_st.site = &yf_emlrtRSI;
    g_st.site = &ag_emlrtRSI;
    if (mmj < 1) {
      iy = -1;
    } else {
      g_st.site = &cg_emlrtRSI;
      h_st.site = &ig_emlrtRSI;
      n_t = (ptrdiff_t)(mmj);
      h_st.site = &jg_emlrtRSI;
      m_t = (ptrdiff_t)(1);
      Aia0_t = (double *)(&b_A->data[b_b]);
      m_t = idamax(&n_t, Aia0_t, &m_t);
      iy = (int32_T)(m_t) - 1;
    }

    e_st.site = &ef_emlrtRSI;
    e_st.site = &df_emlrtRSI;
    if (b_A->data[b_b + iy] != 0.0) {
      if (iy != 0) {
        e_st.site = &cf_emlrtRSI;
        ipiv->data[j] = (j + iy) + 1;
        e_st.site = &bf_emlrtRSI;
        e_st.site = &af_emlrtRSI;
        e_st.site = &ye_emlrtRSI;
        f_st.site = &kg_emlrtRSI;
        g_st.site = &lg_emlrtRSI;
        ix = j;
        iy += j;
        h_st.site = &qg_emlrtRSI;
        h_st.site = &pg_emlrtRSI;
        h_st.site = &og_emlrtRSI;
        if (1 > A->size[1]) {
          overflow = FALSE;
        } else {
          overflow = (A->size[1] > 2147483646);
        }

        if (overflow) {
          i_st.site = &hb_emlrtRSI;
          check_forloop_overflow_error(&i_st);
        }

        for (i = 1; i <= A->size[1]; i++) {
          temp = b_A->data[ix];
          b_A->data[ix] = b_A->data[iy];
          b_A->data[iy] = temp;
          h_st.site = &ng_emlrtRSI;
          ix += A->size[1];
          h_st.site = &mg_emlrtRSI;
          iy += A->size[1];
        }
      }

      e_st.site = &xe_emlrtRSI;
      e_st.site = &xe_emlrtRSI;
      iy = b_b + mmj;
      e_st.site = &xe_emlrtRSI;
      f_st.site = &gb_emlrtRSI;
      if (b_b + 2 > iy) {
        c_b = FALSE;
      } else {
        c_b = (iy > 2147483646);
      }

      if (c_b) {
        f_st.site = &hb_emlrtRSI;
        check_forloop_overflow_error(&f_st);
      }

      for (i = b_b + 1; i + 1 <= iy; i++) {
        e_st.site = &we_emlrtRSI;
        b_A->data[i] /= b_A->data[b_b];
      }
    } else {
      info = j + 1;
    }

    e_st.site = &te_emlrtRSI;
    e_st.site = &ve_emlrtRSI;
    e_st.site = &ue_emlrtRSI;
    e_st.site = &te_emlrtRSI;
    eml_xgeru(mmj - 1, (A->size[1] - j) - 1, b_b + 2, (b_b + A->size[1]) + 1,
              A->size[1], b_A, (b_b + A->size[1]) + 2, A->size[1]);
  }

  if ((info == 0) && (!(b_A->data[(A->size[1] + b_A->size[0] * (A->size[1] - 1))
                        - 1] != 0.0))) {
    info = A->size[1];
  }

  if (info > 0) {
    b_st.site = &ie_emlrtRSI;
    warn_singular(&b_st);
  }

  b_st.site = &ke_emlrtRSI;
  c_st.site = &gb_emlrtRSI;
  if (1 > A->size[1]) {
    overflow = FALSE;
  } else {
    overflow = (A->size[1] > 2147483646);
  }

  if (overflow) {
    c_st.site = &hb_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  for (i = 0; i + 1 <= A->size[1]; i++) {
    if (ipiv->data[i] != i + 1) {
      temp = B->data[i];
      B->data[i] = B->data[ipiv->data[i] - 1];
      B->data[ipiv->data[i] - 1] = temp;
    }
  }

  emxFree_int32_T(&ipiv);
  b_st.site = &le_emlrtRSI;
  c_st.site = &lh_emlrtRSI;
  d_st.site = &oh_emlrtRSI;
  d_st.site = &mh_emlrtRSI;
  temp = 1.0;
  DIAGA = 'U';
  TRANSA = 'N';
  UPLO = 'L';
  SIDE = 'L';
  e_st.site = &ai_emlrtRSI;
  m_t = (ptrdiff_t)(A->size[1]);
  e_st.site = &bi_emlrtRSI;
  n_t = (ptrdiff_t)(1);
  e_st.site = &ci_emlrtRSI;
  lda_t = (ptrdiff_t)(A->size[1]);
  e_st.site = &di_emlrtRSI;
  ldb_t = (ptrdiff_t)(A->size[1]);
  Aia0_t = (double *)(&b_A->data[0]);
  Bib0_t = (double *)(&B->data[0]);
  alpha1_t = (double *)(&temp);
  dtrsm(&SIDE, &UPLO, &TRANSA, &DIAGA, &m_t, &n_t, alpha1_t, Aia0_t, &lda_t,
        Bib0_t, &ldb_t);
  b_st.site = &me_emlrtRSI;
  c_st.site = &lh_emlrtRSI;
  d_st.site = &oh_emlrtRSI;
  d_st.site = &mh_emlrtRSI;
  temp = 1.0;
  DIAGA = 'N';
  TRANSA = 'N';
  UPLO = 'U';
  SIDE = 'L';
  e_st.site = &ai_emlrtRSI;
  m_t = (ptrdiff_t)(A->size[1]);
  e_st.site = &bi_emlrtRSI;
  n_t = (ptrdiff_t)(1);
  e_st.site = &ci_emlrtRSI;
  lda_t = (ptrdiff_t)(A->size[1]);
  e_st.site = &di_emlrtRSI;
  ldb_t = (ptrdiff_t)(A->size[1]);
  Aia0_t = (double *)(&b_A->data[0]);
  Bib0_t = (double *)(&B->data[0]);
  alpha1_t = (double *)(&temp);
  dtrsm(&SIDE, &UPLO, &TRANSA, &DIAGA, &m_t, &n_t, alpha1_t, Aia0_t, &lda_t,
        Bib0_t, &ldb_t);
  emxFree_real_T(&b_A);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

static void eml_matlab_zlarf(const emlrtStack *sp, int32_T m, int32_T n, int32_T
  iv0, real_T tau, emxArray_real_T *C, int32_T ic0, int32_T ldc, emxArray_real_T
  *work)
{
  int32_T lastv;
  int32_T lastc;
  boolean_T exitg2;
  int32_T coltop;
  int32_T colbottom;
  boolean_T b_coltop;
  int32_T exitg1;
  real_T alpha1;
  real_T beta1;
  char_T TRANSA;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t lda_t;
  ptrdiff_t incx_t;
  ptrdiff_t incy_t;
  double * alpha1_t;
  double * beta1_t;
  double * yiy0_t;
  double * Aia0_t;
  double * xix0_t;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  if (tau != 0.0) {
    lastv = m;
    st.site = &vl_emlrtRSI;
    st.site = &vl_emlrtRSI;
    st.site = &vl_emlrtRSI;
    lastc = iv0 + m;
    while ((lastv > 0) && (C->data[lastc - 2] == 0.0)) {
      st.site = &wl_emlrtRSI;
      lastv--;
      st.site = &xl_emlrtRSI;
      lastc--;
    }

    st.site = &yl_emlrtRSI;
    lastc = n;
    exitg2 = FALSE;
    while ((exitg2 == FALSE) && (lastc > 0)) {
      b_st.site = &rl_emlrtRSI;
      b_st.site = &rl_emlrtRSI;
      b_st.site = &rl_emlrtRSI;
      coltop = ic0 + (lastc - 1) * ldc;
      b_st.site = &sl_emlrtRSI;
      b_st.site = &sl_emlrtRSI;
      colbottom = (coltop + lastv) - 1;
      b_st.site = &tl_emlrtRSI;
      if (coltop > colbottom) {
        b_coltop = FALSE;
      } else {
        b_coltop = (colbottom > 2147483646);
      }

      if (b_coltop) {
        c_st.site = &hb_emlrtRSI;
        check_forloop_overflow_error(&c_st);
      }

      do {
        exitg1 = 0;
        if (coltop <= colbottom) {
          if (C->data[coltop - 1] != 0.0) {
            exitg1 = 1;
          } else {
            coltop++;
          }
        } else {
          b_st.site = &ul_emlrtRSI;
          lastc--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = TRUE;
      }
    }
  } else {
    lastv = 0;
    lastc = 0;
  }

  if (lastv > 0) {
    st.site = &am_emlrtRSI;
    b_st.site = &cm_emlrtRSI;
    if (lastc < 1) {
    } else {
      alpha1 = 1.0;
      beta1 = 0.0;
      TRANSA = 'C';
      m_t = (ptrdiff_t)(lastv);
      n_t = (ptrdiff_t)(lastc);
      lda_t = (ptrdiff_t)(ldc);
      incx_t = (ptrdiff_t)(1);
      incy_t = (ptrdiff_t)(1);
      alpha1_t = (double *)(&alpha1);
      beta1_t = (double *)(&beta1);
      yiy0_t = (double *)(&work->data[0]);
      Aia0_t = (double *)(&C->data[ic0 - 1]);
      xix0_t = (double *)(&C->data[iv0 - 1]);
      dgemv(&TRANSA, &m_t, &n_t, alpha1_t, Aia0_t, &lda_t, xix0_t, &incx_t,
            beta1_t, yiy0_t, &incy_t);
    }

    st.site = &bm_emlrtRSI;
    alpha1 = -tau;
    b_st.site = &um_emlrtRSI;
    if (lastc < 1) {
    } else {
      m_t = (ptrdiff_t)(lastv);
      n_t = (ptrdiff_t)(lastc);
      incx_t = (ptrdiff_t)(1);
      incy_t = (ptrdiff_t)(1);
      lda_t = (ptrdiff_t)(ldc);
      alpha1_t = (double *)(&alpha1);
      Aia0_t = (double *)(&C->data[ic0 - 1]);
      beta1_t = (double *)(&C->data[iv0 - 1]);
      yiy0_t = (double *)(&work->data[0]);
      dger(&m_t, &n_t, alpha1_t, beta1_t, &incx_t, yiy0_t, &incy_t, Aia0_t,
           &lda_t);
    }
  }
}

static real_T eml_matlab_zlarfg(const emlrtStack *sp, int32_T n, real_T *alpha1,
  emxArray_real_T *x, int32_T ix0)
{
  real_T tau;
  real_T xnorm;
  int32_T knt;
  real_T d0;
  boolean_T b9;
  int32_T k;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  tau = 0.0;
  if (n <= 0) {
  } else {
    st.site = &qk_emlrtRSI;
    xnorm = eml_xnrm2(n - 1, x, ix0);
    if (xnorm != 0.0) {
      st.site = &rk_emlrtRSI;
      xnorm = muDoubleScalarHypot(*alpha1, xnorm);
      if (*alpha1 >= 0.0) {
        xnorm = -xnorm;
      }

      st.site = &sk_emlrtRSI;
      if (muDoubleScalarAbs(xnorm) < 1.0020841800044864E-292) {
        knt = 0;
        do {
          st.site = &tk_emlrtRSI;
          knt++;
          st.site = &uk_emlrtRSI;
          eml_xscal(n - 1, 9.9792015476736E+291, x, ix0);
          st.site = &vk_emlrtRSI;
          xnorm *= 9.9792015476736E+291;
          st.site = &wk_emlrtRSI;
          *alpha1 *= 9.9792015476736E+291;
          st.site = &xk_emlrtRSI;
        } while (!(muDoubleScalarAbs(xnorm) >= 1.0020841800044864E-292));

        st.site = &yk_emlrtRSI;
        xnorm = eml_xnrm2(n - 1, x, ix0);
        st.site = &al_emlrtRSI;
        xnorm = muDoubleScalarHypot(*alpha1, xnorm);
        if (*alpha1 >= 0.0) {
          xnorm = -xnorm;
        }

        st.site = &bl_emlrtRSI;
        tau = (xnorm - *alpha1) / xnorm;
        st.site = &cl_emlrtRSI;
        d0 = 1.0 / (*alpha1 - xnorm);
        st.site = &dl_emlrtRSI;
        eml_xscal(n - 1, d0, x, ix0);
        st.site = &el_emlrtRSI;
        if (1 > knt) {
          b9 = FALSE;
        } else {
          b9 = (knt > 2147483646);
        }

        if (b9) {
          b_st.site = &hb_emlrtRSI;
          check_forloop_overflow_error(&b_st);
        }

        for (k = 1; k <= knt; k++) {
          st.site = &fl_emlrtRSI;
          xnorm *= 1.0020841800044864E-292;
        }

        *alpha1 = xnorm;
      } else {
        st.site = &gl_emlrtRSI;
        tau = (xnorm - *alpha1) / xnorm;
        st.site = &hl_emlrtRSI;
        d0 = 1.0 / (*alpha1 - xnorm);
        st.site = &il_emlrtRSI;
        eml_xscal(n - 1, d0, x, ix0);
        *alpha1 = xnorm;
      }
    }
  }

  return tau;
}

static void eml_qrsolve(const emlrtStack *sp, const emxArray_real_T *A,
  emxArray_real_T *B, emxArray_real_T *Y)
{
  emxArray_real_T *b_A;
  emxArray_real_T *work;
  int32_T mn;
  int32_T ix;
  int32_T iy;
  emxArray_real_T *tau;
  emxArray_int32_T *jpvt;
  int32_T m;
  int32_T n;
  int32_T b_mn;
  emxArray_real_T *vn1;
  emxArray_real_T *vn2;
  int32_T k;
  boolean_T overflow;
  boolean_T b7;
  int32_T i;
  int32_T i_i;
  int32_T nmi;
  int32_T mmi;
  int32_T pvt;
  boolean_T b8;
  real_T atmp;
  boolean_T b_i;
  real_T temp1;
  ptrdiff_t n_t;
  ptrdiff_t incx_t;
  double * xix0_t;
  boolean_T exitg1;
  const mxArray *y;
  static const int32_T iv19[2] = { 1, 8 };

  const mxArray *m8;
  char_T cv26[8];
  static const char_T cv27[8] = { '%', '%', '%', 'd', '.', '%', 'd', 'e' };

  char_T cv28[14];
  uint32_T unnamed_idx_0;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  emlrtStack g_st;
  emlrtStack h_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &e_st;
  f_st.tls = e_st.tls;
  g_st.prev = &f_st;
  g_st.tls = f_st.tls;
  h_st.prev = &g_st;
  h_st.tls = g_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &b_A, 2, &u_emlrtRTEI, TRUE);
  c_emxInit_real_T(sp, &work, 1, &y_emlrtRTEI, TRUE);
  st.site = &pi_emlrtRSI;
  b_st.site = &ne_emlrtRSI;
  c_st.site = &oe_emlrtRSI;
  d_st.site = &pe_emlrtRSI;
  e_st.site = &qe_emlrtRSI;
  mn = (int32_T)muDoubleScalarMin(A->size[0], A->size[1]);
  st.site = &qi_emlrtRSI;
  b_st.site = &bj_emlrtRSI;
  c_st.site = &cj_emlrtRSI;
  ix = b_A->size[0] * b_A->size[1];
  b_A->size[0] = A->size[0];
  b_A->size[1] = A->size[1];
  emxEnsureCapacity(&c_st, (emxArray__common *)b_A, ix, (int32_T)sizeof(real_T),
                    &u_emlrtRTEI);
  iy = A->size[0] * A->size[1];
  for (ix = 0; ix < iy; ix++) {
    b_A->data[ix] = A->data[ix];
  }

  c_emxInit_real_T(&c_st, &tau, 1, &u_emlrtRTEI, TRUE);
  b_emxInit_int32_T(&c_st, &jpvt, 2, &u_emlrtRTEI, TRUE);
  m = b_A->size[0];
  n = b_A->size[1];
  d_st.site = &dj_emlrtRSI;
  e_st.site = &ne_emlrtRSI;
  f_st.site = &oe_emlrtRSI;
  g_st.site = &pe_emlrtRSI;
  b_mn = muIntScalarMin_sint32(b_A->size[0], b_A->size[1]);
  ix = tau->size[0];
  tau->size[0] = b_mn;
  emxEnsureCapacity(&c_st, (emxArray__common *)tau, ix, (int32_T)sizeof(real_T),
                    &v_emlrtRTEI);
  d_st.site = &ej_emlrtRSI;
  e_st.site = &mf_emlrtRSI;
  f_st.site = &nf_emlrtRSI;
  g_st.site = &of_emlrtRSI;
  g_st.site = &pf_emlrtRSI;
  eml_signed_integer_colon(&g_st, b_A->size[1], jpvt);
  if (b_A->size[0] == 0) {
  } else {
    iy = b_A->size[1];
    ix = work->size[0];
    work->size[0] = iy;
    emxEnsureCapacity(&c_st, (emxArray__common *)work, ix, (int32_T)sizeof
                      (real_T), &u_emlrtRTEI);
    for (ix = 0; ix < iy; ix++) {
      work->data[ix] = 0.0;
    }

    c_emxInit_real_T(&c_st, &vn1, 1, &w_emlrtRTEI, TRUE);
    c_emxInit_real_T(&c_st, &vn2, 1, &x_emlrtRTEI, TRUE);
    iy = b_A->size[1];
    ix = vn1->size[0];
    vn1->size[0] = iy;
    emxEnsureCapacity(&c_st, (emxArray__common *)vn1, ix, (int32_T)sizeof(real_T),
                      &w_emlrtRTEI);
    ix = vn2->size[0];
    vn2->size[0] = iy;
    emxEnsureCapacity(&c_st, (emxArray__common *)vn2, ix, (int32_T)sizeof(real_T),
                      &x_emlrtRTEI);
    k = 1;
    d_st.site = &fj_emlrtRSI;
    e_st.site = &gb_emlrtRSI;
    if (1 > b_A->size[1]) {
      overflow = FALSE;
    } else {
      overflow = (b_A->size[1] > 2147483646);
    }

    if (overflow) {
      e_st.site = &hb_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }

    for (iy = 0; iy + 1 <= b_A->size[1]; iy++) {
      d_st.site = &gj_emlrtRSI;
      vn1->data[iy] = eml_xnrm2(b_A->size[0], b_A, k);
      vn2->data[iy] = vn1->data[iy];
      d_st.site = &hj_emlrtRSI;
      k += b_A->size[0];
    }

    d_st.site = &ij_emlrtRSI;
    e_st.site = &gb_emlrtRSI;
    if (1 > b_mn) {
      b7 = FALSE;
    } else {
      b7 = (b_mn > 2147483646);
    }

    if (b7) {
      e_st.site = &hb_emlrtRSI;
      check_forloop_overflow_error(&e_st);
    }

    for (i = 1; i <= b_mn; i++) {
      d_st.site = &jj_emlrtRSI;
      d_st.site = &kj_emlrtRSI;
      d_st.site = &lj_emlrtRSI;
      d_st.site = &lj_emlrtRSI;
      i_i = (i + (i - 1) * m) - 1;
      d_st.site = &mj_emlrtRSI;
      nmi = n - i;
      d_st.site = &nj_emlrtRSI;
      mmi = m - i;
      d_st.site = &oj_emlrtRSI;
      d_st.site = &pj_emlrtRSI;
      d_st.site = &qj_emlrtRSI;
      iy = eml_ixamax(1 + nmi, vn1, i);
      d_st.site = &qj_emlrtRSI;
      pvt = (i + iy) - 2;
      if (pvt + 1 != i) {
        d_st.site = &rj_emlrtRSI;
        d_st.site = &rj_emlrtRSI;
        d_st.site = &rj_emlrtRSI;
        d_st.site = &sj_emlrtRSI;
        d_st.site = &sj_emlrtRSI;
        d_st.site = &tj_emlrtRSI;
        e_st.site = &kg_emlrtRSI;
        f_st.site = &lg_emlrtRSI;
        ix = m * pvt;
        iy = m * (i - 1);
        g_st.site = &og_emlrtRSI;
        if (1 > m) {
          b8 = FALSE;
        } else {
          b8 = (m > 2147483646);
        }

        if (b8) {
          h_st.site = &hb_emlrtRSI;
          check_forloop_overflow_error(&h_st);
        }

        for (k = 1; k <= m; k++) {
          atmp = b_A->data[ix];
          b_A->data[ix] = b_A->data[iy];
          b_A->data[iy] = atmp;
          g_st.site = &ng_emlrtRSI;
          ix++;
          g_st.site = &mg_emlrtRSI;
          iy++;
        }

        iy = jpvt->data[pvt];
        jpvt->data[pvt] = jpvt->data[i - 1];
        jpvt->data[i - 1] = iy;
        vn1->data[pvt] = vn1->data[i - 1];
        vn2->data[pvt] = vn2->data[i - 1];
      }

      if (i < m) {
        d_st.site = &uj_emlrtRSI;
        atmp = b_A->data[i_i];
        d_st.site = &uj_emlrtRSI;
        tau->data[i - 1] = eml_matlab_zlarfg(&d_st, mmi + 1, &atmp, b_A, i_i + 2);
      } else {
        atmp = b_A->data[i_i];
        d_st.site = &vj_emlrtRSI;
        tau->data[i - 1] = b_eml_matlab_zlarfg();
      }

      b_A->data[i_i] = atmp;
      if (i < n) {
        atmp = b_A->data[i_i];
        b_A->data[i_i] = 1.0;
        d_st.site = &wj_emlrtRSI;
        d_st.site = &wj_emlrtRSI;
        d_st.site = &xj_emlrtRSI;
        eml_matlab_zlarf(&d_st, mmi + 1, nmi, i_i + 1, tau->data[i - 1], b_A, i
                         + i * m, m, work);
        b_A->data[i_i] = atmp;
      }

      d_st.site = &yj_emlrtRSI;
      e_st.site = &gb_emlrtRSI;
      if (i + 1 > n) {
        b_i = FALSE;
      } else {
        b_i = (n > 2147483646);
      }

      if (b_i) {
        e_st.site = &hb_emlrtRSI;
        check_forloop_overflow_error(&e_st);
      }

      for (iy = i; iy + 1 <= n; iy++) {
        d_st.site = &ak_emlrtRSI;
        d_st.site = &ak_emlrtRSI;
        d_st.site = &ak_emlrtRSI;
        if (vn1->data[iy] != 0.0) {
          d_st.site = &bk_emlrtRSI;
          e_st.site = &ed_emlrtRSI;
          temp1 = muDoubleScalarAbs(b_A->data[(i + b_A->size[0] * iy) - 1]) /
            vn1->data[iy];
          d_st.site = &ck_emlrtRSI;
          atmp = temp1 * temp1;
          temp1 = 1.0 - temp1 * temp1;
          if (1.0 - atmp < 0.0) {
            temp1 = 0.0;
          }

          atmp = vn1->data[iy] / vn2->data[iy];
          d_st.site = &dk_emlrtRSI;
          d_st.site = &dk_emlrtRSI;
          if (temp1 * (atmp * atmp) <= 1.4901161193847656E-8) {
            if (i < m) {
              d_st.site = &ek_emlrtRSI;
              d_st.site = &ek_emlrtRSI;
              e_st.site = &gk_emlrtRSI;
              f_st.site = &hk_emlrtRSI;
              if (mmi < 1) {
                atmp = 0.0;
              } else {
                f_st.site = &jk_emlrtRSI;
                g_st.site = &ok_emlrtRSI;
                n_t = (ptrdiff_t)(mmi);
                g_st.site = &pk_emlrtRSI;
                incx_t = (ptrdiff_t)(1);
                xix0_t = (double *)(&b_A->data[i + m * iy]);
                atmp = dnrm2(&n_t, xix0_t, &incx_t);
              }

              vn1->data[iy] = atmp;
              vn2->data[iy] = vn1->data[iy];
            } else {
              vn1->data[iy] = 0.0;
              vn2->data[iy] = 0.0;
            }
          } else {
            d_st.site = &fk_emlrtRSI;
            e_st.site = &tb_emlrtRSI;
            d_st.site = &fk_emlrtRSI;
            vn1->data[iy] *= muDoubleScalarSqrt(temp1);
          }
        }
      }
    }

    emxFree_real_T(&vn2);
    emxFree_real_T(&vn1);
  }

  temp1 = 0.0;
  if (mn > 0) {
    st.site = &ri_emlrtRSI;
    b_st.site = &vm_emlrtRSI;
    c_st.site = &oe_emlrtRSI;
    d_st.site = &pe_emlrtRSI;
    e_st.site = &qe_emlrtRSI;
    st.site = &ri_emlrtRSI;
    b_st.site = &hg_emlrtRSI;
    c_st.site = &ed_emlrtRSI;
    b_st.site = &hg_emlrtRSI;
    c_st.site = &ed_emlrtRSI;
    st.site = &ri_emlrtRSI;
    st.site = &ri_emlrtRSI;
    atmp = muDoubleScalarMax(A->size[0], A->size[1]) * muDoubleScalarAbs
      (b_A->data[0]) * 2.2204460492503131E-16;
    k = 0;
    exitg1 = FALSE;
    while ((exitg1 == FALSE) && (k <= mn - 1)) {
      st.site = &si_emlrtRSI;
      b_st.site = &hg_emlrtRSI;
      c_st.site = &ed_emlrtRSI;
      b_st.site = &hg_emlrtRSI;
      c_st.site = &ed_emlrtRSI;
      if (muDoubleScalarAbs(b_A->data[k + b_A->size[0] * k]) <= atmp) {
        st.site = &ti_emlrtRSI;
        y = NULL;
        m8 = mxCreateCharArray(2, iv19);
        for (i = 0; i < 8; i++) {
          cv26[i] = cv27[i];
        }

        emlrtInitCharArrayR2013a(&st, 8, m8, cv26);
        emlrtAssign(&y, m8);
        b_st.site = &qq_emlrtRSI;
        emlrt_marshallIn(&b_st, c_sprintf(&b_st, b_sprintf(&b_st, y,
          emlrt_marshallOut(14.0), emlrt_marshallOut(6.0), &t_emlrtMCI),
          emlrt_marshallOut(atmp), &u_emlrtMCI), "sprintf", cv28);
        st.site = &ui_emlrtRSI;
        b_eml_warning(&st, temp1, cv28);
        exitg1 = TRUE;
      } else {
        temp1++;
        k++;
      }
    }
  }

  unnamed_idx_0 = (uint32_T)A->size[1];
  ix = Y->size[0];
  Y->size[0] = (int32_T)unnamed_idx_0;
  emxEnsureCapacity(sp, (emxArray__common *)Y, ix, (int32_T)sizeof(real_T),
                    &u_emlrtRTEI);
  iy = (int32_T)unnamed_idx_0;
  for (ix = 0; ix < iy; ix++) {
    Y->data[ix] = 0.0;
  }

  for (iy = 0; iy < mn; iy++) {
    if (tau->data[iy] != 0.0) {
      atmp = B->data[iy];
      ix = A->size[0] + (int32_T)(1.0 - ((1.0 + (real_T)iy) + 1.0));
      emlrtForLoopVectorCheckR2012b((1.0 + (real_T)iy) + 1.0, 1.0, A->size[0],
        mxDOUBLE_CLASS, ix, &nb_emlrtRTEI, sp);
      for (i = 0; i < ix; i++) {
        unnamed_idx_0 = ((uint32_T)iy + i) + 2U;
        st.site = &vi_emlrtRSI;
        atmp += b_A->data[((int32_T)unnamed_idx_0 + b_A->size[0] * iy) - 1] *
          B->data[(int32_T)unnamed_idx_0 - 1];
      }

      st.site = &wi_emlrtRSI;
      atmp *= tau->data[iy];
      if (atmp != 0.0) {
        B->data[iy] -= atmp;
        ix = A->size[0] + (int32_T)(1.0 - ((1.0 + (real_T)iy) + 1.0));
        emlrtForLoopVectorCheckR2012b((1.0 + (real_T)iy) + 1.0, 1.0, A->size[0],
          mxDOUBLE_CLASS, ix, &ob_emlrtRTEI, sp);
        for (i = 0; i < ix; i++) {
          unnamed_idx_0 = ((uint32_T)iy + i) + 2U;
          st.site = &xi_emlrtRSI;
          B->data[(int32_T)unnamed_idx_0 - 1] -= b_A->data[((int32_T)
            unnamed_idx_0 + b_A->size[0] * iy) - 1] * atmp;
        }
      }
    }
  }

  emxFree_real_T(&tau);
  emlrtForLoopVectorCheckR2012b(1.0, 1.0, temp1, mxDOUBLE_CLASS, (int32_T)temp1,
    &pb_emlrtRTEI, sp);
  for (i = 0; i < (int32_T)temp1; i++) {
    Y->data[jpvt->data[i] - 1] = B->data[i];
  }

  emlrtForLoopVectorCheckR2012b(temp1, -1.0, 1.0, mxDOUBLE_CLASS, (int32_T)-(1.0
    + (-1.0 - temp1)), &qb_emlrtRTEI, sp);
  for (iy = 0; iy < (int32_T)-(1.0 + (-1.0 - temp1)); iy++) {
    atmp = temp1 + -(real_T)iy;
    st.site = &yi_emlrtRSI;
    Y->data[jpvt->data[(int32_T)atmp - 1] - 1] /= b_A->data[((int32_T)atmp +
      b_A->size[0] * ((int32_T)atmp - 1)) - 1];
    for (i = 0; i < (int32_T)(atmp - 1.0); i++) {
      st.site = &aj_emlrtRSI;
      Y->data[jpvt->data[i] - 1] -= Y->data[jpvt->data[(int32_T)atmp - 1] - 1] *
        b_A->data[i + b_A->size[0] * ((int32_T)atmp - 1)];
    }
  }

  emxFree_int32_T(&jpvt);
  emxFree_real_T(&work);
  emxFree_real_T(&b_A);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

static void eml_xgeru(int32_T m, int32_T n, int32_T ix0, int32_T iy0, int32_T
                      incy, emxArray_real_T *A, int32_T ia0, int32_T lda)
{
  real_T alpha1;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t incx_t;
  ptrdiff_t incy_t;
  ptrdiff_t lda_t;
  double * alpha1_t;
  double * Aia0_t;
  double * Aix0_t;
  double * Aiy0_t;
  if ((m < 1) || (n < 1)) {
  } else {
    alpha1 = -1.0;
    m_t = (ptrdiff_t)(m);
    n_t = (ptrdiff_t)(n);
    incx_t = (ptrdiff_t)(1);
    incy_t = (ptrdiff_t)(incy);
    lda_t = (ptrdiff_t)(lda);
    alpha1_t = (double *)(&alpha1);
    Aia0_t = (double *)(&A->data[ia0 - 1]);
    Aix0_t = (double *)(&A->data[ix0 - 1]);
    Aiy0_t = (double *)(&A->data[iy0 - 1]);
    dger(&m_t, &n_t, alpha1_t, Aix0_t, &incx_t, Aiy0_t, &incy_t, Aia0_t, &lda_t);
  }
}

static real_T eml_xnrm2(int32_T n, const emxArray_real_T *x, int32_T ix0)
{
  real_T y;
  ptrdiff_t n_t;
  ptrdiff_t incx_t;
  double * xix0_t;
  if (n < 1) {
    y = 0.0;
  } else {
    n_t = (ptrdiff_t)(n);
    incx_t = (ptrdiff_t)(1);
    xix0_t = (double *)(&x->data[ix0 - 1]);
    y = dnrm2(&n_t, xix0_t, &incx_t);
  }

  return y;
}

static void eml_xscal(int32_T n, real_T a, emxArray_real_T *x, int32_T ix0)
{
  ptrdiff_t n_t;
  ptrdiff_t incx_t;
  double * xix0_t;
  double * a_t;
  if (n < 1) {
  } else {
    n_t = (ptrdiff_t)(n);
    incx_t = (ptrdiff_t)(1);
    xix0_t = (double *)(&x->data[ix0 - 1]);
    a_t = (double *)(&a);
    dscal(&n_t, a_t, xix0_t, &incx_t);
  }
}

static void emlrt_marshallIn(const emlrtStack *sp, const mxArray *d_sprintf,
  const char_T *identifier, char_T y[14])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  b_emlrt_marshallIn(sp, emlrtAlias(d_sprintf), &thisId, y);
  emlrtDestroyArray(&d_sprintf);
}

static void i_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, char_T ret[14])
{
  int32_T iv30[2];
  int32_T i7;
  for (i7 = 0; i7 < 2; i7++) {
    iv30[i7] = 1 + 13 * i7;
  }

  emlrtCheckBuiltInR2012b(sp, msgId, src, "char", FALSE, 2U, iv30);
  emlrtImportCharArray(src, ret, 14);
  emlrtDestroyArray(&src);
}

static void warn_singular(const emlrtStack *sp)
{
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &kh_emlrtRSI;
  eml_warning(&st);
}

void mldivide(const emlrtStack *sp, const emxArray_real_T *A, const
              emxArray_real_T *B, emxArray_real_T *Y)
{
  const mxArray *y;
  static const int32_T iv17[2] = { 1, 21 };

  const mxArray *m6;
  char_T cv22[21];
  int32_T i;
  static const char_T cv23[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'd', 'i', 'm', 'a', 'g', 'r', 'e', 'e' };

  emxArray_real_T *b_B;
  emxArray_real_T *r4;
  uint32_T unnamed_idx_0;
  int32_T loop_ub;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  if (B->size[0] == A->size[0]) {
  } else {
    y = NULL;
    m6 = mxCreateCharArray(2, iv17);
    for (i = 0; i < 21; i++) {
      cv22[i] = cv23[i];
    }

    emlrtInitCharArrayR2013a(sp, 21, m6, cv22);
    emlrtAssign(&y, m6);
    st.site = &ge_emlrtRSI;
    b_error(&st, b_message(&st, y, &o_emlrtMCI), &o_emlrtMCI);
  }

  c_emxInit_real_T(sp, &b_B, 1, &q_emlrtRTEI, TRUE);
  c_emxInit_real_T(sp, &r4, 1, &q_emlrtRTEI, TRUE);
  if ((A->size[0] == 0) || (B->size[0] == 0)) {
    unnamed_idx_0 = (uint32_T)A->size[1];
    i = Y->size[0];
    Y->size[0] = (int32_T)unnamed_idx_0;
    emxEnsureCapacity(sp, (emxArray__common *)Y, i, (int32_T)sizeof(real_T),
                      &q_emlrtRTEI);
    loop_ub = (int32_T)unnamed_idx_0;
    for (i = 0; i < loop_ub; i++) {
      Y->data[i] = 0.0;
    }
  } else if (A->size[0] == A->size[1]) {
    i = Y->size[0];
    Y->size[0] = B->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)Y, i, (int32_T)sizeof(real_T),
                      &q_emlrtRTEI);
    loop_ub = B->size[0];
    for (i = 0; i < loop_ub; i++) {
      Y->data[i] = B->data[i];
    }

    st.site = &ge_emlrtRSI;
    eml_lusolve(&st, A, Y);
  } else {
    i = b_B->size[0];
    b_B->size[0] = B->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)b_B, i, (int32_T)sizeof(real_T),
                      &q_emlrtRTEI);
    loop_ub = B->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_B->data[i] = B->data[i];
    }

    st.site = &ge_emlrtRSI;
    eml_qrsolve(&st, A, b_B, r4);
    i = Y->size[0];
    Y->size[0] = r4->size[0];
    emxEnsureCapacity(sp, (emxArray__common *)Y, i, (int32_T)sizeof(real_T),
                      &q_emlrtRTEI);
    loop_ub = r4->size[0];
    for (i = 0; i < loop_ub; i++) {
      Y->data[i] = r4->data[i];
    }
  }

  emxFree_real_T(&r4);
  emxFree_real_T(&b_B);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (mldivide.c) */
