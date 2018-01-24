/*
 * mmfMex.c
 *
 * Code generation for function 'mmfMex'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "mpower.h"
#include "cat.h"
#include "mmfMex_emxutil.h"
#include "prod.h"
#include "diff.h"
#include "erfc.h"
#include "rdivide.h"
#include "repmat.h"
#include "bsxfun.h"
#include "eml_int_forloop_overflow_check.h"
#include "eml_error.h"
#include "repeatEntries_mex.h"
#include "mldivide.h"
#include "mmfMex_mexutil.h"
#include "mmfMex_data.h"

/* Variable Definitions */
static emlrtRSInfo emlrtRSI = { 22, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo b_emlrtRSI = { 24, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo c_emlrtRSI = { 26, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo d_emlrtRSI = { 35, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo e_emlrtRSI = { 59, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo f_emlrtRSI = { 66, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo g_emlrtRSI = { 70, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo h_emlrtRSI = { 71, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo i_emlrtRSI = { 74, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo j_emlrtRSI = { 90, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo k_emlrtRSI = { 94, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo l_emlrtRSI = { 99, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo m_emlrtRSI = { 107, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo n_emlrtRSI = { 109, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo o_emlrtRSI = { 113, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo p_emlrtRSI = { 117, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo q_emlrtRSI = { 120, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo r_emlrtRSI = { 122, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRSInfo kb_emlrtRSI = { 42, "power",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/power.m" };

static emlrtRSInfo lb_emlrtRSI = { 58, "power",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/power.m" };

static emlrtRSInfo dc_emlrtRSI = { 1, "mrdivide",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mrdivide.p" };

static emlrtRSInfo wm_emlrtRSI = { 64, "mtimes",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo xm_emlrtRSI = { 21, "mtimes",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo un_emlrtRSI = { 51, "power",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/power.m" };

static emlrtMCInfo v_emlrtMCI = { 94, 13, "mtimes",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo w_emlrtMCI = { 93, 23, "mtimes",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo x_emlrtMCI = { 99, 13, "mtimes",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtMCInfo y_emlrtMCI = { 98, 23, "mtimes",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRTEInfo emlrtRTEI = { 1, 49, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRTEInfo d_emlrtRTEI = { 20, 1, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRTEInfo e_emlrtRTEI = { 22, 1, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRTEInfo f_emlrtRTEI = { 24, 1, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRTEInfo g_emlrtRTEI = { 35, 5, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRTEInfo h_emlrtRTEI = { 70, 5, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRTEInfo i_emlrtRTEI = { 90, 1, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtRTEInfo j_emlrtRTEI = { 59, 5, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtBCInfo emlrtBCI = { -1, -1, 26, 39, "sigma2", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtBCInfo b_emlrtBCI = { -1, -1, 35, 24, "X", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtECInfo emlrtECI = { 2, 59, 19, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtECInfo b_emlrtECI = { 2, 66, 24, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtECInfo c_emlrtECI = { 2, 66, 40, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtBCInfo c_emlrtBCI = { -1, -1, 74, 23, "gaussListOutput", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtECInfo d_emlrtECI = { -1, 74, 5, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtECInfo e_emlrtECI = { -1, 94, 8, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtBCInfo d_emlrtBCI = { -1, -1, 96, 7, "newAmpBg", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtECInfo f_emlrtECI = { -1, 96, 7, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtBCInfo e_emlrtBCI = { -1, -1, 107, 48, "gaussListOutput", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtBCInfo f_emlrtBCI = { -1, -1, 117, 47, "jacobian", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtBCInfo g_emlrtBCI = { -1, -1, 117, 64, "amp", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtBCInfo h_emlrtBCI = { -1, -1, 117, 20, "jacobian", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtECInfo g_emlrtECI = { -1, 117, 9, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtBCInfo i_emlrtBCI = { -1, -1, 120, 47, "jacobian", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtBCInfo j_emlrtBCI = { -1, -1, 120, 63, "amp", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtBCInfo k_emlrtBCI = { -1, -1, 120, 20, "jacobian", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtECInfo h_emlrtECI = { -1, 120, 9, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtBCInfo l_emlrtBCI = { -1, -1, 122, 44, "jacobian", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtBCInfo m_emlrtBCI = { -1, -1, 122, 103, "X", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtECInfo i_emlrtECI = { 2, 122, 32, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtBCInfo n_emlrtBCI = { -1, -1, 122, 16, "jacobian", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtECInfo j_emlrtECI = { -1, 122, 5, "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m" };

static emlrtBCInfo o_emlrtBCI = { -1, -1, 95, 6, "newAmpBg", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtBCInfo p_emlrtBCI = { -1, -1, 115, 8, "amp", "mmfMex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/mmfMex.m", 0 };

static emlrtRSInfo bq_emlrtRSI = { 93, "mtimes",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo cq_emlrtRSI = { 98, "mtimes",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo nq_emlrtRSI = { 94, "mtimes",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mtimes.m" };

static emlrtRSInfo oq_emlrtRSI = { 99, "mtimes",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/mtimes.m" };

/* Function Declarations */
static void eml_xgemm(int32_T m, int32_T k, const emxArray_real_T *A, int32_T
                      lda, const emxArray_real_T *B, int32_T ldb,
                      emxArray_real_T *C, int32_T ldc);

/* Function Definitions */
static void eml_xgemm(int32_T m, int32_T k, const emxArray_real_T *A, int32_T
                      lda, const emxArray_real_T *B, int32_T ldb,
                      emxArray_real_T *C, int32_T ldc)
{
  real_T alpha1;
  real_T beta1;
  char_T TRANSB;
  char_T TRANSA;
  ptrdiff_t m_t;
  ptrdiff_t n_t;
  ptrdiff_t k_t;
  ptrdiff_t lda_t;
  ptrdiff_t ldb_t;
  ptrdiff_t ldc_t;
  double * alpha1_t;
  double * Aia0_t;
  double * Bib0_t;
  double * beta1_t;
  double * Cic0_t;
  if (m < 1) {
  } else {
    alpha1 = 1.0;
    beta1 = 0.0;
    TRANSB = 'N';
    TRANSA = 'N';
    m_t = (ptrdiff_t)(m);
    n_t = (ptrdiff_t)(1);
    k_t = (ptrdiff_t)(k);
    lda_t = (ptrdiff_t)(lda);
    ldb_t = (ptrdiff_t)(ldb);
    ldc_t = (ptrdiff_t)(ldc);
    alpha1_t = (double *)(&alpha1);
    Aia0_t = (double *)(&A->data[0]);
    Bib0_t = (double *)(&B->data[0]);
    beta1_t = (double *)(&beta1);
    Cic0_t = (double *)(&C->data[0]);
    dgemm(&TRANSA, &TRANSB, &m_t, &n_t, &k_t, alpha1_t, Aia0_t, &lda_t, Bib0_t,
          &ldb_t, beta1_t, Cic0_t, &ldc_t);
  }
}

void mmfMex(const emlrtStack *sp, const emxArray_real_T *maskAmp, const
            emxArray_real_T *maskCoord, const emxArray_real_T *X, const real_T
            sigma_data[3], const int32_T sigma_size[2], emxArray_real_T
            *jacobian, emxArray_real_T *resi, emxArray_real_T *amp, real_T *bg,
            real_T *nGauss, real_T *nDim)
{
  emxArray_real_T *gaussListOutput;
  int32_T i;
  int32_T i0;
  int32_T loop_ub;
  emxArray_real_T *sigma2;
  emxArray_real_T *r0;
  real_T b_maskCoord[2];
  emxArray_real_T b_sigma_data;
  emxArray_real_T *tmp;
  real_T normN;
  real_T sigmaSq_data[3];
  int32_T sigma2_size[2];
  real_T upperCol;
  emxArray_real_T *center2;
  emxArray_real_T *gaussList;
  emxArray_real_T *newAmpBg;
  emxArray_real_T *coordList2;
  emxArray_real_T *b_gaussList;
  emxArray_int32_T *r1;
  emxArray_real_T *r2;
  emxArray_real_T *b_coordList2;
  emxArray_real_T *b_center2;
  emxArray_real_T *c_center2;
  emxArray_real_T *c_maskCoord;
  real_T X_data[3];
  int32_T i1;
  int32_T X_size[2];
  real_T d_maskCoord[2];
  emxArray_real_T b_X_data;
  int32_T d_center2[2];
  int32_T iv0[1];
  const mxArray *y;
  static const int32_T iv1[2] = { 1, 45 };

  const mxArray *m0;
  char_T cv0[45];
  static const char_T cv1[45] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'm', 't', 'i', 'm', 'e', 's', '_', 'n', 'o', 'D',
    'y', 'n', 'a', 'm', 'i', 'c', 'S', 'c', 'a', 'l', 'a', 'r', 'E', 'x', 'p',
    'a', 'n', 's', 'i', 'o', 'n' };

  const mxArray *b_y;
  static const int32_T iv2[2] = { 1, 21 };

  char_T cv2[21];
  static const char_T cv3[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'i', 'n', 'n', 'e', 'r', 'd', 'i', 'm' };

  int32_T b_loop_ub;
  uint32_T ysize[2];
  int32_T sigmaSq_size[2];
  emxArray_real_T *b_gaussListOutput;
  emxArray_real_T *varargin_1;
  uint32_T sz1[2];
  const mxArray *c_y;
  static const int32_T iv3[2] = { 1, 39 };

  char_T cv4[39];
  static const char_T cv5[39] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'c', 'a', 't', 'e', 'n', 'a', 't', 'e', '_', 'd', 'i',
    'm', 'e', 'n', 's', 'i', 'o', 'n', 'M', 'i', 's', 'm', 'a', 't', 'c', 'h' };

  const mxArray *d_y;
  static const int32_T iv4[2] = { 1, 39 };

  int32_T iy;
  boolean_T b0;
  boolean_T b1;
  int32_T cJacobian;
  emxArray_int32_T *r3;
  emxArray_real_T *e_center2;
  int32_T i2;
  real_T b_newAmpBg;
  int32_T iv5[2];
  int32_T iv6[2];
  real_T c_X_data[3];
  int32_T b_X_size[2];
  int32_T iv7[2];
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &b_st;
  d_st.tls = b_st.tls;
  e_st.prev = &st;
  e_st.tls = st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  emxInit_real_T(sp, &gaussListOutput, 2, &d_emlrtRTEI, TRUE);

  /* MMFMEX Summary of this function goes here */
  /*    Detailed explanation goes here */
  /* % INITIAL ALLOCATION */
  /* To get bg from fitting, fit [G|1] * (a,bg)' */
  /* Ones -> To obtain the bg column from the start. gaussList is 1-n */
  /* gauss + one col of ones for bg. */
  /* gaussList = ones(size(maskAmp,1),size(X,1)+1); */
  /* Number of gaussian is equal to the number of spots; */
  /* nGauss = nSpot; */
  *nGauss = X->size[0];
  *nDim = X->size[1];

  /* gaussList = GaussListND_mexCode_(maskCoord,sigma,X); */
  /* nGauss = size(centerInput,1); */
  i = maskCoord->size[0];
  i0 = gaussListOutput->size[0] * gaussListOutput->size[1];
  gaussListOutput->size[0] = i;
  emxEnsureCapacity(sp, (emxArray__common *)gaussListOutput, i0, (int32_T)sizeof
                    (real_T), &emlrtRTEI);
  i = X->size[0] + 1;
  i0 = gaussListOutput->size[0] * gaussListOutput->size[1];
  gaussListOutput->size[1] = i;
  emxEnsureCapacity(sp, (emxArray__common *)gaussListOutput, i0, (int32_T)sizeof
                    (real_T), &emlrtRTEI);
  loop_ub = maskCoord->size[0] * (X->size[0] + 1);
  for (i0 = 0; i0 < loop_ub; i0++) {
    gaussListOutput->data[i0] = 1.0;
  }

  emxInit_real_T(sp, &sigma2, 2, &e_emlrtRTEI, TRUE);
  emxInit_real_T(sp, &r0, 2, &emlrtRTEI, TRUE);
  b_maskCoord[0] = maskCoord->size[0];
  b_maskCoord[1] = 1.0;
  b_sigma_data.data = (real_T *)sigma_data;
  b_sigma_data.size = (int32_T *)sigma_size;
  b_sigma_data.allocatedSize = 3;
  b_sigma_data.numDimensions = 2;
  b_sigma_data.canFreeData = FALSE;
  st.site = &emlrtRSI;
  repmat(&st, &b_sigma_data, b_maskCoord, r0);
  i0 = sigma2->size[0] * sigma2->size[1];
  sigma2->size[0] = r0->size[0];
  sigma2->size[1] = r0->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)sigma2, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = r0->size[0] * r0->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    sigma2->data[i0] = r0->data[i0];
  }

  emxInit_real_T(sp, &tmp, 2, &f_emlrtRTEI, TRUE);
  st.site = &b_emlrtRSI;
  rdivide(&st, sigma2, tmp);
  i0 = sigma2->size[0];
  emlrtDynamicBoundsCheckFastR2012b(1, 1, i0, &emlrtBCI, sp);
  st.site = &c_emlrtRSI;
  st.site = &c_emlrtRSI;
  normN = mpower(0.5 * (real_T)X->size[1]);
  loop_ub = sigma2->size[1];
  sigma2_size[0] = 1;
  sigma2_size[1] = loop_ub;
  for (i0 = 0; i0 < loop_ub; i0++) {
    sigmaSq_data[i0] = sigma2->data[sigma2->size[0] * i0];
  }

  st.site = &c_emlrtRSI;
  upperCol = prod(&st, sigmaSq_data, sigma2_size);
  st.site = &c_emlrtRSI;
  normN *= upperCol;
  i = 0;
  emxInit_real_T(sp, &center2, 2, &g_emlrtRTEI, TRUE);
  b_emxInit_real_T(sp, &gaussList, 3, &h_emlrtRTEI, TRUE);
  c_emxInit_real_T(sp, &newAmpBg, 1, &i_emlrtRTEI, TRUE);
  b_emxInit_real_T(sp, &coordList2, 3, &j_emlrtRTEI, TRUE);
  b_emxInit_real_T(sp, &b_gaussList, 3, &h_emlrtRTEI, TRUE);
  emxInit_int32_T(sp, &r1, 1, &emlrtRTEI, TRUE);
  b_emxInit_real_T(sp, &r2, 3, &emlrtRTEI, TRUE);
  b_emxInit_real_T(sp, &b_coordList2, 3, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &b_center2, 2, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &c_center2, 2, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &c_maskCoord, 2, &emlrtRTEI, TRUE);
  while (i <= X->size[0] - 1) {
    loop_ub = X->size[1];
    i0 = X->size[0];
    i1 = 1 + i;
    i0 = emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &b_emlrtBCI, sp);
    X_size[0] = 1;
    X_size[1] = loop_ub;
    for (i1 = 0; i1 < loop_ub; i1++) {
      X_data[X_size[0] * i1] = X->data[(i0 + X->size[0] * i1) - 1];
    }

    d_maskCoord[0] = maskCoord->size[0];
    d_maskCoord[1] = 1.0;
    b_X_data.data = (real_T *)&X_data;
    b_X_data.size = (int32_T *)&X_size;
    b_X_data.allocatedSize = 3;
    b_X_data.numDimensions = 2;
    b_X_data.canFreeData = FALSE;
    st.site = &d_emlrtRSI;
    repmat(&st, &b_X_data, d_maskCoord, r0);
    i0 = center2->size[0] * center2->size[1];
    center2->size[0] = r0->size[0];
    center2->size[1] = r0->size[1];
    emxEnsureCapacity(sp, (emxArray__common *)center2, i0, (int32_T)sizeof
                      (real_T), &emlrtRTEI);
    loop_ub = r0->size[0] * r0->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      center2->data[i0] = r0->data[i0];
    }

    /* ====================== */
    /*  CALC GAUSSLIST */
    /* ====================== */
    /*  instead of calculating Gauss-values for very complicated geometries, we */
    /*  make a coordinate transformation so that we can use sigma=1 in all */
    /*  dimensions */
    /*  0.5*erfc(-(x+0.5)/sqrt(2))-0.5*erfc(-(x-0.5)/sqrt(2)) gives the integral on the */
    /*  pixel at 1 of a Gaussian with mean 0 and sigma 1 */
    /* center3 = center2(:,1); */
    /* clear center2 */
    /*  convert coordList to 0/1 */
    /* coordList2 = (coordList(1:nCoords,1:nDims) - center2(1:nCoords,1:nDims))./sigma2(1:nCoords,1:nDims); */
    for (i0 = 0; i0 < 2; i0++) {
      sigma2_size[i0] = maskCoord->size[i0];
    }

    for (i0 = 0; i0 < 2; i0++) {
      d_center2[i0] = center2->size[i0];
    }

    emlrtSizeEqCheck2DFastR2012b(sigma2_size, d_center2, &emlrtECI, sp);
    i0 = c_maskCoord->size[0] * c_maskCoord->size[1];
    c_maskCoord->size[0] = maskCoord->size[0];
    c_maskCoord->size[1] = maskCoord->size[1];
    emxEnsureCapacity(sp, (emxArray__common *)c_maskCoord, i0, (int32_T)sizeof
                      (real_T), &emlrtRTEI);
    loop_ub = maskCoord->size[0] * maskCoord->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_maskCoord->data[i0] = maskCoord->data[i0] - center2->data[i0];
    }

    st.site = &e_emlrtRSI;
    b_rdivide(&st, c_maskCoord, sigma2, center2);

    /* clear coordList center3 */
    /*  double coordList as preparation for erfc */
    /* fixed bug: must divide the 0.5 by sigma - KJ */
    /* coordList2 = cat(3,coordList2-0.5./sigma2(1:nCoords,1:nDims), coordList2+0.5./sigma2(1:nCoords,1:nDims)); */
    for (i0 = 0; i0 < 2; i0++) {
      d_center2[i0] = center2->size[i0];
    }

    for (i0 = 0; i0 < 2; i0++) {
      sigma2_size[i0] = tmp->size[i0];
    }

    emlrtSizeEqCheck2DFastR2012b(d_center2, sigma2_size, &b_emlrtECI, sp);
    for (i0 = 0; i0 < 2; i0++) {
      d_center2[i0] = center2->size[i0];
    }

    for (i0 = 0; i0 < 2; i0++) {
      sigma2_size[i0] = tmp->size[i0];
    }

    emlrtSizeEqCheck2DFastR2012b(d_center2, sigma2_size, &c_emlrtECI, sp);
    i0 = b_center2->size[0] * b_center2->size[1];
    b_center2->size[0] = center2->size[0];
    b_center2->size[1] = center2->size[1];
    emxEnsureCapacity(sp, (emxArray__common *)b_center2, i0, (int32_T)sizeof
                      (real_T), &emlrtRTEI);
    loop_ub = center2->size[0] * center2->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_center2->data[i0] = center2->data[i0] - tmp->data[i0];
    }

    i0 = c_center2->size[0] * c_center2->size[1];
    c_center2->size[0] = center2->size[0];
    c_center2->size[1] = center2->size[1];
    emxEnsureCapacity(sp, (emxArray__common *)c_center2, i0, (int32_T)sizeof
                      (real_T), &emlrtRTEI);
    loop_ub = center2->size[0] * center2->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_center2->data[i0] = center2->data[i0] + tmp->data[i0];
    }

    st.site = &f_emlrtRSI;
    cat(&st, b_center2, c_center2, coordList2);

    /*  calculate gaussList */
    /* Jonas was missing the minus sign in erfc. I corrected that - KJ */
    st.site = &g_emlrtRSI;
    b_st.site = &dc_emlrtRSI;
    c_st.site = &ib_emlrtRSI;
    i0 = b_coordList2->size[0] * b_coordList2->size[1] * b_coordList2->size[2];
    b_coordList2->size[0] = coordList2->size[0];
    b_coordList2->size[1] = coordList2->size[1];
    b_coordList2->size[2] = 2;
    emxEnsureCapacity(sp, (emxArray__common *)b_coordList2, i0, (int32_T)sizeof
                      (real_T), &emlrtRTEI);
    loop_ub = coordList2->size[0] * coordList2->size[1] * coordList2->size[2];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_coordList2->data[i0] = -coordList2->data[i0] / 1.4142135623730951;
    }

    st.site = &g_emlrtRSI;
    b_erfc(&st, b_coordList2, coordList2);
    st.site = &g_emlrtRSI;
    i0 = r2->size[0] * r2->size[1] * r2->size[2];
    r2->size[0] = coordList2->size[0];
    r2->size[1] = coordList2->size[1];
    r2->size[2] = 2;
    emxEnsureCapacity(sp, (emxArray__common *)r2, i0, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    loop_ub = coordList2->size[0] * coordList2->size[1] * coordList2->size[2];
    for (i0 = 0; i0 < loop_ub; i0++) {
      r2->data[i0] = 0.5 * coordList2->data[i0];
    }

    st.site = &g_emlrtRSI;
    diff(&st, r2, gaussList);
    st.site = &h_emlrtRSI;
    b_prod(&st, gaussList, b_gaussList);

    /*  norm gaussList */
    loop_ub = gaussListOutput->size[0];
    i0 = r1->size[0];
    r1->size[0] = loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r1, i0, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i0 = 0; i0 < loop_ub; i0++) {
      r1->data[i0] = i0;
    }

    i0 = gaussListOutput->size[1];
    i1 = i + 1;
    emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &c_emlrtBCI, sp);
    st.site = &i_emlrtRSI;
    loop_ub = b_gaussList->size[0];
    i0 = newAmpBg->size[0];
    newAmpBg->size[0] = loop_ub;
    emxEnsureCapacity(&st, (emxArray__common *)newAmpBg, i0, (int32_T)sizeof
                      (real_T), &emlrtRTEI);
    for (i0 = 0; i0 < loop_ub; i0++) {
      newAmpBg->data[i0] = b_gaussList->data[i0] * normN;
    }

    iv0[0] = r1->size[0];
    emlrtSubAssignSizeCheckR2012b(iv0, 1, *(int32_T (*)[1])newAmpBg->size, 1,
      &d_emlrtECI, sp);
    loop_ub = newAmpBg->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      gaussListOutput->data[r1->data[i0] + gaussListOutput->size[0] * i] =
        newAmpBg->data[i0];
    }

    i++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, sp);
  }

  emxFree_real_T(&c_maskCoord);
  emxFree_real_T(&c_center2);
  emxFree_real_T(&b_center2);
  emxFree_real_T(&b_coordList2);
  emxFree_real_T(&r2);
  emxFree_real_T(&r0);
  emxFree_real_T(&b_gaussList);
  emxFree_real_T(&coordList2);
  emxFree_real_T(&gaussList);
  emxFree_real_T(&tmp);
  emxFree_real_T(&sigma2);

  /* % GAUSS LIST CALCULATION */
  /*  for gaussIdx = 1:nGauss */
  /*      tmp = GaussListND_mexCode(maskCoord,sigma,X(gaussIdx,:)); */
  /*      gaussList(:,gaussIdx) = tmp(:,1); */
  /*      %gaussList(:,gaussIdx) = GaussListND_mexCode_mex(maskCoord,sigma,X(gaussIdx,:)); */
  /*  end */
  /* % GAUSS FITTING */
  /* I = (G1,G2,G3,....,1) *(a1;a2;a3;...;bg); */
  /* To fit -> a1;a2..;bg = (G1,G2,G3,...,1) \ I */
  st.site = &j_emlrtRSI;
  mldivide(&st, gaussListOutput, maskAmp, newAmpBg);

  /* % RESIDUAL CALCULATION */
  /* Resi : I - I' -> I' = (G1,G2,G3,G4...1) * (a1;a2;a3;...;bg); */
  st.site = &k_emlrtRSI;
  b_st.site = &xm_emlrtRSI;
  if (!(gaussListOutput->size[1] == newAmpBg->size[0])) {
    if (((gaussListOutput->size[0] == 1) && (gaussListOutput->size[1] == 1)) ||
        (newAmpBg->size[0] == 1)) {
      y = NULL;
      m0 = mxCreateCharArray(2, iv1);
      for (i = 0; i < 45; i++) {
        cv0[i] = cv1[i];
      }

      emlrtInitCharArrayR2013a(&b_st, 45, m0, cv0);
      emlrtAssign(&y, m0);
      c_st.site = &bq_emlrtRSI;
      d_st.site = &nq_emlrtRSI;
      b_error(&c_st, b_message(&d_st, y, &v_emlrtMCI), &w_emlrtMCI);
    } else {
      b_y = NULL;
      m0 = mxCreateCharArray(2, iv2);
      for (i = 0; i < 21; i++) {
        cv2[i] = cv3[i];
      }

      emlrtInitCharArrayR2013a(&b_st, 21, m0, cv2);
      emlrtAssign(&b_y, m0);
      c_st.site = &cq_emlrtRSI;
      d_st.site = &oq_emlrtRSI;
      b_error(&c_st, b_message(&d_st, b_y, &x_emlrtMCI), &y_emlrtMCI);
    }
  }

  if ((gaussListOutput->size[1] == 1) || (newAmpBg->size[0] == 1)) {
    i0 = resi->size[0];
    resi->size[0] = gaussListOutput->size[0];
    emxEnsureCapacity(&st, (emxArray__common *)resi, i0, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    loop_ub = gaussListOutput->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      resi->data[i0] = 0.0;
      b_loop_ub = gaussListOutput->size[1];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        resi->data[i0] += gaussListOutput->data[i0 + gaussListOutput->size[0] *
          i1] * newAmpBg->data[i1];
      }
    }
  } else {
    ysize[0] = (uint32_T)gaussListOutput->size[0];
    i0 = resi->size[0];
    resi->size[0] = (int32_T)ysize[0];
    emxEnsureCapacity(&st, (emxArray__common *)resi, i0, (int32_T)sizeof(real_T),
                      &emlrtRTEI);
    loop_ub = (int32_T)ysize[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      resi->data[i0] = 0.0;
    }

    b_st.site = &wm_emlrtRSI;
    eml_xgemm(gaussListOutput->size[0], gaussListOutput->size[1],
              gaussListOutput, gaussListOutput->size[0], newAmpBg,
              gaussListOutput->size[1], resi, gaussListOutput->size[0]);
  }

  i0 = maskAmp->size[0];
  i1 = resi->size[0];
  emlrtSizeEqCheck1DFastR2012b(i0, i1, &e_emlrtECI, sp);
  i0 = resi->size[0];
  resi->size[0] = maskAmp->size[0];
  emxEnsureCapacity(sp, (emxArray__common *)resi, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  loop_ub = maskAmp->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    resi->data[i0] = maskAmp->data[i0] - resi->data[i0];
  }

  i0 = newAmpBg->size[0];
  i1 = newAmpBg->size[0];
  *bg = newAmpBg->data[emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &o_emlrtBCI,
    sp) - 1];
  if (1 > newAmpBg->size[0] - 1) {
    loop_ub = 0;
  } else {
    i0 = newAmpBg->size[0];
    i1 = newAmpBg->size[0] - 1;
    loop_ub = emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &d_emlrtBCI, sp);
  }

  emlrtVectorVectorIndexCheckR2012b(newAmpBg->size[0], 1, 1, loop_ub,
    &f_emlrtECI, sp);
  i0 = amp->size[0];
  amp->size[0] = loop_ub;
  emxEnsureCapacity(sp, (emxArray__common *)amp, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  for (i0 = 0; i0 < loop_ub; i0++) {
    amp->data[i0] = newAmpBg->data[i0];
  }

  st.site = &l_emlrtRSI;
  b_st.site = &kb_emlrtRSI;
  c_st.site = &un_emlrtRSI;
  for (i0 = 0; i0 < 2; i0++) {
    ysize[i0] = (uint32_T)sigma_size[i0];
  }

  sigmaSq_size[0] = 1;
  sigmaSq_size[1] = (int32_T)ysize[1];
  for (i = 0; i < (int32_T)ysize[1]; i++) {
    c_st.site = &lb_emlrtRSI;
    sigmaSq_data[i] = sigma_data[i] * sigma_data[i];
  }

  /* % ASSEMBLY of JACOBIAN MATRIX */
  /* Compute gradiant of function --> aG(x,y,z)+bg = I; */
  /* Preallocation : jacobian is */
  /* [x1,y1,z1,x2,y2,z2...,xn,yn,zn,a1,a2..,an,bg]; */
  /* jacobian = zeros(size(gaussList,1), ( nGauss*nDim + nGauss + 1) ); */
  /* Increased gaussList */
  if (1 > gaussListOutput->size[1] - 1) {
    b_loop_ub = 0;
  } else {
    i0 = gaussListOutput->size[1];
    i1 = gaussListOutput->size[1] - 1;
    b_loop_ub = emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &e_emlrtBCI, sp);
  }

  emxInit_real_T(sp, &b_gaussListOutput, 2, &emlrtRTEI, TRUE);
  i = gaussListOutput->size[0];
  i0 = b_gaussListOutput->size[0] * b_gaussListOutput->size[1];
  b_gaussListOutput->size[0] = b_loop_ub;
  b_gaussListOutput->size[1] = i;
  emxEnsureCapacity(sp, (emxArray__common *)b_gaussListOutput, i0, (int32_T)
                    sizeof(real_T), &emlrtRTEI);
  for (i0 = 0; i0 < i; i0++) {
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_gaussListOutput->data[i1 + b_gaussListOutput->size[0] * i0] =
        gaussListOutput->data[i0 + gaussListOutput->size[0] * i1];
    }
  }

  emxInit_real_T(sp, &varargin_1, 2, &emlrtRTEI, TRUE);
  st.site = &m_emlrtRSI;
  repeatEntries_mex(&st, b_gaussListOutput, X->size[1], varargin_1);
  i0 = jacobian->size[0] * jacobian->size[1];
  jacobian->size[0] = varargin_1->size[1];
  jacobian->size[1] = varargin_1->size[0];
  emxEnsureCapacity(sp, (emxArray__common *)jacobian, i0, (int32_T)sizeof(real_T),
                    &emlrtRTEI);
  b_loop_ub = varargin_1->size[0];
  emxFree_real_T(&b_gaussListOutput);
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    i = varargin_1->size[1];
    for (i1 = 0; i1 < i; i1++) {
      jacobian->data[i1 + jacobian->size[0] * i0] = varargin_1->data[i0 +
        varargin_1->size[0] * i1];
    }
  }

  /* dI/da = G and dI/dbg = 1 */
  st.site = &n_emlrtRSI;
  i0 = varargin_1->size[0] * varargin_1->size[1];
  varargin_1->size[0] = jacobian->size[0];
  varargin_1->size[1] = jacobian->size[1];
  emxEnsureCapacity(&st, (emxArray__common *)varargin_1, i0, (int32_T)sizeof
                    (real_T), &emlrtRTEI);
  b_loop_ub = jacobian->size[0] * jacobian->size[1];
  for (i0 = 0; i0 < b_loop_ub; i0++) {
    varargin_1->data[i0] = jacobian->data[i0];
  }

  b_st.site = &xb_emlrtRSI;
  for (i0 = 0; i0 < 2; i0++) {
    sz1[i0] = (uint32_T)jacobian->size[i0];
  }

  for (i = 0; i < 2; i++) {
    ysize[i] = sz1[i];
  }

  ysize[1] += gaussListOutput->size[1];
  i0 = jacobian->size[0] * jacobian->size[1];
  jacobian->size[0] = (int32_T)ysize[0];
  jacobian->size[1] = (int32_T)ysize[1];
  emxEnsureCapacity(&st, (emxArray__common *)jacobian, i0, (int32_T)sizeof
                    (real_T), &c_emlrtRTEI);
  if ((varargin_1->size[0] == 0) || (varargin_1->size[1] == 0)) {
    b_st.site = &yb_emlrtRSI;
    eml_error(&b_st);
  }

  b_st.site = &ac_emlrtRSI;
  if (isconsistent(jacobian, varargin_1)) {
  } else {
    c_y = NULL;
    m0 = mxCreateCharArray(2, iv3);
    for (i = 0; i < 39; i++) {
      cv4[i] = cv5[i];
    }

    emlrtInitCharArrayR2013a(&st, 39, m0, cv4);
    emlrtAssign(&c_y, m0);
    b_st.site = &ac_emlrtRSI;
    e_st.site = &iq_emlrtRSI;
    b_error(&b_st, b_message(&e_st, c_y, &m_emlrtMCI), &n_emlrtMCI);
  }

  if (gaussListOutput->size[0] == 0) {
    b_st.site = &yb_emlrtRSI;
    eml_error(&b_st);
  }

  b_st.site = &ac_emlrtRSI;
  if (isconsistent(jacobian, gaussListOutput)) {
  } else {
    d_y = NULL;
    m0 = mxCreateCharArray(2, iv4);
    for (i = 0; i < 39; i++) {
      cv4[i] = cv5[i];
    }

    emlrtInitCharArrayR2013a(&st, 39, m0, cv4);
    emlrtAssign(&d_y, m0);
    b_st.site = &ac_emlrtRSI;
    e_st.site = &iq_emlrtRSI;
    b_error(&b_st, b_message(&e_st, d_y, &m_emlrtMCI), &n_emlrtMCI);
  }

  iy = -1;
  b_loop_ub = varargin_1->size[0] * varargin_1->size[1];
  b_st.site = &bc_emlrtRSI;
  c_st.site = &gb_emlrtRSI;
  if (1 > b_loop_ub) {
    b0 = FALSE;
  } else {
    b0 = (b_loop_ub > 2147483646);
  }

  if (b0) {
    c_st.site = &hb_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  for (i = 1; i <= b_loop_ub; i++) {
    b_st.site = &cc_emlrtRSI;
    iy++;
    jacobian->data[iy] = varargin_1->data[i - 1];
  }

  b_loop_ub = gaussListOutput->size[0] * gaussListOutput->size[1];
  b_st.site = &bc_emlrtRSI;
  c_st.site = &gb_emlrtRSI;
  if (1 > b_loop_ub) {
    b1 = FALSE;
  } else {
    b1 = (b_loop_ub > 2147483646);
  }

  if (b1) {
    c_st.site = &hb_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  for (i = 1; i <= b_loop_ub; i++) {
    b_st.site = &cc_emlrtRSI;
    iy++;
    jacobian->data[iy] = gaussListOutput->data[i - 1];
  }

  emxFree_real_T(&gaussListOutput);

  /* For each coordinate, calculate the first partial derivative--> for */
  /* x,y,z is a*g*(x-c)^2/2*sigma^2 * a*-(x-c)/sigma^2 */
  cJacobian = 0;
  emxInit_int32_T(sp, &r3, 1, &emlrtRTEI, TRUE);
  emxInit_real_T(sp, &e_center2, 2, &emlrtRTEI, TRUE);
  while (cJacobian <= X->size[0] - 1) {
    st.site = &o_emlrtRSI;
    normN = ((1.0 + (real_T)cJacobian) - 1.0) * (real_T)X->size[1];

    /* x..y..z */
    upperCol = ((normN + 1.0) + (real_T)X->size[1]) - 1.0;
    i0 = 1 + cJacobian;
    if (newAmpBg->data[emlrtDynamicBoundsCheckFastR2012b(i0, 1, loop_ub,
         &p_emlrtBCI, sp) - 1] < 0.0) {
      /* jacobian(:,col:upperCol) = bsxfun(@times,jacobian(:,col:upperCol),-amp(cJacobian)); */
      if (normN + 1.0 > upperCol) {
        i0 = 0;
        i1 = 0;
      } else {
        i0 = jacobian->size[1];
        i1 = (int32_T)(normN + 1.0);
        i0 = emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &f_emlrtBCI, sp) - 1;
        i1 = jacobian->size[1];
        i2 = (int32_T)upperCol;
        i1 = emlrtDynamicBoundsCheckFastR2012b(i2, 1, i1, &f_emlrtBCI, sp);
      }

      if (normN + 1.0 > upperCol) {
        i2 = 0;
        i = 0;
      } else {
        i2 = jacobian->size[1];
        i = (int32_T)(normN + 1.0);
        i2 = emlrtDynamicBoundsCheckFastR2012b(i, 1, i2, &h_emlrtBCI, sp) - 1;
        i = jacobian->size[1];
        iy = (int32_T)upperCol;
        i = emlrtDynamicBoundsCheckFastR2012b(iy, 1, i, &h_emlrtBCI, sp);
      }

      b_loop_ub = jacobian->size[0];
      iy = r1->size[0];
      r1->size[0] = b_loop_ub;
      emxEnsureCapacity(sp, (emxArray__common *)r1, iy, (int32_T)sizeof(int32_T),
                        &emlrtRTEI);
      for (iy = 0; iy < b_loop_ub; iy++) {
        r1->data[iy] = iy;
      }

      iy = r3->size[0];
      r3->size[0] = i - i2;
      emxEnsureCapacity(sp, (emxArray__common *)r3, iy, (int32_T)sizeof(int32_T),
                        &emlrtRTEI);
      b_loop_ub = i - i2;
      for (i = 0; i < b_loop_ub; i++) {
        r3->data[i] = i2 + i;
      }

      st.site = &p_emlrtRSI;
      i2 = cJacobian + 1;
      emlrtDynamicBoundsCheckFastR2012b(i2, 1, loop_ub, &g_emlrtBCI, &st);
      b_loop_ub = jacobian->size[0];
      b_newAmpBg = -newAmpBg->data[cJacobian];
      i2 = varargin_1->size[0] * varargin_1->size[1];
      varargin_1->size[0] = b_loop_ub;
      varargin_1->size[1] = i1 - i0;
      emxEnsureCapacity(&st, (emxArray__common *)varargin_1, i2, (int32_T)sizeof
                        (real_T), &emlrtRTEI);
      i = i1 - i0;
      for (i1 = 0; i1 < i; i1++) {
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          varargin_1->data[i2 + varargin_1->size[0] * i1] = jacobian->data[i2 +
            jacobian->size[0] * (i0 + i1)] * b_newAmpBg;
        }
      }

      iv5[0] = r1->size[0];
      iv5[1] = r3->size[0];
      emlrtSubAssignSizeCheckR2012b(iv5, 2, *(int32_T (*)[2])varargin_1->size, 2,
        &g_emlrtECI, sp);
      b_loop_ub = varargin_1->size[1];
      for (i0 = 0; i0 < b_loop_ub; i0++) {
        i = varargin_1->size[0];
        for (i1 = 0; i1 < i; i1++) {
          jacobian->data[r1->data[i1] + jacobian->size[0] * r3->data[i0]] =
            varargin_1->data[i1 + varargin_1->size[0] * i0];
        }
      }
    } else {
      /* jacobian(:,col:upperCol) = bsxfun(@times,jacobian(:,col:upperCol),amp(cJacobian)); */
      if (normN + 1.0 > upperCol) {
        i0 = 0;
        i1 = 0;
      } else {
        i0 = jacobian->size[1];
        i1 = (int32_T)(normN + 1.0);
        i0 = emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &i_emlrtBCI, sp) - 1;
        i1 = jacobian->size[1];
        i2 = (int32_T)upperCol;
        i1 = emlrtDynamicBoundsCheckFastR2012b(i2, 1, i1, &i_emlrtBCI, sp);
      }

      if (normN + 1.0 > upperCol) {
        i2 = 0;
        i = 0;
      } else {
        i2 = jacobian->size[1];
        i = (int32_T)(normN + 1.0);
        i2 = emlrtDynamicBoundsCheckFastR2012b(i, 1, i2, &k_emlrtBCI, sp) - 1;
        i = jacobian->size[1];
        iy = (int32_T)upperCol;
        i = emlrtDynamicBoundsCheckFastR2012b(iy, 1, i, &k_emlrtBCI, sp);
      }

      b_loop_ub = jacobian->size[0];
      iy = r1->size[0];
      r1->size[0] = b_loop_ub;
      emxEnsureCapacity(sp, (emxArray__common *)r1, iy, (int32_T)sizeof(int32_T),
                        &emlrtRTEI);
      for (iy = 0; iy < b_loop_ub; iy++) {
        r1->data[iy] = iy;
      }

      iy = r3->size[0];
      r3->size[0] = i - i2;
      emxEnsureCapacity(sp, (emxArray__common *)r3, iy, (int32_T)sizeof(int32_T),
                        &emlrtRTEI);
      b_loop_ub = i - i2;
      for (i = 0; i < b_loop_ub; i++) {
        r3->data[i] = i2 + i;
      }

      st.site = &q_emlrtRSI;
      i2 = cJacobian + 1;
      emlrtDynamicBoundsCheckFastR2012b(i2, 1, loop_ub, &j_emlrtBCI, &st);
      b_loop_ub = jacobian->size[0];
      b_newAmpBg = newAmpBg->data[cJacobian];
      i2 = varargin_1->size[0] * varargin_1->size[1];
      varargin_1->size[0] = b_loop_ub;
      varargin_1->size[1] = i1 - i0;
      emxEnsureCapacity(&st, (emxArray__common *)varargin_1, i2, (int32_T)sizeof
                        (real_T), &emlrtRTEI);
      i = i1 - i0;
      for (i1 = 0; i1 < i; i1++) {
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          varargin_1->data[i2 + varargin_1->size[0] * i1] = jacobian->data[i2 +
            jacobian->size[0] * (i0 + i1)] * b_newAmpBg;
        }
      }

      iv6[0] = r1->size[0];
      iv6[1] = r3->size[0];
      emlrtSubAssignSizeCheckR2012b(iv6, 2, *(int32_T (*)[2])varargin_1->size, 2,
        &h_emlrtECI, sp);
      b_loop_ub = varargin_1->size[1];
      for (i0 = 0; i0 < b_loop_ub; i0++) {
        i = varargin_1->size[0];
        for (i1 = 0; i1 < i; i1++) {
          jacobian->data[r1->data[i1] + jacobian->size[0] * r3->data[i0]] =
            varargin_1->data[i1 + varargin_1->size[0] * i0];
        }
      }
    }

    if (normN + 1.0 > upperCol) {
      i0 = 0;
      i1 = 0;
    } else {
      i0 = jacobian->size[1];
      i1 = (int32_T)(normN + 1.0);
      i0 = emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &l_emlrtBCI, sp) - 1;
      i1 = jacobian->size[1];
      i2 = (int32_T)upperCol;
      i1 = emlrtDynamicBoundsCheckFastR2012b(i2, 1, i1, &l_emlrtBCI, sp);
    }

    b_loop_ub = jacobian->size[0];
    i2 = varargin_1->size[0] * varargin_1->size[1];
    varargin_1->size[0] = b_loop_ub;
    varargin_1->size[1] = i1 - i0;
    emxEnsureCapacity(sp, (emxArray__common *)varargin_1, i2, (int32_T)sizeof
                      (real_T), &emlrtRTEI);
    i = i1 - i0;
    for (i1 = 0; i1 < i; i1++) {
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        varargin_1->data[i2 + varargin_1->size[0] * i1] = -jacobian->data[i2 +
          jacobian->size[0] * (i0 + i1)];
      }
    }

    b_loop_ub = X->size[1];
    i0 = X->size[0];
    i1 = 1 + cJacobian;
    i0 = emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &m_emlrtBCI, sp);
    b_X_size[0] = 1;
    b_X_size[1] = b_loop_ub;
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      c_X_data[i1] = X->data[(i0 + X->size[0] * i1) - 1];
    }

    st.site = &r_emlrtRSI;
    bsxfun(&st, maskCoord, c_X_data, b_X_size, center2);
    i0 = e_center2->size[0] * e_center2->size[1];
    e_center2->size[0] = center2->size[0];
    e_center2->size[1] = center2->size[1];
    emxEnsureCapacity(sp, (emxArray__common *)e_center2, i0, (int32_T)sizeof
                      (real_T), &emlrtRTEI);
    b_loop_ub = center2->size[0] * center2->size[1];
    for (i0 = 0; i0 < b_loop_ub; i0++) {
      e_center2->data[i0] = center2->data[i0];
    }

    st.site = &r_emlrtRSI;
    b_bsxfun(&st, e_center2, sigmaSq_data, sigmaSq_size, center2);
    for (i0 = 0; i0 < 2; i0++) {
      sigma2_size[i0] = varargin_1->size[i0];
    }

    for (i0 = 0; i0 < 2; i0++) {
      d_center2[i0] = center2->size[i0];
    }

    emlrtSizeEqCheck2DFastR2012b(sigma2_size, d_center2, &i_emlrtECI, sp);
    if (normN + 1.0 > upperCol) {
      i0 = 0;
      i1 = 0;
    } else {
      i0 = jacobian->size[1];
      i1 = (int32_T)(normN + 1.0);
      i0 = emlrtDynamicBoundsCheckFastR2012b(i1, 1, i0, &n_emlrtBCI, sp) - 1;
      i1 = jacobian->size[1];
      i2 = (int32_T)upperCol;
      i1 = emlrtDynamicBoundsCheckFastR2012b(i2, 1, i1, &n_emlrtBCI, sp);
    }

    b_loop_ub = jacobian->size[0];
    i2 = r1->size[0];
    r1->size[0] = b_loop_ub;
    emxEnsureCapacity(sp, (emxArray__common *)r1, i2, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    for (i2 = 0; i2 < b_loop_ub; i2++) {
      r1->data[i2] = i2;
    }

    i2 = r3->size[0];
    r3->size[0] = i1 - i0;
    emxEnsureCapacity(sp, (emxArray__common *)r3, i2, (int32_T)sizeof(int32_T),
                      &emlrtRTEI);
    b_loop_ub = i1 - i0;
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      r3->data[i1] = i0 + i1;
    }

    i0 = center2->size[0] * center2->size[1];
    center2->size[0] = varargin_1->size[0];
    center2->size[1] = varargin_1->size[1];
    emxEnsureCapacity(sp, (emxArray__common *)center2, i0, (int32_T)sizeof
                      (real_T), &emlrtRTEI);
    b_loop_ub = varargin_1->size[0] * varargin_1->size[1];
    for (i0 = 0; i0 < b_loop_ub; i0++) {
      center2->data[i0] *= varargin_1->data[i0];
    }

    iv7[0] = r1->size[0];
    iv7[1] = r3->size[0];
    emlrtSubAssignSizeCheckR2012b(iv7, 2, *(int32_T (*)[2])center2->size, 2,
      &j_emlrtECI, sp);
    b_loop_ub = center2->size[1];
    for (i0 = 0; i0 < b_loop_ub; i0++) {
      i = center2->size[0];
      for (i1 = 0; i1 < i; i1++) {
        jacobian->data[r1->data[i1] + jacobian->size[0] * r3->data[i0]] =
          center2->data[i1 + center2->size[0] * i0];
      }
    }

    /* tmp = -jacobian(:,col:upperCol) .* ((maskCoord - X(cJacobian,:)) ./ sigmaSq); */
    /* jacobian(:,col:upperCol) = tmp(:,1); */
    cJacobian++;
    emlrtBreakCheckFastR2012b(emlrtBreakCheckR2012bFlagVar, sp);
  }

  emxFree_real_T(&e_center2);
  emxFree_real_T(&varargin_1);
  emxFree_int32_T(&r3);
  emxFree_int32_T(&r1);
  emxFree_real_T(&newAmpBg);
  emxFree_real_T(&center2);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (mmfMex.c) */
