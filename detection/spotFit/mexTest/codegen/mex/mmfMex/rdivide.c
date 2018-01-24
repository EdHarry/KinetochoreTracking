/*
 * rdivide.c
 *
 * Code generation for function 'rdivide'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "rdivide.h"
#include "mmfMex_emxutil.h"
#include "mmfMex_mexutil.h"
#include "mmfMex_data.h"

/* Variable Definitions */
static emlrtRSInfo ub_emlrtRSI = { 13, "rdivide",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/rdivide.m" };

static emlrtMCInfo k_emlrtMCI = { 14, 5, "rdivide",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/rdivide.m" };

static emlrtMCInfo l_emlrtMCI = { 13, 15, "rdivide",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/rdivide.m" };

static emlrtRTEInfo l_emlrtRTEI = { 1, 14, "rdivide",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/rdivide.m" };

static emlrtRSInfo jq_emlrtRSI = { 14, "rdivide",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/rdivide.m" };

/* Function Definitions */
void b_rdivide(const emlrtStack *sp, const emxArray_real_T *x, const
               emxArray_real_T *y, emxArray_real_T *z)
{
  uint32_T varargin_1[2];
  int32_T i;
  uint32_T varargin_2[2];
  boolean_T p;
  boolean_T b_p;
  boolean_T exitg1;
  const mxArray *b_y;
  static const int32_T iv14[2] = { 1, 15 };

  const mxArray *m4;
  char_T cv18[15];
  static const char_T cv19[15] = { 'M', 'A', 'T', 'L', 'A', 'B', ':', 'd', 'i',
    'm', 'a', 'g', 'r', 'e', 'e' };

  int32_T loop_ub;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = sp;
  b_st.tls = sp->tls;
  st.site = &ub_emlrtRSI;
  for (i = 0; i < 2; i++) {
    varargin_1[i] = (uint32_T)x->size[i];
  }

  for (i = 0; i < 2; i++) {
    varargin_2[i] = (uint32_T)y->size[i];
  }

  p = FALSE;
  b_p = TRUE;
  i = 0;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (i < 2)) {
    if (!(varargin_1[i] == varargin_2[i])) {
      b_p = FALSE;
      exitg1 = TRUE;
    } else {
      i++;
    }
  }

  if (!b_p) {
  } else {
    p = TRUE;
  }

  if (p) {
  } else {
    b_y = NULL;
    m4 = mxCreateCharArray(2, iv14);
    for (i = 0; i < 15; i++) {
      cv18[i] = cv19[i];
    }

    emlrtInitCharArrayR2013a(sp, 15, m4, cv18);
    emlrtAssign(&b_y, m4);
    st.site = &ub_emlrtRSI;
    b_st.site = &jq_emlrtRSI;
    b_error(&st, b_message(&b_st, b_y, &k_emlrtMCI), &l_emlrtMCI);
  }

  st.site = &ib_emlrtRSI;
  i = z->size[0] * z->size[1];
  z->size[0] = x->size[0];
  z->size[1] = x->size[1];
  emxEnsureCapacity(&st, (emxArray__common *)z, i, (int32_T)sizeof(real_T),
                    &l_emlrtRTEI);
  loop_ub = x->size[0] * x->size[1];
  for (i = 0; i < loop_ub; i++) {
    z->data[i] = x->data[i] / y->data[i];
  }
}

void rdivide(const emlrtStack *sp, const emxArray_real_T *y, emxArray_real_T *z)
{
  int32_T i3;
  int32_T loop_ub;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &ib_emlrtRSI;
  i3 = z->size[0] * z->size[1];
  z->size[0] = y->size[0];
  z->size[1] = y->size[1];
  emxEnsureCapacity(&st, (emxArray__common *)z, i3, (int32_T)sizeof(real_T),
                    &l_emlrtRTEI);
  loop_ub = y->size[0] * y->size[1];
  for (i3 = 0; i3 < loop_ub; i3++) {
    z->data[i3] = 0.5 / y->data[i3];
  }
}

/* End of code generation (rdivide.c) */
