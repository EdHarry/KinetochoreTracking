/*
 * diff.c
 *
 * Code generation for function 'diff'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "diff.h"
#include "mmfMex_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "mmfMex_data.h"

/* Variable Definitions */
static emlrtRSInfo md_emlrtRSI = { 100, "diff",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/diff.m" };

static emlrtRSInfo nd_emlrtRSI = { 84, "diff",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/diff.m" };

static emlrtRSInfo od_emlrtRSI = { 81, "diff",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/diff.m" };

static emlrtRSInfo pd_emlrtRSI = { 77, "diff",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/diff.m" };

static emlrtRSInfo qd_emlrtRSI = { 76, "diff",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/diff.m" };

static emlrtRSInfo rd_emlrtRSI = { 74, "diff",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/diff.m" };

static emlrtRTEInfo n_emlrtRTEI = { 72, 1, "diff",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/diff.m" };

/* Function Definitions */
void diff(const emlrtStack *sp, const emxArray_real_T *x, emxArray_real_T *y)
{
  uint32_T ySize[3];
  int32_T ix;
  int32_T stride;
  int32_T iy;
  boolean_T b5;
  int32_T s;
  real_T work;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  for (ix = 0; ix < 3; ix++) {
    ySize[ix] = (uint32_T)x->size[ix];
  }

  ix = y->size[0] * y->size[1] * y->size[2];
  y->size[0] = (int32_T)ySize[0];
  y->size[1] = (int32_T)ySize[1];
  y->size[2] = 1;
  emxEnsureCapacity(sp, (emxArray__common *)y, ix, (int32_T)sizeof(real_T),
                    &n_emlrtRTEI);
  if (!(((int32_T)ySize[0] == 0) || ((int32_T)ySize[1] == 0))) {
    st.site = &rd_emlrtRSI;
    stride = 1;
    for (ix = 0; ix < 2; ix++) {
      stride *= x->size[ix];
    }

    st.site = &qd_emlrtRSI;
    st.site = &pd_emlrtRSI;
    st.site = &od_emlrtRSI;
    ix = 0;
    iy = 0;
    st.site = &nd_emlrtRSI;
    if (1 > stride) {
      b5 = FALSE;
    } else {
      b5 = (stride > 2147483646);
    }

    if (b5) {
      b_st.site = &hb_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }

    for (s = 1; s <= stride; s++) {
      work = x->data[ix];
      st.site = &md_emlrtRSI;
      work = x->data[ix + stride] - work;
      y->data[iy] = work;
      ix++;
      iy++;
    }
  }
}

/* End of code generation (diff.c) */
