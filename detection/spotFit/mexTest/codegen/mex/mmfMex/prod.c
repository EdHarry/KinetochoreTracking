/*
 * prod.c
 *
 * Code generation for function 'prod'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "prod.h"
#include "mmfMex_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "mmfMex_mexutil.h"
#include "mmfMex_data.h"

/* Variable Definitions */
static emlrtRSInfo nb_emlrtRSI = { 17, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo ob_emlrtRSI = { 59, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo pb_emlrtRSI = { 60, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo ud_emlrtRSI = { 64, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo vd_emlrtRSI = { 68, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo wd_emlrtRSI = { 70, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo xd_emlrtRSI = { 71, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo yd_emlrtRSI = { 74, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo ae_emlrtRSI = { 75, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo be_emlrtRSI = { 76, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo ce_emlrtRSI = { 78, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtMCInfo g_emlrtMCI = { 18, 9, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtMCInfo h_emlrtMCI = { 17, 19, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtMCInfo i_emlrtMCI = { 23, 9, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtMCInfo j_emlrtMCI = { 20, 19, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRTEInfo o_emlrtRTEI = { 53, 1, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRTEInfo p_emlrtRTEI = { 1, 14, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo yp_emlrtRSI = { 20, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo kq_emlrtRSI = { 23, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

static emlrtRSInfo lq_emlrtRSI = { 18, "prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/datafun/prod.m" };

/* Function Definitions */
void b_prod(const emlrtStack *sp, const emxArray_real_T *x, emxArray_real_T *y)
{
  uint32_T sz[3];
  int32_T vstride;
  int32_T iy;
  int32_T ixstart;
  int32_T j;
  int32_T ix;
  real_T p;
  boolean_T overflow;
  int32_T k;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  for (vstride = 0; vstride < 3; vstride++) {
    sz[vstride] = (uint32_T)x->size[vstride];
  }

  vstride = y->size[0] * y->size[1] * y->size[2];
  y->size[0] = (int32_T)sz[0];
  y->size[1] = 1;
  y->size[2] = 1;
  emxEnsureCapacity(sp, (emxArray__common *)y, vstride, (int32_T)sizeof(real_T),
                    &o_emlrtRTEI);
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    vstride = y->size[0] * y->size[1] * y->size[2];
    y->size[0] = (int32_T)sz[0];
    y->size[1] = 1;
    y->size[2] = 1;
    emxEnsureCapacity(sp, (emxArray__common *)y, vstride, (int32_T)sizeof(real_T),
                      &p_emlrtRTEI);
    iy = (int32_T)sz[0];
    for (vstride = 0; vstride < iy; vstride++) {
      y->data[vstride] = 1.0;
    }
  } else {
    st.site = &ud_emlrtRSI;
    vstride = x->size[0];
    iy = -1;
    st.site = &vd_emlrtRSI;
    ixstart = -1;
    st.site = &wd_emlrtRSI;
    if (vstride > 2147483646) {
      b_st.site = &hb_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }

    for (j = 1; j <= vstride; j++) {
      st.site = &xd_emlrtRSI;
      ixstart++;
      ix = ixstart;
      p = x->data[ixstart];
      st.site = &yd_emlrtRSI;
      if (2 > x->size[1]) {
        overflow = FALSE;
      } else {
        overflow = (x->size[1] > 2147483646);
      }

      if (overflow) {
        b_st.site = &hb_emlrtRSI;
        check_forloop_overflow_error(&b_st);
      }

      for (k = 2; k <= x->size[1]; k++) {
        st.site = &ae_emlrtRSI;
        ix += vstride;
        st.site = &be_emlrtRSI;
        p *= x->data[ix];
      }

      st.site = &ce_emlrtRSI;
      iy++;
      y->data[iy] = p;
    }
  }
}

real_T prod(const emlrtStack *sp, const real_T x_data[3], const int32_T x_size[2])
{
  real_T y;
  boolean_T p;
  boolean_T b_p;
  int32_T i;
  int32_T exitg1;
  const mxArray *b_y;
  static const int32_T iv12[2] = { 1, 31 };

  const mxArray *m3;
  char_T cv14[31];
  static const char_T cv15[31] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'p', 'r', 'o', 'd', '_', 's', 'p', 'e', 'c', 'i',
    'a', 'l', 'E', 'm', 'p', 't', 'y' };

  const mxArray *c_y;
  static const int32_T iv13[2] = { 1, 36 };

  char_T cv16[36];
  static const char_T cv17[36] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'a', 'u', 't', 'o', 'D', 'i', 'm', 'I', 'n', 'c',
    'o', 'm', 'p', 'a', 't', 'i', 'b', 'i', 'l', 'i', 't', 'y' };

  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = sp;
  b_st.tls = sp->tls;
  st.site = &nb_emlrtRSI;
  p = FALSE;
  b_p = FALSE;
  i = 0;
  do {
    exitg1 = 0;
    if (i < 2) {
      if (x_size[i] != 0) {
        exitg1 = 1;
      } else {
        i++;
      }
    } else {
      b_p = TRUE;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  if (!b_p) {
  } else {
    p = TRUE;
  }

  if (!p) {
  } else {
    b_y = NULL;
    m3 = mxCreateCharArray(2, iv12);
    for (i = 0; i < 31; i++) {
      cv14[i] = cv15[i];
    }

    emlrtInitCharArrayR2013a(sp, 31, m3, cv14);
    emlrtAssign(&b_y, m3);
    st.site = &nb_emlrtRSI;
    b_st.site = &lq_emlrtRSI;
    b_error(&st, b_message(&b_st, b_y, &g_emlrtMCI), &h_emlrtMCI);
  }

  if ((x_size[1] == 1) || (x_size[1] != 1)) {
    p = TRUE;
  } else {
    p = FALSE;
  }

  if (p) {
  } else {
    c_y = NULL;
    m3 = mxCreateCharArray(2, iv13);
    for (i = 0; i < 36; i++) {
      cv16[i] = cv17[i];
    }

    emlrtInitCharArrayR2013a(sp, 36, m3, cv16);
    emlrtAssign(&c_y, m3);
    st.site = &yp_emlrtRSI;
    b_st.site = &kq_emlrtRSI;
    b_error(&st, b_message(&b_st, c_y, &i_emlrtMCI), &j_emlrtMCI);
  }

  if (x_size[1] == 0) {
    y = 1.0;
  } else {
    y = x_data[0];
    st.site = &ob_emlrtRSI;
    for (i = 2; i <= x_size[1]; i++) {
      st.site = &pb_emlrtRSI;
      y *= x_data[i - 1];
    }
  }

  return y;
}

/* End of code generation (prod.c) */
