/*
 * colon.c
 *
 * Code generation for function 'colon'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "colon.h"
#include "mmfMex_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "mmfMex_data.h"

/* Variable Definitions */
static emlrtRSInfo qf_emlrtRSI = { 157, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo rf_emlrtRSI = { 156, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo sf_emlrtRSI = { 151, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRTEInfo s_emlrtRTEI = { 152, 1, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

/* Function Definitions */
void eml_signed_integer_colon(const emlrtStack *sp, int32_T b, emxArray_int32_T *
  y)
{
  int32_T yk;
  boolean_T b6;
  int32_T k;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  st.site = &sf_emlrtRSI;
  yk = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = b;
  emxEnsureCapacity(sp, (emxArray__common *)y, yk, (int32_T)sizeof(int32_T),
                    &s_emlrtRTEI);
  y->data[0] = 1;
  yk = 1;
  st.site = &rf_emlrtRSI;
  if (2 > b) {
    b6 = FALSE;
  } else {
    b6 = (b > 2147483646);
  }

  if (b6) {
    b_st.site = &hb_emlrtRSI;
    check_forloop_overflow_error(&b_st);
  }

  for (k = 2; k <= b; k++) {
    st.site = &qf_emlrtRSI;
    yk++;
    y->data[k - 1] = yk;
  }
}

/* End of code generation (colon.c) */
