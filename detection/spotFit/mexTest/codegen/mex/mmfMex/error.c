/*
 * error.c
 *
 * Code generation for function 'error'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "error.h"
#include "mmfMex_mexutil.h"

/* Variable Definitions */
static emlrtMCInfo ab_emlrtMCI = { 16, 1, "error",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/lang/error.m" };

static emlrtRSInfo wp_emlrtRSI = { 16, "error",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/lang/error.m" };

/* Function Definitions */
void error(const emlrtStack *sp)
{
  const mxArray *y;
  static const int32_T iv23[2] = { 1, 41 };

  const mxArray *m11;
  char_T cv33[41];
  int32_T i;
  static const char_T cv34[41] = { 't', 'w', 'o', ' ', 'n', 'o', 'n', '-', 'e',
    'm', 'p', 't', 'y', ' ', 'i', 'n', 'p', 'u', 't', ' ', 'a', 'r', 'g', 'u',
    'm', 'e', 'n', 't', 's', ' ', 'a', 'r', 'e', ' ', 'n', 'e', 'e', 'd', 'e',
    'd', '!' };

  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  y = NULL;
  m11 = mxCreateCharArray(2, iv23);
  for (i = 0; i < 41; i++) {
    cv33[i] = cv34[i];
  }

  emlrtInitCharArrayR2013a(sp, 41, m11, cv33);
  emlrtAssign(&y, m11);
  st.site = &wp_emlrtRSI;
  b_error(&st, y, &ab_emlrtMCI);
}

/* End of code generation (error.c) */
