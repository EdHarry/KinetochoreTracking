/*
 * eml_warning.c
 *
 * Code generation for function 'eml_warning'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "eml_warning.h"
#include "repmat.h"
#include "mmfMex_mexutil.h"

/* Variable Definitions */
static emlrtMCInfo r_emlrtMCI = { 16, 13, "eml_warning",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_warning.m" };

static emlrtMCInfo s_emlrtMCI = { 16, 5, "eml_warning",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_warning.m" };

static emlrtRSInfo hq_emlrtRSI = { 16, "eml_warning",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_warning.m" };

/* Function Declarations */
static void warning(const emlrtStack *sp, const mxArray *b, emlrtMCInfo
                    *location);

/* Function Definitions */
static void warning(const emlrtStack *sp, const mxArray *b, emlrtMCInfo
                    *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(sp, 0, NULL, 1, &pArray, "warning", TRUE, location);
}

void b_eml_warning(const emlrtStack *sp, real_T varargin_2, const char_T
                   varargin_3[14])
{
  const mxArray *y;
  static const int32_T iv20[2] = { 1, 32 };

  const mxArray *m9;
  char_T cv29[32];
  int32_T i;
  static const char_T cv30[32] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'r', 'a', 'n', 'k', 'D', 'e', 'f', 'i', 'c', 'i', 'e',
    'n', 't', 'M', 'a', 't', 'r', 'i', 'x' };

  const mxArray *b_y;
  static const int32_T iv21[2] = { 1, 14 };

  char_T b_varargin_3[14];
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  y = NULL;
  m9 = mxCreateCharArray(2, iv20);
  for (i = 0; i < 32; i++) {
    cv29[i] = cv30[i];
  }

  emlrtInitCharArrayR2013a(sp, 32, m9, cv29);
  emlrtAssign(&y, m9);
  b_y = NULL;
  m9 = mxCreateCharArray(2, iv21);
  for (i = 0; i < 14; i++) {
    b_varargin_3[i] = varargin_3[i];
  }

  emlrtInitCharArrayR2013a(sp, 14, m9, b_varargin_3);
  emlrtAssign(&b_y, m9);
  st.site = &hq_emlrtRSI;
  warning(&st, message(&st, y, emlrt_marshallOut(varargin_2), b_y, &r_emlrtMCI),
          &s_emlrtMCI);
}

void eml_warning(const emlrtStack *sp)
{
  const mxArray *y;
  static const int32_T iv18[2] = { 1, 27 };

  const mxArray *m7;
  char_T cv24[27];
  int32_T i;
  static const char_T cv25[27] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 's', 'i', 'n', 'g', 'u', 'l', 'a', 'r', 'M', 'a', 't',
    'r', 'i', 'x' };

  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  y = NULL;
  m7 = mxCreateCharArray(2, iv18);
  for (i = 0; i < 27; i++) {
    cv24[i] = cv25[i];
  }

  emlrtInitCharArrayR2013a(sp, 27, m7, cv24);
  emlrtAssign(&y, m7);
  st.site = &hq_emlrtRSI;
  warning(&st, b_message(&st, y, &r_emlrtMCI), &s_emlrtMCI);
}

/* End of code generation (eml_warning.c) */
