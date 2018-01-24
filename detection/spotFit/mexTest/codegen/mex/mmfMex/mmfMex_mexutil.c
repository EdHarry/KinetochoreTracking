/*
 * mmfMex_mexutil.c
 *
 * Code generation for function 'mmfMex_mexutil'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "mmfMex_mexutil.h"

/* Function Definitions */
void b_error(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location)
{
  const mxArray *pArray;
  pArray = b;
  emlrtCallMATLABR2012b(sp, 0, NULL, 1, &pArray, "error", TRUE, location);
}

const mxArray *b_message(const emlrtStack *sp, const mxArray *b, emlrtMCInfo
  *location)
{
  const mxArray *pArray;
  const mxArray *m18;
  pArray = b;
  return emlrtCallMATLABR2012b(sp, 1, &m18, 1, &pArray, "message", TRUE,
    location);
}

const mxArray *emlrt_marshallOut(real_T u)
{
  const mxArray *y;
  const mxArray *m14;
  y = NULL;
  m14 = mxCreateDoubleScalar(u);
  emlrtAssign(&y, m14);
  return y;
}

const mxArray *message(const emlrtStack *sp, const mxArray *b, const mxArray *c,
  const mxArray *d, emlrtMCInfo *location)
{
  const mxArray *pArrays[3];
  const mxArray *m17;
  pArrays[0] = b;
  pArrays[1] = c;
  pArrays[2] = d;
  return emlrtCallMATLABR2012b(sp, 1, &m17, 3, pArrays, "message", TRUE,
    location);
}

/* End of code generation (mmfMex_mexutil.c) */
