/*
 * GaussListND_mexCode_terminate.c
 *
 * Code generation for function 'GaussListND_mexCode_terminate'
 *
 * C source code generated on: Sun Dec  2 22:57:02 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GaussListND_mexCode.h"
#include "GaussListND_mexCode_terminate.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void GaussListND_mexCode_atexit(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void GaussListND_mexCode_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (GaussListND_mexCode_terminate.c) */
