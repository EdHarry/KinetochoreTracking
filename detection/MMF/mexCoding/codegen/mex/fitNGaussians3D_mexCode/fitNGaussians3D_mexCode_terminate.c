/*
 * fitNGaussians3D_mexCode_terminate.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode_terminate'
 *
 * C source code generated on: Tue Nov 19 11:16:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "fitNGaussians3D_mexCode_terminate.h"

/* Function Definitions */
void fitNGaussians3D_mexCode_atexit(void)
{
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, &emlrtContextGlobal, NULL, 1);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

void fitNGaussians3D_mexCode_terminate(void)
{
  emlrtLeaveRtStackR2012b(emlrtRootTLSGlobal);
  emlrtDestroyRootTLS(&emlrtRootTLSGlobal);
}

/* End of code generation (fitNGaussians3D_mexCode_terminate.c) */
