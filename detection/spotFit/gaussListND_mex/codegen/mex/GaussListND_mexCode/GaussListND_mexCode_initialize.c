/*
 * GaussListND_mexCode_initialize.c
 *
 * Code generation for function 'GaussListND_mexCode_initialize'
 *
 * C source code generated on: Sun Dec  2 22:57:02 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GaussListND_mexCode.h"
#include "GaussListND_mexCode_initialize.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
static const volatile char_T *emlrtBreakCheckR2012bFlagVar;

/* Function Declarations */

/* Function Definitions */
void GaussListND_mexCode_initialize(emlrtContext *aContext)
{
  emlrtBreakCheckR2012bFlagVar = emlrtGetBreakCheckFlagAddressR2012b();
  emlrtCreateRootTLS(&emlrtRootTLSGlobal, aContext, NULL, 1);
  emlrtClearAllocCountR2012b(emlrtRootTLSGlobal, FALSE, 0U, 0);
  emlrtEnterRtStackR2012b(emlrtRootTLSGlobal);
  emlrtFirstTimeR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (GaussListND_mexCode_initialize.c) */
