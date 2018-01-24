/*
 * fitNGaussians3D_mexCode_J_terminate.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode_J_terminate'
 *
 * C source code generated on: Mon May 21 19:42:21 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_J.h"
#include "fitNGaussians3D_mexCode_J_terminate.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

void fitNGaussians3D_mexCode_J_atexit(void)
{
    emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void fitNGaussians3D_mexCode_J_terminate(void)
{
    emlrtLeaveRtStack(&emlrtContextGlobal);
}
/* End of code generation (fitNGaussians3D_mexCode_J_terminate.c) */
