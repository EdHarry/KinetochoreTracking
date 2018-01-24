/*
 * fitNGaussians3D_mexCode_F_terminate.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode_F_terminate'
 *
 * C source code generated on: Mon May 21 17:33:35 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_F.h"
#include "fitNGaussians3D_mexCode_F_terminate.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

void fitNGaussians3D_mexCode_F_atexit(void)
{
    emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void fitNGaussians3D_mexCode_F_terminate(void)
{
    emlrtLeaveRtStack(&emlrtContextGlobal);
}
/* End of code generation (fitNGaussians3D_mexCode_F_terminate.c) */
