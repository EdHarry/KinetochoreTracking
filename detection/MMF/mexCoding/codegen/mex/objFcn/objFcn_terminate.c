/*
 * objFcn_terminate.c
 *
 * Code generation for function 'objFcn_terminate'
 *
 * C source code generated on: Fri May 25 21:48:51 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "objFcn.h"
#include "objFcn_terminate.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

void objFcn_atexit(void)
{
    emlrtExitTimeCleanup(&emlrtContextGlobal);
}

void objFcn_terminate(void)
{
    emlrtLeaveRtStack(&emlrtContextGlobal);
}
/* End of code generation (objFcn_terminate.c) */
