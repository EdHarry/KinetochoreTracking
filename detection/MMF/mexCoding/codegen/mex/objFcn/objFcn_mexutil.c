/*
 * objFcn_mexutil.c
 *
 * Code generation for function 'objFcn_mexutil'
 *
 * C source code generated on: Fri May 25 21:48:50 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "objFcn.h"
#include "objFcn_mexutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

int32_T _s32_d_(real_T b)
{
    int32_T a;
    a = (int32_T)b;
    if ((real_T)a != (b < 0.0 ? muDoubleScalarCeil((real_T)b) : muDoubleScalarFloor((real_T)b))) {
        emlrtIntegerOverflowErrorR2008a(0);
    }
    return a;
}
/* End of code generation (objFcn_mexutil.c) */
