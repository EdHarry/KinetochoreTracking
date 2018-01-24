/*
 * exp.c
 *
 * Code generation for function 'exp'
 *
 * C source code generated on: Mon May 21 19:42:21 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_J.h"
#include "exp.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 * 
 */
void b_exp(emxArray_real_T *x)
{
    int32_T loop_ub;
    uint32_T k;
    loop_ub = x->size[0];
    for (k = 1U; k <= (uint32_T)loop_ub; k++) {
        x->data[(int32_T)k - 1] = muDoubleScalarExp(x->data[(int32_T)k - 1]);
    }
}
/* End of code generation (exp.c) */
