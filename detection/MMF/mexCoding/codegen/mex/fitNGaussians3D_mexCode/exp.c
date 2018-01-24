/*
 * exp.c
 *
 * Code generation for function 'exp'
 *
 * C source code generated on: Tue Nov 19 11:16:20 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "exp.h"

/* Function Definitions */
void b_exp(emxArray_real_T *x)
{
  int32_T i6;
  int32_T k;
  i6 = x->size[0];
  for (k = 0; k < i6; k++) {
    x->data[k] = muDoubleScalarExp(x->data[k]);
  }
}

/* End of code generation (exp.c) */
