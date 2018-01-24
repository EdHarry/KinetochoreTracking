/*
 * exp.c
 *
 * Code generation for function 'exp'
 *
 * C source code generated on: Thu May  3 12:56:20 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "exp.h"
#include "fitNGaussians3D_mexCode_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_exp(emxArray_real_T *x)
{
  emxArray_real_T *b_x;
  int32_T i5;
  int32_T k;
  c_emxInit_real_T(&b_x, 1);
  i5 = x->size[0];
  for (k = 0; k <= i5 - 1; k++) {
    x->data[(int32_T)(1.0 + (real_T)k) - 1] = exp(x->data[(int32_T)(1.0 +
      (real_T)k) - 1]);
  }

  emxFree_real_T(&b_x);
}

/* End of code generation (exp.c) */
