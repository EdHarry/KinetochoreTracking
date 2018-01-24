/*
 * rdivide.c
 *
 * Code generation for function 'rdivide'
 *
 * C source code generated on: Tue Nov 19 11:16:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "rdivide.h"
#include "fitNGaussians3D_mexCode_emxutil.h"

/* Function Definitions */
void rdivide(const emxArray_real_T *x, real_T y, emxArray_real_T *z)
{
  int32_T i2;
  int32_T loop_ub;
  i2 = z->size[0];
  z->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)z, i2, (int32_T)sizeof(real_T));
  loop_ub = x->size[0];
  for (i2 = 0; i2 < loop_ub; i2++) {
    z->data[i2] = x->data[i2] / y;
  }
}

/* End of code generation (rdivide.c) */
