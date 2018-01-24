/*
 * rdivide.c
 *
 * Code generation for function 'rdivide'
 *
 * C source code generated on: Mon May 28 14:49:14 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_J.h"
#include "rdivide.h"
#include "fitNGaussians3D_mexCode_J_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 *
 */
void rdivide(const emxArray_real_T *x, real_T y, emxArray_real_T *z)
{
  int32_T i1;
  int32_T loop_ub;
  i1 = z->size[0];
  z->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)z, i1, (int32_T)sizeof(real_T));
  loop_ub = x->size[0] - 1;
  for (i1 = 0; i1 <= loop_ub; i1++) {
    z->data[i1] = x->data[i1] / y;
  }
}

/* End of code generation (rdivide.c) */
