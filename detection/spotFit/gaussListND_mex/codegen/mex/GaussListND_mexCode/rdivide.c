/*
 * rdivide.c
 *
 * Code generation for function 'rdivide'
 *
 * C source code generated on: Sun Dec  2 22:57:02 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GaussListND_mexCode.h"
#include "rdivide.h"
#include "GaussListND_mexCode_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void rdivide(real_T x, const emxArray_real_T *y, emxArray_real_T *z)
{
  int32_T i1;
  int32_T loop_ub;
  i1 = z->size[0] * z->size[1];
  z->size[0] = y->size[0];
  z->size[1] = y->size[1];
  emxEnsureCapacity((emxArray__common *)z, i1, (int32_T)sizeof(real_T));
  loop_ub = y->size[0] * y->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    z->data[i1] = 0.5 / y->data[i1];
  }
}

/* End of code generation (rdivide.c) */
