/*
 * squeeze.c
 *
 * Code generation for function 'squeeze'
 *
 * C source code generated on: Thu May  3 12:56:20 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "squeeze.h"
#include "fitNGaussians3D_mexCode_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void squeeze(const emxArray_real_T *a, emxArray_real_T *b)
{
  int32_T sqsz[3];
  int32_T k;
  sqsz[0] = a->size[0];
  sqsz[1] = 1;
  sqsz[2] = 1;
  k = b->size[0] * b->size[1];
  b->size[0] = sqsz[0];
  b->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)b, k, (int32_T)sizeof(real_T));
  for (k = 0; k + 1 <= a->size[0]; k++) {
    b->data[k] = a->data[k];
  }
}

/* End of code generation (squeeze.c) */
