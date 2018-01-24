/*
 * power.c
 *
 * Code generation for function 'power'
 *
 * C source code generated on: Mon May 28 14:49:15 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_J.h"
#include "power.h"
#include "GaussListND_mexCode.h"
#include "fitNGaussians3D_mexCode_J_emxutil.h"
#include "fitNGaussians3D_mexCode_J_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 *
 */
void power(const emxArray_real_T *a, emxArray_real_T *y)
{
  uint32_T unnamed_idx_0;
  int32_T i2;
  int32_T k;
  unnamed_idx_0 = (uint32_T)a->size[0];
  i2 = y->size[0];
  y->size[0] = (int32_T)unnamed_idx_0;
  emxEnsureCapacity((emxArray__common *)y, i2, (int32_T)sizeof(real_T));
  i2 = y->size[0];
  for (k = 0; k <= i2 - 1; k++) {
    y->data[k] = rt_powd_snf(a->data[k], 2.0);
  }
}

/* End of code generation (power.c) */
