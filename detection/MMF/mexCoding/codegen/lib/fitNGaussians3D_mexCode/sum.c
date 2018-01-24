/*
 * sum.c
 *
 * Code generation for function 'sum'
 *
 * C source code generated on: Thu May  3 12:56:20 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "sum.h"
#include "fitNGaussians3D_mexCode_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void sum(const emxArray_real_T *x, emxArray_real_T *y)
{
  uint32_T sz[2];
  int32_T ixstart;
  int32_T k;
  int32_T vlen;
  int32_T npages;
  int32_T ix;
  int32_T iy;
  int32_T i;
  real_T s;
  for (ixstart = 0; ixstart < 2; ixstart++) {
    sz[ixstart] = (uint32_T)x->size[ixstart];
  }

  sz[0] = 1U;
  ixstart = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int32_T)sz[1];
  emxEnsureCapacity((emxArray__common *)y, ixstart, (int32_T)sizeof(real_T));
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    ixstart = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)y, ixstart, (int32_T)sizeof(real_T));
    k = y->size[1] - 1;
    for (ixstart = 0; ixstart <= k; ixstart++) {
      y->data[y->size[0] * ixstart] = 0.0;
    }
  } else {
    vlen = x->size[0];
    npages = x->size[1];
    ix = -1;
    iy = -1;
    for (i = 1; i <= npages; i++) {
      ixstart = ix + 1;
      ix++;
      s = x->data[ixstart];
      for (k = 2; k <= vlen; k++) {
        ix++;
        s += x->data[ix];
      }

      iy++;
      y->data[iy] = s;
    }
  }
}

/* End of code generation (sum.c) */
