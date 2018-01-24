/*
 * cat.c
 *
 * Code generation for function 'cat'
 *
 * C source code generated on: Tue Nov 19 11:16:20 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "cat.h"
#include "fitNGaussians3D_mexCode_emxutil.h"

/* Function Definitions */
void cat(const emxArray_real_T *varargin_1, const emxArray_real_T *varargin_2,
         emxArray_real_T *y)
{
  uint32_T ysize[3];
  int32_T j;
  uint32_T sz1[2];
  int32_T iy;
  for (j = 0; j < 3; j++) {
    ysize[j] = 1U;
  }

  sz1[0] = (uint32_T)varargin_1->size[0];
  sz1[1] = 1U;
  for (j = 0; j < 2; j++) {
    ysize[j] = sz1[j];
  }

  j = y->size[0] * y->size[1] * y->size[2];
  y->size[0] = (int32_T)ysize[0];
  y->size[1] = (int32_T)ysize[1];
  y->size[2] = 2;
  emxEnsureCapacity((emxArray__common *)y, j, (int32_T)sizeof(real_T));
  iy = -1;
  for (j = 1; j <= varargin_1->size[0]; j++) {
    iy++;
    y->data[iy] = varargin_1->data[j - 1];
  }

  for (j = 1; j <= varargin_2->size[0]; j++) {
    iy++;
    y->data[iy] = varargin_2->data[j - 1];
  }
}

/* End of code generation (cat.c) */
