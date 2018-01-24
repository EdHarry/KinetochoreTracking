/*
 * multiGaussListND_mexCode.c
 *
 * Code generation for function 'multiGaussListND_mexCode'
 *
 * C source code generated on: Sun Dec  2 23:59:22 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "multiGaussListND_mexCode.h"
#include "GaussListND_mexCode.h"
#include "multiGaussListND_mexCode_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void multiGaussListND_mexCode(const emxArray_real_T *maskCoord, const real_T
  X_data[1500], const int32_T X_size[2], const real_T sigma_data[3], const
  int32_T sigma_size[2], emxArray_real_T *gaussList)
{
  int32_T maskCoord_idx_0;
  int32_T i0;
  int32_T gaussIdx;
  emxArray_real_T *tmp;
  real_T b_X_data[3];
  int32_T b_X_size[2];
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);

  /* function [ gaussList, newAmpBg ] = multiGaussListND_mexCode( maskAmp, maskCoord, X, sigma ) */
  maskCoord_idx_0 = maskCoord->size[0];
  i0 = gaussList->size[0] * gaussList->size[1];
  gaussList->size[0] = maskCoord_idx_0;
  emxEnsureCapacity((emxArray__common *)gaussList, i0, (int32_T)sizeof(real_T));
  i0 = gaussList->size[0] * gaussList->size[1];
  gaussList->size[1] = X_size[0] + 1;
  emxEnsureCapacity((emxArray__common *)gaussList, i0, (int32_T)sizeof(real_T));
  maskCoord_idx_0 = maskCoord->size[0] * (X_size[0] + 1);
  for (i0 = 0; i0 < maskCoord_idx_0; i0++) {
    gaussList->data[i0] = 1.0;
  }

  /* % GAUSS LIST CALCULATION */
  gaussIdx = 0;
  emxInit_real_T(&tmp, 3, TRUE);
  while (gaussIdx <= X_size[0] - 1) {
    maskCoord_idx_0 = X_size[1];
    b_X_size[0] = 1;
    b_X_size[1] = X_size[1];
    for (i0 = 0; i0 < maskCoord_idx_0; i0++) {
      b_X_data[i0] = X_data[gaussIdx + X_size[0] * i0];
    }

    GaussListND_mexCode(maskCoord, sigma_data, sigma_size, b_X_data, b_X_size,
                        tmp);
    maskCoord_idx_0 = tmp->size[0] - 1;
    for (i0 = 0; i0 <= maskCoord_idx_0; i0++) {
      gaussList->data[i0 + gaussList->size[0] * gaussIdx] = tmp->data[i0];
    }

    gaussIdx++;
  }

  emxFree_real_T(&tmp);

  /* % GAUSS FITTING */
  /* I = (G1,G2,G3,....,1) *(a1;a2;a3;...;bg); */
  /* To fit -> a1;a2..;bg = (G1,G2,G3,...,1) \ I */
  /* newAmpBg = gaussList \ maskAmp; */
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (multiGaussListND_mexCode.c) */
