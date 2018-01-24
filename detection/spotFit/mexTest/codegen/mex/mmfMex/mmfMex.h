/*
 * mmfMex.h
 *
 * Code generation for function 'mmfMex'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

#ifndef __MMFMEX_H__
#define __MMFMEX_H__
/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blas.h"
#include "rtwtypes.h"
#include "mmfMex_types.h"

/* Function Declarations */
extern void mmfMex(const emlrtStack *sp, const emxArray_real_T *maskAmp, const emxArray_real_T *maskCoord, const emxArray_real_T *X, const real_T sigma_data[3], const int32_T sigma_size[2], emxArray_real_T *jacobian, emxArray_real_T *resi, emxArray_real_T *amp, real_T *bg, real_T *nGauss, real_T *nDim);
#endif
/* End of code generation (mmfMex.h) */
