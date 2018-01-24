/*
 * fitNGaussians3D_mexCode_emxAPI.h
 *
 * Code generation for function 'fitNGaussians3D_mexCode_emxAPI'
 *
 * C source code generated on: Thu May  3 13:06:48 2012
 *
 */

#ifndef __FITNGAUSSIANS3D_MEXCODE_EMXAPI_H__
#define __FITNGAUSSIANS3D_MEXCODE_EMXAPI_H__
/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"

#include "rtwtypes.h"
#include "fitNGaussians3D_mexCode_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern emxArray_real_T *emxCreateND_real_T(int32_T numDimensions, int32_T *size);
extern emxArray_real_T *emxCreateWrapperND_real_T(real_T *data, int32_T numDimensions, int32_T *size);
extern emxArray_real_T *emxCreateWrapper_real_T(real_T *data, int32_T rows, int32_T cols);
extern emxArray_real_T *emxCreate_real_T(int32_T rows, int32_T cols);
extern void emxDestroyArray_real_T(emxArray_real_T *emxArray);
#endif
/* End of code generation (fitNGaussians3D_mexCode_emxAPI.h) */
