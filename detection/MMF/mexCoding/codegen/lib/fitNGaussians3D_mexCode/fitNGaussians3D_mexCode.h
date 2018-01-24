/*
 * fitNGaussians3D_mexCode.h
 *
 * Code generation for function 'fitNGaussians3D_mexCode'
 *
 * C source code generated on: Thu May  3 12:56:20 2012
 *
 */

#ifndef __FITNGAUSSIANS3D_MEXCODE_H__
#define __FITNGAUSSIANS3D_MEXCODE_H__
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
extern void fitNGaussians3D_mexCode(emxArray_real_T *x0, const emxArray_real_T *image, const emxArray_real_T *b_index, const real_T psfSigma[2], emxArray_real_T *F, emxArray_real_T *J);
#endif
/* End of code generation (fitNGaussians3D_mexCode.h) */
