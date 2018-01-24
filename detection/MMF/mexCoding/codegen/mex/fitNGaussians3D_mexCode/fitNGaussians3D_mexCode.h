/*
 * fitNGaussians3D_mexCode.h
 *
 * Code generation for function 'fitNGaussians3D_mexCode'
 *
 * C source code generated on: Tue Nov 19 11:16:19 2013
 *
 */

#ifndef __FITNGAUSSIANS3D_MEXCODE_H__
#define __FITNGAUSSIANS3D_MEXCODE_H__
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
#include "fitNGaussians3D_mexCode_types.h"

/* Function Declarations */
extern void fitNGaussians3D_mexCode(emxArray_real_T *x0, const emxArray_real_T *image, const emxArray_real_T *b_index, const real_T psfSigma[2], emxArray_real_T *F, emxArray_real_T *J);
#endif
/* End of code generation (fitNGaussians3D_mexCode.h) */
