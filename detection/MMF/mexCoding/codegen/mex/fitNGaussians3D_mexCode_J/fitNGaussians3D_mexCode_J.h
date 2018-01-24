/*
 * fitNGaussians3D_mexCode_J.h
 *
 * Code generation for function 'fitNGaussians3D_mexCode_J'
 *
 * C source code generated on: Mon May 21 19:42:21 2012
 *
 */

#ifndef __FITNGAUSSIANS3D_MEXCODE_J_H__
#define __FITNGAUSSIANS3D_MEXCODE_J_H__
/* Include files */
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "mwmathutil.h"

#include "tmwtypes.h"
#include "mex.h"
#include "emlrt.h"
#include "blascompat32.h"
#include "rtwtypes.h"
#include "fitNGaussians3D_mexCode_J_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void fitNGaussians3D_mexCode_J(emxArray_real_T *x0, const emxArray_real_T *image, const emxArray_real_T *b_index, const real_T psfSigma[2], emxArray_real_T *J);
#endif
/* End of code generation (fitNGaussians3D_mexCode_J.h) */
