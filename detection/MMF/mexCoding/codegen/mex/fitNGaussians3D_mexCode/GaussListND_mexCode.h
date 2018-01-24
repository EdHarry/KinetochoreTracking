/*
 * GaussListND_mexCode.h
 *
 * Code generation for function 'GaussListND_mexCode'
 *
 * C source code generated on: Tue Nov 19 11:16:20 2013
 *
 */

#ifndef __GAUSSLISTND_MEXCODE_H__
#define __GAUSSLISTND_MEXCODE_H__
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
extern void GaussListND_mexCode(const emxArray_real_T *coordList, real_T sigma, real_T center, emxArray_real_T *gaussList);
#endif
/* End of code generation (GaussListND_mexCode.h) */
