/*
 * GaussListND_mexCode.h
 *
 * Code generation for function 'GaussListND_mexCode'
 *
 * C source code generated on: Mon May 21 17:33:35 2012
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
#include "blascompat32.h"
#include "rtwtypes.h"
#include "fitNGaussians3D_mexCode_F_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void GaussListND_mexCode(const emxArray_real_T *coordList, real_T sigma, real_T center, emxArray_real_T *gaussList);
#endif
/* End of code generation (GaussListND_mexCode.h) */
