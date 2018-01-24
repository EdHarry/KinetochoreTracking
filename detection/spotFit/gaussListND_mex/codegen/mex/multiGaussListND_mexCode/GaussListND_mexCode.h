/*
 * GaussListND_mexCode.h
 *
 * Code generation for function 'GaussListND_mexCode'
 *
 * C source code generated on: Sun Dec  2 23:59:22 2012
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
#include "multiGaussListND_mexCode_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void GaussListND_mexCode(const emxArray_real_T *coordList, const real_T sigma_data[3], const int32_T sigma_size[2], const real_T center_data[3], const int32_T center_size[2], emxArray_real_T *gaussList);
#endif
/* End of code generation (GaussListND_mexCode.h) */
