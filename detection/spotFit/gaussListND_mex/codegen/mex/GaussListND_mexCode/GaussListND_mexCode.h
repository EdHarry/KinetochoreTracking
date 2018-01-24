/*
 * GaussListND_mexCode.h
 *
 * Code generation for function 'GaussListND_mexCode'
 *
 * C source code generated on: Sun Dec  2 22:57:02 2012
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
#include "GaussListND_mexCode_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void GaussListND_mexCode(const emxArray_real_T *coordList, const emxArray_real_T *sigma, const emxArray_real_T *center, emxArray_real_T *gaussList);
#endif
/* End of code generation (GaussListND_mexCode.h) */
