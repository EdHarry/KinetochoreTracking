/*
 * GaussListND_mexCode.h
 *
 * Code generation for function 'GaussListND_mexCode'
 *
 * C source code generated on: Fri May 25 21:48:51 2012
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
#include "objFcn_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void GaussListND_mexCode(const emxArray_real_T *coordList, real_T sigma, real_T center, emxArray_real_T *gaussList);
#endif
/* End of code generation (GaussListND_mexCode.h) */
