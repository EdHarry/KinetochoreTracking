/*
 * GaussListND_mexCode_emxutil.h
 *
 * Code generation for function 'GaussListND_mexCode_emxutil'
 *
 * C source code generated on: Sun Dec  2 22:57:02 2012
 *
 */

#ifndef __GAUSSLISTND_MEXCODE_EMXUTIL_H__
#define __GAUSSLISTND_MEXCODE_EMXUTIL_H__
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
extern void b_emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
extern void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel, int32_T elementSize);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
#endif
/* End of code generation (GaussListND_mexCode_emxutil.h) */
