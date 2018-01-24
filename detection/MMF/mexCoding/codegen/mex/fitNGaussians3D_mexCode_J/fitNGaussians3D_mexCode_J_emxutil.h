/*
 * fitNGaussians3D_mexCode_J_emxutil.h
 *
 * Code generation for function 'fitNGaussians3D_mexCode_J_emxutil'
 *
 * C source code generated on: Mon May 21 19:42:21 2012
 *
 */

#ifndef __FITNGAUSSIANS3D_MEXCODE_J_EMXUTIL_H__
#define __FITNGAUSSIANS3D_MEXCODE_J_EMXUTIL_H__
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
extern void b_emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
extern void c_emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
extern void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel, int32_T elementSize);
extern void emxFree_boolean_T(emxArray_boolean_T **pEmxArray);
extern void emxFree_int32_T(emxArray_int32_T **pEmxArray);
extern void emxFree_real_T(emxArray_real_T **pEmxArray);
extern void emxInit_boolean_T(emxArray_boolean_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
extern void emxInit_int32_T(emxArray_int32_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
extern void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush);
#endif
/* End of code generation (fitNGaussians3D_mexCode_J_emxutil.h) */
