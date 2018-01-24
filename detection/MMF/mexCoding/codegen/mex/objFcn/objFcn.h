/*
 * objFcn.h
 *
 * Code generation for function 'objFcn'
 *
 * C source code generated on: Fri May 25 21:48:51 2012
 *
 */

#ifndef __OBJFCN_H__
#define __OBJFCN_H__
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
extern void objFcn(const int32_T *mode, int32_T m, int32_T n, int32_T ldfj, int32_T needfi, emxArray_real_T *x, emxArray_real_T *fjac, int32_T nstate, const struct_T *user, emxArray_real_T *f);
#endif
/* End of code generation (objFcn.h) */
