/*
 * repmat.h
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Sun Dec  2 22:57:02 2012
 *
 */

#ifndef __REPMAT_H__
#define __REPMAT_H__
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
extern void b_repmat(const emxArray_real_T *a, const real_T m[2], emxArray_real_T *b);
extern void repmat(const emxArray_real_T *a, const real_T m[2], emxArray_real_T *b);
#endif
/* End of code generation (repmat.h) */
