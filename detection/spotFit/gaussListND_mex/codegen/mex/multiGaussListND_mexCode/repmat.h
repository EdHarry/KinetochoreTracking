/*
 * repmat.h
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Sun Dec  2 23:59:22 2012
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
#include "multiGaussListND_mexCode_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void repmat(const real_T a_data[3], const int32_T a_size[2], const real_T m[2], emxArray_real_T *b);
#endif
/* End of code generation (repmat.h) */
