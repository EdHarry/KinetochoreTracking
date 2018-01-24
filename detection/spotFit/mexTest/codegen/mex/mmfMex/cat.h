/*
 * cat.h
 *
 * Code generation for function 'cat'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

#ifndef __CAT_H__
#define __CAT_H__
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
#include "mmfMex_types.h"

/* Function Declarations */
extern void cat(const emlrtStack *sp, const emxArray_real_T *varargin_1, const emxArray_real_T *varargin_2, emxArray_real_T *y);
extern boolean_T isconsistent(const emxArray_real_T *y, const emxArray_real_T *x);
#endif
/* End of code generation (cat.h) */
