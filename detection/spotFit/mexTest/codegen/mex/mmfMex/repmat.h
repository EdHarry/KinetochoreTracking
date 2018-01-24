/*
 * repmat.h
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
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
#include "blas.h"
#include "rtwtypes.h"
#include "mmfMex_types.h"

/* Function Declarations */
extern void eml_assert_valid_size_arg(const emlrtStack *sp, const real_T varargin_1[2]);
extern void repmat(const emlrtStack *sp, const emxArray_real_T *a, const real_T m[2], emxArray_real_T *b);
#endif
/* End of code generation (repmat.h) */
