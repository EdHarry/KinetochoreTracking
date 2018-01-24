/*
 * bsxfun.h
 *
 * Code generation for function 'bsxfun'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

#ifndef __BSXFUN_H__
#define __BSXFUN_H__
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
extern void b_bsxfun(const emlrtStack *sp, const emxArray_real_T *a, const real_T b_data[3], const int32_T b_size[2], emxArray_real_T *c);
extern void bsxfun(const emlrtStack *sp, const emxArray_real_T *a, const real_T b_data[3], const int32_T b_size[2], emxArray_real_T *c);
#endif
/* End of code generation (bsxfun.h) */
