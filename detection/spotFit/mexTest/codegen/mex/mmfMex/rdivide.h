/*
 * rdivide.h
 *
 * Code generation for function 'rdivide'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

#ifndef __RDIVIDE_H__
#define __RDIVIDE_H__
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
extern void b_rdivide(const emlrtStack *sp, const emxArray_real_T *x, const emxArray_real_T *y, emxArray_real_T *z);
extern void rdivide(const emlrtStack *sp, const emxArray_real_T *y, emxArray_real_T *z);
#endif
/* End of code generation (rdivide.h) */
