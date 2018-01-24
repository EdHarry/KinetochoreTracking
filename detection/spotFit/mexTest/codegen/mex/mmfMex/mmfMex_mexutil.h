/*
 * mmfMex_mexutil.h
 *
 * Code generation for function 'mmfMex_mexutil'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

#ifndef __MMFMEX_MEXUTIL_H__
#define __MMFMEX_MEXUTIL_H__
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
extern void b_error(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location);
extern const mxArray *b_message(const emlrtStack *sp, const mxArray *b, emlrtMCInfo *location);
extern const mxArray *emlrt_marshallOut(real_T u);
extern const mxArray *message(const emlrtStack *sp, const mxArray *b, const mxArray *c, const mxArray *d, emlrtMCInfo *location);
#endif
/* End of code generation (mmfMex_mexutil.h) */
