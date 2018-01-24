/*
 * repmat.h
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Tue Nov 19 11:16:20 2013
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
#include "fitNGaussians3D_mexCode_types.h"

/* Function Declarations */
extern void b_repmat(const emxArray_real_T *a, real_T n, emxArray_real_T *b);
extern void c_repmat(real_T a, real_T m, emxArray_real_T *b);
extern void d_repmat(const emxArray_real_T *a, real_T m, emxArray_real_T *b);
extern void repmat(real_T a, const real_T m[2], emxArray_real_T *b);
#endif
/* End of code generation (repmat.h) */
