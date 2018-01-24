/*
 * repmat.h
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Thu May  3 13:06:48 2012
 *
 */

#ifndef __REPMAT_H__
#define __REPMAT_H__
/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"

#include "rtwtypes.h"
#include "fitNGaussians3D_mexCode_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void b_repmat(real_T a, real_T m, emxArray_real_T *b);
extern void c_repmat(const emxArray_real_T *a, real_T m, emxArray_real_T *b);
extern void repmat(const emxArray_real_T *a, real_T n, emxArray_real_T *b);
#endif
/* End of code generation (repmat.h) */
