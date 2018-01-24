/*
 * length.c
 *
 * Code generation for function 'length'
 *
 * C source code generated on: Mon May 21 17:33:35 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_F.h"
#include "length.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 * 
 */
real_T length(const emxArray_real_T *x)
{
    real_T n;
    if (x->size[0] == 0) {
        n = 0.0;
    } else if (x->size[0] > 1) {
        n = (real_T)x->size[0];
    } else {
        n = 1.0;
    }
    return n;
}
/* End of code generation (length.c) */
