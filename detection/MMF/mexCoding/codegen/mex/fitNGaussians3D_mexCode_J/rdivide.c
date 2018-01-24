/*
 * rdivide.c
 *
 * Code generation for function 'rdivide'
 *
 * C source code generated on: Mon May 21 19:42:21 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_J.h"
#include "rdivide.h"
#include "fitNGaussians3D_mexCode_J_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 * 
 */
void rdivide(const emxArray_real_T *x, real_T y, emxArray_real_T *z)
{
    int32_T i2;
    int32_T loop_ub;
    i2 = z->size[0];
    z->size[0] = x->size[0];
    emxEnsureCapacity((emxArray__common *)z, i2, (int32_T)sizeof(real_T));
    loop_ub = x->size[0] - 1;
    for (i2 = 0; i2 <= loop_ub; i2++) {
        z->data[i2] = x->data[i2] / y;
    }
}
/* End of code generation (rdivide.c) */
