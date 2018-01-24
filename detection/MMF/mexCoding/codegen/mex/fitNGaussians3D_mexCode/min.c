/*
 * min.c
 *
 * Code generation for function 'min'
 *
 * C source code generated on: Sun May  6 15:41:06 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "min.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 * 
 */
real_T b_min(const emxArray_real_T *varargin_1)
{
    real_T minval;
    int32_T n;
    int32_T itmp;
    int32_T ix;
    boolean_T guard1 = FALSE;
    boolean_T searchingForNonNaN;
    int32_T k;
    boolean_T exitg1;
    n = varargin_1->size[0];
    minval = varargin_1->data[0];
    itmp = 1;
    if (n == 1) {
    } else {
        ix = 0;
        guard1 = FALSE;
        if (muDoubleScalarIsNaN(varargin_1->data[0])) {
            searchingForNonNaN = TRUE;
            k = 2;
            exitg1 = 0U;
            while ((exitg1 == 0U) && (k <= n)) {
                ix++;
                if (!muDoubleScalarIsNaN(varargin_1->data[ix])) {
                    minval = varargin_1->data[ix];
                    itmp = k;
                    searchingForNonNaN = FALSE;
                    exitg1 = 1U;
                } else {
                    k++;
                }
            }
            if (searchingForNonNaN) {
            } else {
                guard1 = TRUE;
            }
        } else {
            guard1 = TRUE;
        }
        if (guard1 == TRUE) {
            while (itmp + 1 <= n) {
                if (varargin_1->data[itmp] < minval) {
                    minval = varargin_1->data[itmp];
                }
                itmp++;
            }
        }
    }
    return minval;
}
/* End of code generation (min.c) */
