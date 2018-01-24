/*
 * mrdivide.c
 *
 * Code generation for function 'mrdivide'
 *
 * C source code generated on: Fri May 25 21:48:51 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "objFcn.h"
#include "mrdivide.h"
#include "repmat.h"
#include "objFcn_mexutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 * 
 */
int32_T mrdivide(int32_T A, real_T B)
{
    real_T d0;
    d0 = (real_T)A / B;
    if ((d0 < 4.503599627370496E+15) && (d0 > -4.503599627370496E+15)) {
        d0 = d0 < 0.0 ? muDoubleScalarCeil(d0 - 0.5) : muDoubleScalarFloor(d0 + 0.5);
    }
    return _s32_d_(d0);
}
/* End of code generation (mrdivide.c) */
