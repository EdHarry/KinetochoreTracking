/*
 * power.c
 *
 * Code generation for function 'power'
 *
 * C source code generated on: Fri May 25 21:48:51 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "objFcn.h"
#include "power.h"
#include "objFcn_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 * 
 */
void power(const emxArray_real_T *a, emxArray_real_T *y)
{
    uint32_T unnamed_idx_0;
    int32_T loop_ub;
    unnamed_idx_0 = (uint32_T)a->size[0];
    loop_ub = y->size[0];
    y->size[0] = (int32_T)unnamed_idx_0;
    emxEnsureCapacity((emxArray__common *)y, loop_ub, (int32_T)sizeof(real_T));
    loop_ub = y->size[0];
    for (unnamed_idx_0 = 1U; unnamed_idx_0 <= (uint32_T)loop_ub; unnamed_idx_0++) {
        y->data[(int32_T)unnamed_idx_0 - 1] = muDoubleScalarPower(a->data[(int32_T)unnamed_idx_0 - 1], 2.0);
    }
}
/* End of code generation (power.c) */
