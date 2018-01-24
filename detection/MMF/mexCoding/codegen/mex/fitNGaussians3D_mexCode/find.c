/*
 * find.c
 *
 * Code generation for function 'find'
 *
 * C source code generated on: Sun May  6 15:41:06 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "find.h"
#include "fitNGaussians3D_mexCode_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 * 
 */
void find(const emxArray_boolean_T *x, emxArray_real_T *i)
{
    int32_T nx;
    int32_T idx;
    int32_T i4;
    uint32_T ii;
    boolean_T exitg1;
    boolean_T guard1 = FALSE;
    emxArray_real_T *b_i;
    emlrtHeapReferenceStackEnterFcn();
    nx = x->size[0];
    idx = 0;
    i4 = i->size[0];
    i->size[0] = nx;
    emxEnsureCapacity((emxArray__common *)i, i4, (int32_T)sizeof(real_T));
    ii = 1U;
    exitg1 = 0U;
    while ((exitg1 == 0U) && (ii <= (uint32_T)nx)) {
        guard1 = FALSE;
        if (x->data[(int32_T)ii - 1]) {
            idx++;
            i->data[idx - 1] = (real_T)ii;
            if (idx >= nx) {
                exitg1 = 1U;
            } else {
                guard1 = TRUE;
            }
        } else {
            guard1 = TRUE;
        }
        if (guard1 == TRUE) {
            ii++;
        }
    }
    if (nx == 1) {
        if (idx == 0) {
            i4 = i->size[0];
            i->size[0] = 0;
            emxEnsureCapacity((emxArray__common *)i, i4, (int32_T)sizeof(real_T));
        }
    } else {
        if (1 > idx) {
            idx = 0;
        }
        c_emxInit_real_T(&b_i, 1, TRUE);
        i4 = b_i->size[0];
        b_i->size[0] = idx;
        emxEnsureCapacity((emxArray__common *)b_i, i4, (int32_T)sizeof(real_T));
        nx = idx - 1;
        for (i4 = 0; i4 <= nx; i4++) {
            b_i->data[i4] = i->data[i4];
        }
        i4 = i->size[0];
        i->size[0] = b_i->size[0];
        emxEnsureCapacity((emxArray__common *)i, i4, (int32_T)sizeof(real_T));
        nx = b_i->size[0] - 1;
        for (i4 = 0; i4 <= nx; i4++) {
            i->data[i4] = b_i->data[i4];
        }
        emxFree_real_T(&b_i);
    }
    emlrtHeapReferenceStackLeaveFcn();
}
/* End of code generation (find.c) */
