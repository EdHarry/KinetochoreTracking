/*
 * squeeze.c
 *
 * Code generation for function 'squeeze'
 *
 * C source code generated on: Mon May 21 17:33:35 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_F.h"
#include "squeeze.h"
#include "fitNGaussians3D_mexCode_F_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 * 
 */
void squeeze(const emxArray_real_T *a, emxArray_real_T *b)
{
    int32_T k;
    int32_T sqsz[3];
    int32_T loop_ub;
    k = 3;
    while ((k > 2) && (a->size[2] == 1)) {
        k = 2;
    }
    if (k <= 2) {
        sqsz[0] = a->size[0];
        sqsz[1] = 1;
        sqsz[2] = 1;
        k = b->size[0] * b->size[1];
        b->size[0] = sqsz[0];
        b->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)b, k, (int32_T)sizeof(real_T));
        loop_ub = a->size[0] * a->size[2];
        for (k = 0; k + 1 <= loop_ub; k++) {
            b->data[k] = a->data[k];
        }
    } else {
        for (k = 0; k < 3; k++) {
            sqsz[k] = 1;
        }
        k = 0;
        if (a->size[0] != 1) {
            sqsz[0] = a->size[0];
            k = 1;
        }
        if (a->size[2] != 1) {
            sqsz[k] = a->size[2];
        }
        sqsz[2] = 1;
        k = b->size[0] * b->size[1];
        b->size[0] = sqsz[0];
        b->size[1] = sqsz[1];
        emxEnsureCapacity((emxArray__common *)b, k, (int32_T)sizeof(real_T));
        loop_ub = a->size[0] * a->size[2];
        for (k = 0; k + 1 <= loop_ub; k++) {
            b->data[k] = a->data[k];
        }
    }
}
/* End of code generation (squeeze.c) */
