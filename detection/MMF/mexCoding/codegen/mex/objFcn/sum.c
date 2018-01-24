/*
 * sum.c
 *
 * Code generation for function 'sum'
 *
 * C source code generated on: Fri May 25 21:48:51 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "objFcn.h"
#include "sum.h"
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
void sum(const emxArray_real_T *x, emxArray_real_T *y)
{
    int32_T vlen;
    uint32_T sz[2];
    int32_T npages;
    int32_T ixstart;
    int32_T iy;
    int32_T i;
    real_T s;
    int32_T k;
    for (vlen = 0; vlen < 2; vlen++) {
        sz[vlen] = (uint32_T)x->size[vlen];
    }
    sz[0] = 1U;
    vlen = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = (int32_T)sz[1];
    emxEnsureCapacity((emxArray__common *)y, vlen, (int32_T)sizeof(real_T));
    if ((x->size[0] == 0) || (x->size[1] == 0)) {
        vlen = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, vlen, (int32_T)sizeof(real_T));
        npages = y->size[1] - 1;
        for (vlen = 0; vlen <= npages; vlen++) {
            y->data[y->size[0] * vlen] = 0.0;
        }
    } else {
        vlen = x->size[0];
        npages = x->size[1];
        ixstart = -1;
        iy = -1;
        for (i = 1; i <= npages; i++) {
            ixstart++;
            s = x->data[ixstart];
            for (k = 2; k <= vlen; k++) {
                ixstart++;
                s += x->data[ixstart];
            }
            iy++;
            y->data[iy] = s;
        }
    }
}
/* End of code generation (sum.c) */
