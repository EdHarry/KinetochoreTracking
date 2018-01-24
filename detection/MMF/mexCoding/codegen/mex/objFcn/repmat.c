/*
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Fri May 25 21:48:51 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "objFcn.h"
#include "repmat.h"
#include "objFcn_emxutil.h"
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
void b_repmat(const emxArray_real_T *a, int32_T m, emxArray_real_T *b)
{
    int32_T mv[2];
    int32_T ncols;
    int32_T iv0[2];
    int32_T ia;
    int32_T ib;
    int32_T iacol;
    int32_T jcol;
    int32_T itilerow;
    mv[0] = m;
    mv[1] = 1;
    for (ncols = 0; ncols < 2; ncols++) {
        iv0[ncols] = a->size[ncols] * mv[ncols];
    }
    ncols = b->size[0] * b->size[1];
    b->size[0] = iv0[0];
    b->size[1] = iv0[1];
    emxEnsureCapacity((emxArray__common *)b, ncols, (int32_T)sizeof(real_T));
    if ((b->size[0] == 0) || (b->size[1] == 0)) {
    } else {
        ncols = a->size[1];
        ia = -1;
        ib = 0;
        iacol = 0;
        for (jcol = 1; jcol <= ncols; jcol++) {
            for (itilerow = 1; itilerow <= mv[0]; itilerow++) {
                b->data[ib] = a->data[iacol];
                ia = iacol;
                ib++;
            }
            iacol = ia + 1;
        }
    }
}

/*
 * 
 */
void c_repmat(const emxArray_real_T *a, int32_T n, emxArray_real_T *b)
{
    int32_T mv[2];
    int32_T b_a[2];
    int32_T nrows;
    int32_T iv1[2];
    int32_T ib;
    int32_T jtilecol;
    int32_T ia;
    int32_T k;
    mv[0] = 1;
    mv[1] = n;
    b_a[0] = a->size[0];
    b_a[1] = 1;
    for (nrows = 0; nrows < 2; nrows++) {
        iv1[nrows] = b_a[nrows] * mv[nrows];
    }
    nrows = b->size[0] * b->size[1];
    b->size[0] = iv1[0];
    b->size[1] = iv1[1];
    emxEnsureCapacity((emxArray__common *)b, nrows, (int32_T)sizeof(real_T));
    if ((b->size[0] == 0) || (b->size[1] == 0)) {
    } else {
        nrows = a->size[0];
        ib = 0;
        for (jtilecol = 1; jtilecol <= mv[1]; jtilecol++) {
            ia = 0;
            for (k = 1; k <= nrows; k++) {
                b->data[ib] = a->data[ia];
                ia++;
                ib++;
            }
        }
    }
}

/*
 * 
 */
void d_repmat(real_T a, int32_T m, emxArray_real_T *b)
{
    int32_T b_m[2];
    int32_T i2;
    int32_T outsize[2];
    int32_T loop_ub;
    b_m[0] = m;
    b_m[1] = 1;
    for (i2 = 0; i2 < 2; i2++) {
        outsize[i2] = b_m[i2];
    }
    i2 = b->size[0];
    b->size[0] = outsize[0];
    emxEnsureCapacity((emxArray__common *)b, i2, (int32_T)sizeof(real_T));
    loop_ub = outsize[0] - 1;
    for (i2 = 0; i2 <= loop_ub; i2++) {
        b->data[i2] = a;
    }
}

/*
 * 
 */
void repmat(real_T a, const real_T m[2], emxArray_real_T *b)
{
    int32_T i0;
    real_T d1;
    int32_T outsize[2];
    int32_T loop_ub;
    for (i0 = 0; i0 < 2; i0++) {
        d1 = m[i0];
        if ((d1 < 4.503599627370496E+15) && (d1 > -4.503599627370496E+15)) {
            d1 = d1 < 0.0 ? muDoubleScalarCeil(d1 - 0.5) : muDoubleScalarFloor(d1 + 0.5);
        }
        outsize[i0] = _s32_d_(d1);
    }
    i0 = b->size[0] * b->size[1];
    b->size[0] = outsize[0];
    emxEnsureCapacity((emxArray__common *)b, i0, (int32_T)sizeof(real_T));
    i0 = b->size[0] * b->size[1];
    b->size[1] = outsize[1];
    emxEnsureCapacity((emxArray__common *)b, i0, (int32_T)sizeof(real_T));
    loop_ub = outsize[0] * outsize[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        b->data[i0] = a;
    }
}
/* End of code generation (repmat.c) */
