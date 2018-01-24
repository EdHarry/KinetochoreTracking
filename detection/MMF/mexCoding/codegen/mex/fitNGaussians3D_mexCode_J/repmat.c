/*
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Mon May 21 19:42:21 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_J.h"
#include "repmat.h"
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
void repmat(const emxArray_real_T *a, real_T m, emxArray_real_T *b)
{
    real_T d0;
    int32_T ncols;
    int32_T mv[2];
    int32_T iv1[2];
    int32_T ia;
    int32_T ib;
    int32_T iacol;
    int32_T jcol;
    int32_T itilerow;
    d0 = m;
    d0 = d0 < 0.0 ? muDoubleScalarCeil(d0 - 0.5) : muDoubleScalarFloor(d0 + 0.5);
    if (d0 < 2.147483648E+9) {
        if (d0 >= -2.147483648E+9) {
            ncols = (int32_T)d0;
        } else {
            ncols = MIN_int32_T;
        }
    } else if (d0 >= 2.147483648E+9) {
        ncols = MAX_int32_T;
    } else {
        ncols = 0;
    }
    mv[0] = ncols;
    mv[1] = 1;
    for (ncols = 0; ncols < 2; ncols++) {
        iv1[ncols] = a->size[ncols] * mv[ncols];
    }
    ncols = b->size[0] * b->size[1];
    b->size[0] = iv1[0];
    b->size[1] = iv1[1];
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
/* End of code generation (repmat.c) */
