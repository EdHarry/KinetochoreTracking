/*
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Mon May 28 14:49:15 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_J.h"
#include "repmat.h"
#include "fitNGaussians3D_mexCode_J_emxutil.h"
#include "fitNGaussians3D_mexCode_J_rtwutil.h"

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
  int32_T ia;
  int32_T mv[2];
  int32_T outsize[2];
  int32_T ib;
  int32_T iacol;
  int32_T jcol;
  int32_T itilerow;
  d0 = rt_roundd_snf(m);
  if (d0 < 2.147483648E+9) {
    if (d0 >= -2.147483648E+9) {
      ia = (int32_T)d0;
    } else {
      ia = MIN_int32_T;
    }
  } else if (d0 >= 2.147483648E+9) {
    ia = MAX_int32_T;
  } else {
    ia = 0;
  }

  mv[0] = ia;
  mv[1] = 1;
  for (ia = 0; ia < 2; ia++) {
    outsize[ia] = a->size[ia] * mv[ia];
  }

  ia = b->size[0] * b->size[1];
  b->size[0] = outsize[0];
  b->size[1] = outsize[1];
  emxEnsureCapacity((emxArray__common *)b, ia, (int32_T)sizeof(real_T));
  if ((b->size[0] == 0) || (b->size[1] == 0)) {
  } else {
    ia = 1;
    ib = 0;
    iacol = 1;
    for (jcol = 1; jcol <= a->size[1]; jcol++) {
      for (itilerow = 1; itilerow <= mv[0]; itilerow++) {
        b->data[ib] = a->data[iacol - 1];
        ia = iacol + 1;
        ib++;
      }

      iacol = ia;
    }
  }
}

/* End of code generation (repmat.c) */
