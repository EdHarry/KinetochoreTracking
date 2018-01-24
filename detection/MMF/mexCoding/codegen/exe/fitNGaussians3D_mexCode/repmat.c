/*
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Thu May  3 13:06:48 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "repmat.h"
#include "fitNGaussians3D_mexCode_emxutil.h"
#include "fitNGaussians3D_mexCode_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_repmat(real_T a, real_T m, emxArray_real_T *b)
{
  int32_T b_m[2];
  real_T d1;
  int32_T i4;
  int32_T outsize[2];
  int32_T loop_ub;
  d1 = rt_roundd_snf(m);
  if (d1 < 2.147483648E+9) {
    if (d1 >= -2.147483648E+9) {
      i4 = (int32_T)d1;
    } else {
      i4 = MIN_int32_T;
    }
  } else if (d1 >= 2.147483648E+9) {
    i4 = MAX_int32_T;
  } else {
    i4 = 0;
  }

  b_m[0] = i4;
  b_m[1] = 1;
  for (i4 = 0; i4 < 2; i4++) {
    outsize[i4] = b_m[i4];
  }

  i4 = b->size[0];
  b->size[0] = outsize[0];
  emxEnsureCapacity((emxArray__common *)b, i4, (int32_T)sizeof(real_T));
  loop_ub = outsize[0] - 1;
  for (i4 = 0; i4 <= loop_ub; i4++) {
    b->data[i4] = a;
  }
}

void c_repmat(const emxArray_real_T *a, real_T m, emxArray_real_T *b)
{
  real_T d2;
  int32_T ia;
  int32_T mv[2];
  int32_T outsize[2];
  int32_T ib;
  int32_T iacol;
  int32_T jcol;
  int32_T itilerow;
  d2 = rt_roundd_snf(m);
  if (d2 < 2.147483648E+9) {
    if (d2 >= -2.147483648E+9) {
      ia = (int32_T)d2;
    } else {
      ia = MIN_int32_T;
    }
  } else if (d2 >= 2.147483648E+9) {
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

void repmat(const emxArray_real_T *a, real_T n, emxArray_real_T *b)
{
  int32_T mv[2];
  real_T d0;
  int32_T ib;
  int32_T b_a[2];
  int32_T outsize[2];
  int32_T jtilecol;
  int32_T ia;
  int32_T k;
  mv[0] = 1;
  d0 = rt_roundd_snf(n);
  if (d0 < 2.147483648E+9) {
    if (d0 >= -2.147483648E+9) {
      ib = (int32_T)d0;
    } else {
      ib = MIN_int32_T;
    }
  } else if (d0 >= 2.147483648E+9) {
    ib = MAX_int32_T;
  } else {
    ib = 0;
  }

  mv[1] = ib;
  b_a[0] = a->size[0];
  b_a[1] = 1;
  for (ib = 0; ib < 2; ib++) {
    outsize[ib] = b_a[ib] * mv[ib];
  }

  ib = b->size[0] * b->size[1];
  b->size[0] = outsize[0];
  b->size[1] = outsize[1];
  emxEnsureCapacity((emxArray__common *)b, ib, (int32_T)sizeof(real_T));
  if ((b->size[0] == 0) || (b->size[1] == 0)) {
  } else {
    ib = 0;
    for (jtilecol = 1; jtilecol <= mv[1]; jtilecol++) {
      ia = 0;
      for (k = 1; k <= a->size[0]; k++) {
        b->data[ib] = a->data[ia];
        ia++;
        ib++;
      }
    }
  }
}

/* End of code generation (repmat.c) */
