/*
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Tue Nov 19 11:16:20 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "repmat.h"
#include "fitNGaussians3D_mexCode_emxutil.h"

/* Function Definitions */
void b_repmat(const emxArray_real_T *a, real_T n, emxArray_real_T *b)
{
  int32_T mv[2];
  int32_T b_a[2];
  int32_T iv0[2];
  int32_T ib;
  int32_T jtilecol;
  int32_T ia;
  int32_T k;
  mv[0] = 1;
  mv[1] = (int32_T)n;
  b_a[0] = a->size[0];
  b_a[1] = 1;
  for (ib = 0; ib < 2; ib++) {
    iv0[ib] = b_a[ib] * mv[ib];
  }

  ib = b->size[0] * b->size[1];
  b->size[0] = iv0[0];
  b->size[1] = iv0[1];
  emxEnsureCapacity((emxArray__common *)b, ib, (int32_T)sizeof(real_T));
  if ((iv0[0] == 0) || (iv0[1] == 0)) {
  } else {
    ib = 0;
    for (jtilecol = 1; jtilecol <= (int32_T)n; jtilecol++) {
      ia = 0;
      for (k = 1; k <= a->size[0]; k++) {
        b->data[ib] = a->data[ia];
        ia++;
        ib++;
      }
    }
  }
}

void c_repmat(real_T a, real_T m, emxArray_real_T *b)
{
  int32_T b_m[2];
  int32_T outsize[2];
  int32_T i3;
  int32_T loop_ub;
  b_m[0] = (int32_T)m;
  b_m[1] = 1;
  for (i3 = 0; i3 < 2; i3++) {
    outsize[i3] = b_m[i3];
  }

  i3 = b->size[0];
  b->size[0] = outsize[0];
  emxEnsureCapacity((emxArray__common *)b, i3, (int32_T)sizeof(real_T));
  loop_ub = outsize[0];
  for (i3 = 0; i3 < loop_ub; i3++) {
    b->data[i3] = a;
  }
}

void d_repmat(const emxArray_real_T *a, real_T m, emxArray_real_T *b)
{
  int32_T mv[2];
  int32_T iv1[2];
  int32_T ia;
  int32_T ib;
  int32_T iacol;
  int32_T jcol;
  int32_T itilerow;
  mv[0] = (int32_T)m;
  mv[1] = 1;
  for (ia = 0; ia < 2; ia++) {
    iv1[ia] = a->size[ia] * mv[ia];
  }

  ia = b->size[0] * b->size[1];
  b->size[0] = iv1[0];
  b->size[1] = iv1[1];
  emxEnsureCapacity((emxArray__common *)b, ia, (int32_T)sizeof(real_T));
  if ((iv1[0] == 0) || (iv1[1] == 0)) {
  } else {
    ia = 1;
    ib = 0;
    iacol = 1;
    for (jcol = 1; jcol <= a->size[1]; jcol++) {
      for (itilerow = 1; itilerow <= (int32_T)m; itilerow++) {
        b->data[ib] = a->data[iacol - 1];
        ia = iacol + 1;
        ib++;
      }

      iacol = ia;
    }
  }
}

void repmat(real_T a, const real_T m[2], emxArray_real_T *b)
{
  int32_T outsize[2];
  int32_T i1;
  int32_T loop_ub;
  for (i1 = 0; i1 < 2; i1++) {
    outsize[i1] = (int32_T)m[i1];
  }

  i1 = b->size[0] * b->size[1];
  b->size[0] = outsize[0];
  emxEnsureCapacity((emxArray__common *)b, i1, (int32_T)sizeof(real_T));
  i1 = b->size[0] * b->size[1];
  b->size[1] = outsize[1];
  emxEnsureCapacity((emxArray__common *)b, i1, (int32_T)sizeof(real_T));
  loop_ub = outsize[0] * outsize[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    b->data[i1] = a;
  }
}

/* End of code generation (repmat.c) */
