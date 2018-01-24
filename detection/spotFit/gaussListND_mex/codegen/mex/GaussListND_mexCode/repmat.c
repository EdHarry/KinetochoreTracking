/*
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Sun Dec  2 22:57:02 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GaussListND_mexCode.h"
#include "repmat.h"
#include "GaussListND_mexCode_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_repmat(const emxArray_real_T *a, const real_T m[2], emxArray_real_T *b)
{
  int32_T mv[2];
  int32_T ia;
  int32_T outsize[2];
  int32_T ib;
  int32_T iacol;
  int32_T jcol;
  int32_T itilerow;
  int32_T k;
  for (ia = 0; ia < 2; ia++) {
    mv[ia] = (int32_T)m[ia];
  }

  for (ia = 0; ia < 2; ia++) {
    outsize[ia] = a->size[ia] * mv[ia];
  }

  ia = b->size[0] * b->size[1];
  b->size[0] = outsize[0];
  b->size[1] = outsize[1];
  emxEnsureCapacity((emxArray__common *)b, ia, (int32_T)sizeof(real_T));
  if ((outsize[0] == 0) || (outsize[1] == 0)) {
  } else {
    ia = 1;
    ib = 0;
    iacol = 1;
    for (jcol = 1; jcol <= a->size[1]; jcol++) {
      for (itilerow = 1; itilerow <= mv[0]; itilerow++) {
        ia = iacol;
        for (k = 1; k <= a->size[0]; k++) {
          b->data[ib] = a->data[ia - 1];
          ia++;
          ib++;
        }
      }

      iacol = ia;
    }
  }
}

void repmat(const emxArray_real_T *a, const real_T m[2], emxArray_real_T *b)
{
  int32_T mv[2];
  int32_T ia;
  int32_T outsize[2];
  int32_T ib;
  int32_T iacol;
  int32_T jcol;
  int32_T itilerow;
  for (ia = 0; ia < 2; ia++) {
    mv[ia] = (int32_T)m[ia];
  }

  for (ia = 0; ia < 2; ia++) {
    outsize[ia] = a->size[ia] * mv[ia];
  }

  ia = b->size[0] * b->size[1];
  b->size[0] = outsize[0];
  b->size[1] = outsize[1];
  emxEnsureCapacity((emxArray__common *)b, ia, (int32_T)sizeof(real_T));
  if ((outsize[0] == 0) || (outsize[1] == 0)) {
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
