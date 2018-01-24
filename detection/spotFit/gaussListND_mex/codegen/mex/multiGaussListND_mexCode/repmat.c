/*
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Sun Dec  2 23:59:22 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "multiGaussListND_mexCode.h"
#include "repmat.h"
#include "multiGaussListND_mexCode_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void repmat(const real_T a_data[3], const int32_T a_size[2], const real_T m[2],
            emxArray_real_T *b)
{
  int32_T outsize[2];
  int32_T ia;
  int32_T mv[2];
  int32_T ib;
  int32_T iacol;
  int32_T jcol;
  int32_T itilerow;
  for (ia = 0; ia < 2; ia++) {
    outsize[ia] = a_size[ia] * (int32_T)m[ia];
    mv[ia] = (int32_T)m[ia];
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
    for (jcol = 1; jcol <= a_size[1]; jcol++) {
      for (itilerow = 1; itilerow <= mv[0]; itilerow++) {
        b->data[ib] = a_data[iacol - 1];
        ia = iacol + 1;
        ib++;
      }

      iacol = ia;
    }
  }
}

/* End of code generation (repmat.c) */
