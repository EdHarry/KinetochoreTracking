/*
 * repmat.c
 *
 * Code generation for function 'repmat'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "repmat.h"
#include "mmfMex_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "mmfMex_mexutil.h"
#include "mmfMex_data.h"

/* Variable Definitions */
static emlrtRSInfo u_emlrtRSI = { 60, "repmat",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/repmat.m" };

static emlrtRSInfo y_emlrtRSI = { 41, "eml_assert_valid_size_arg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtRSInfo ab_emlrtRSI = { 56, "eml_assert_valid_size_arg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtMCInfo emlrtMCI = { 42, 9, "eml_assert_valid_size_arg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtMCInfo b_emlrtMCI = { 41, 19, "eml_assert_valid_size_arg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtMCInfo c_emlrtMCI = { 57, 5, "eml_assert_valid_size_arg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtMCInfo d_emlrtMCI = { 56, 15, "eml_assert_valid_size_arg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtRTEInfo k_emlrtRTEI = { 45, 1, "repmat",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/repmat.m" };

static emlrtRSInfo mq_emlrtRSI = { 57, "eml_assert_valid_size_arg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

static emlrtRSInfo pq_emlrtRSI = { 42, "eml_assert_valid_size_arg",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_assert_valid_size_arg.m"
};

/* Function Definitions */
void eml_assert_valid_size_arg(const emlrtStack *sp, const real_T varargin_1[2])
{
  int32_T i;
  boolean_T guard1 = FALSE;
  int32_T exitg1;
  boolean_T p;
  const mxArray *y;
  static const int32_T iv8[2] = { 1, 57 };

  const mxArray *m1;
  char_T cv6[57];
  static const char_T cv7[57] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'e', 'm', 'l', '_', 'a', 's', 's', 'e', 'r', 't',
    '_', 'v', 'a', 'l', 'i', 'd', '_', 's', 'i', 'z', 'e', '_', 'a', 'r', 'g',
    '_', 'i', 'n', 'v', 'a', 'l', 'i', 'd', 'S', 'i', 'z', 'e', 'V', 'e', 'c',
    't', 'o', 'r' };

  const mxArray *b_y;
  const mxArray *c_y;
  real_T a;
  const mxArray *d_y;
  static const int32_T iv9[2] = { 1, 21 };

  char_T cv8[21];
  static const char_T cv9[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'p', 'm', 'a', 'x', 's', 'i', 'z', 'e' };

  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = sp;
  b_st.tls = sp->tls;
  st.site = &y_emlrtRSI;
  i = 0;
  guard1 = FALSE;
  do {
    exitg1 = 0;
    if (i < 2) {
      if (varargin_1[i] != varargin_1[i]) {
        guard1 = TRUE;
        exitg1 = 1;
      } else if (muDoubleScalarIsInf(varargin_1[i])) {
        guard1 = TRUE;
        exitg1 = 1;
      } else {
        i++;
        guard1 = FALSE;
      }
    } else {
      p = TRUE;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  if (guard1 == TRUE) {
    p = FALSE;
  }

  if (p) {
    st.site = &y_emlrtRSI;
    p = TRUE;
  } else {
    p = FALSE;
  }

  if (p) {
  } else {
    y = NULL;
    m1 = mxCreateCharArray(2, iv8);
    for (i = 0; i < 57; i++) {
      cv6[i] = cv7[i];
    }

    emlrtInitCharArrayR2013a(sp, 57, m1, cv6);
    emlrtAssign(&y, m1);
    b_y = NULL;
    m1 = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    *(int32_T *)mxGetData(m1) = MIN_int32_T;
    emlrtAssign(&b_y, m1);
    c_y = NULL;
    m1 = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    *(int32_T *)mxGetData(m1) = MAX_int32_T;
    emlrtAssign(&c_y, m1);
    st.site = &y_emlrtRSI;
    b_st.site = &pq_emlrtRSI;
    b_error(&st, message(&b_st, y, b_y, c_y, &emlrtMCI), &b_emlrtMCI);
  }

  st.site = &ab_emlrtRSI;
  a = 1.0;
  for (i = 0; i < 2; i++) {
    if (varargin_1[i] <= 0.0) {
      a = 0.0;
    } else {
      a *= varargin_1[i];
    }
  }

  st.site = &ab_emlrtRSI;
  if (2.147483647E+9 >= a) {
  } else {
    d_y = NULL;
    m1 = mxCreateCharArray(2, iv9);
    for (i = 0; i < 21; i++) {
      cv8[i] = cv9[i];
    }

    emlrtInitCharArrayR2013a(sp, 21, m1, cv8);
    emlrtAssign(&d_y, m1);
    st.site = &ab_emlrtRSI;
    b_st.site = &mq_emlrtRSI;
    b_error(&st, b_message(&b_st, d_y, &c_emlrtMCI), &d_emlrtMCI);
  }
}

void repmat(const emlrtStack *sp, const emxArray_real_T *a, const real_T m[2],
            emxArray_real_T *b)
{
  int32_T mv[2];
  int32_T ia;
  int32_T outsize[2];
  int32_T ib;
  int32_T iacol;
  boolean_T overflow;
  int32_T jcol;
  boolean_T b2;
  int32_T itilerow;
  emlrtStack st;
  emlrtStack b_st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &s_emlrtRSI;
  b_st.prev = &st;
  b_st.tls = st.tls;
  eml_assert_valid_size_arg(&st, m);
  for (ia = 0; ia < 2; ia++) {
    mv[ia] = (int32_T)m[ia];
  }

  for (ia = 0; ia < 2; ia++) {
    outsize[ia] = a->size[ia] * mv[ia];
  }

  ia = b->size[0] * b->size[1];
  b->size[0] = outsize[0];
  b->size[1] = outsize[1];
  emxEnsureCapacity(sp, (emxArray__common *)b, ia, (int32_T)sizeof(real_T),
                    &k_emlrtRTEI);
  if ((outsize[0] == 0) || (outsize[1] == 0)) {
  } else {
    ia = 1;
    ib = 0;
    st.site = &t_emlrtRSI;
    iacol = 1;
    st.site = &u_emlrtRSI;
    if (1 > a->size[1]) {
      overflow = FALSE;
    } else {
      overflow = (a->size[1] > 2147483646);
    }

    if (overflow) {
      b_st.site = &hb_emlrtRSI;
      check_forloop_overflow_error(&b_st);
    }

    for (jcol = 1; jcol <= a->size[1]; jcol++) {
      st.site = &v_emlrtRSI;
      if (1 > mv[0]) {
        b2 = FALSE;
      } else {
        b2 = (mv[0] > 2147483646);
      }

      if (b2) {
        b_st.site = &hb_emlrtRSI;
        check_forloop_overflow_error(&b_st);
      }

      for (itilerow = 1; itilerow <= mv[0]; itilerow++) {
        b->data[ib] = a->data[iacol - 1];
        st.site = &w_emlrtRSI;
        ia = iacol + 1;
        st.site = &x_emlrtRSI;
        ib++;
      }

      iacol = ia;
    }
  }
}

/* End of code generation (repmat.c) */
