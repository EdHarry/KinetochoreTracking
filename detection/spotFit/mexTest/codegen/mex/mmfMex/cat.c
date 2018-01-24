/*
 * cat.c
 *
 * Code generation for function 'cat'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "cat.h"
#include "mmfMex_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "eml_error.h"
#include "mmfMex_mexutil.h"
#include "mmfMex_data.h"

/* Function Definitions */
void cat(const emlrtStack *sp, const emxArray_real_T *varargin_1, const
         emxArray_real_T *varargin_2, emxArray_real_T *y)
{
  uint32_T ysize[3];
  int32_T i;
  uint32_T sz1[2];
  int32_T exitg2;
  int32_T b_ysize[3];
  boolean_T p;
  const mxArray *b_y;
  static const int32_T iv15[2] = { 1, 39 };

  const mxArray *m5;
  char_T cv20[39];
  static const char_T cv21[39] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'c', 'a', 't', 'e', 'n', 'a', 't', 'e', '_', 'd', 'i',
    'm', 'e', 'n', 's', 'i', 'o', 'n', 'M', 'i', 's', 'm', 'a', 't', 'c', 'h' };

  int32_T exitg1;
  int32_T c_ysize[3];
  const mxArray *c_y;
  static const int32_T iv16[2] = { 1, 39 };

  int32_T iy;
  int32_T b;
  boolean_T b3;
  boolean_T b4;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = sp;
  b_st.tls = sp->tls;
  c_st.prev = &st;
  c_st.tls = st.tls;
  st.site = &xb_emlrtRSI;
  for (i = 0; i < 3; i++) {
    ysize[i] = 1U;
  }

  for (i = 0; i < 2; i++) {
    sz1[i] = (uint32_T)varargin_1->size[i];
  }

  for (i = 0; i < 2; i++) {
    ysize[i] = sz1[i];
  }

  i = y->size[0] * y->size[1] * y->size[2];
  y->size[0] = (int32_T)ysize[0];
  y->size[1] = (int32_T)ysize[1];
  y->size[2] = 2;
  emxEnsureCapacity(sp, (emxArray__common *)y, i, (int32_T)sizeof(real_T),
                    &c_emlrtRTEI);
  if ((varargin_1->size[0] == 0) || (varargin_1->size[1] == 0)) {
    st.site = &yb_emlrtRSI;
    eml_error(&st);
  }

  st.site = &ac_emlrtRSI;
  i = 0;
  do {
    exitg2 = 0;
    if (i < 2) {
      b_ysize[0] = (int32_T)ysize[0];
      b_ysize[1] = (int32_T)ysize[1];
      b_ysize[2] = 2;
      if (b_ysize[i] != varargin_1->size[i]) {
        p = FALSE;
        exitg2 = 1;
      } else {
        i++;
      }
    } else {
      p = TRUE;
      exitg2 = 1;
    }
  } while (exitg2 == 0);

  if (p) {
  } else {
    b_y = NULL;
    m5 = mxCreateCharArray(2, iv15);
    for (i = 0; i < 39; i++) {
      cv20[i] = cv21[i];
    }

    emlrtInitCharArrayR2013a(sp, 39, m5, cv20);
    emlrtAssign(&b_y, m5);
    st.site = &ac_emlrtRSI;
    b_st.site = &iq_emlrtRSI;
    b_error(&st, b_message(&b_st, b_y, &m_emlrtMCI), &n_emlrtMCI);
  }

  if ((varargin_2->size[0] == 0) || (varargin_2->size[1] == 0)) {
    st.site = &yb_emlrtRSI;
    eml_error(&st);
  }

  st.site = &ac_emlrtRSI;
  i = 0;
  do {
    exitg1 = 0;
    if (i < 2) {
      c_ysize[0] = (int32_T)ysize[0];
      c_ysize[1] = (int32_T)ysize[1];
      c_ysize[2] = 2;
      if (c_ysize[i] != varargin_2->size[i]) {
        p = FALSE;
        exitg1 = 1;
      } else {
        i++;
      }
    } else {
      p = TRUE;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  if (p) {
  } else {
    c_y = NULL;
    m5 = mxCreateCharArray(2, iv16);
    for (i = 0; i < 39; i++) {
      cv20[i] = cv21[i];
    }

    emlrtInitCharArrayR2013a(sp, 39, m5, cv20);
    emlrtAssign(&c_y, m5);
    st.site = &ac_emlrtRSI;
    b_st.site = &iq_emlrtRSI;
    b_error(&st, b_message(&b_st, c_y, &m_emlrtMCI), &n_emlrtMCI);
  }

  iy = -1;
  b = varargin_1->size[0] * varargin_1->size[1];
  st.site = &bc_emlrtRSI;
  if (1 > b) {
    b3 = FALSE;
  } else {
    b3 = (b > 2147483646);
  }

  if (b3) {
    c_st.site = &hb_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  for (i = 1; i <= b; i++) {
    st.site = &cc_emlrtRSI;
    iy++;
    y->data[iy] = varargin_1->data[i - 1];
  }

  b = varargin_2->size[0] * varargin_2->size[1];
  st.site = &bc_emlrtRSI;
  if (1 > b) {
    b4 = FALSE;
  } else {
    b4 = (b > 2147483646);
  }

  if (b4) {
    c_st.site = &hb_emlrtRSI;
    check_forloop_overflow_error(&c_st);
  }

  for (i = 1; i <= b; i++) {
    st.site = &cc_emlrtRSI;
    iy++;
    y->data[iy] = varargin_2->data[i - 1];
  }
}

boolean_T isconsistent(const emxArray_real_T *y, const emxArray_real_T *x)
{
  boolean_T p;
  int32_T j;
  int32_T exitg1;
  j = 0;
  do {
    exitg1 = 0;
    if (j < 2) {
      if ((j + 1 != 2) && (y->size[0] != x->size[0])) {
        p = FALSE;
        exitg1 = 1;
      } else {
        j++;
      }
    } else {
      p = TRUE;
      exitg1 = 1;
    }
  } while (exitg1 == 0);

  return p;
}

/* End of code generation (cat.c) */
