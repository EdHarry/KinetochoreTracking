/*
 * bsxfun.c
 *
 * Code generation for function 'bsxfun'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "bsxfun.h"
#include "mmfMex_emxutil.h"
#include "eml_int_forloop_overflow_check.h"
#include "mmfMex_mexutil.h"
#include "mmfMex_data.h"

/* Variable Definitions */
static emlrtRSInfo db_emlrtRSI = { 24, "indexIntRelop",
  "/Applications/MATLAB_R2013b.app/toolbox/shared/coder/coder/+coder/+internal/indexIntRelop.m"
};

static emlrtRSInfo dp_emlrtRSI = { 21, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo ep_emlrtRSI = { 23, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo fp_emlrtRSI = { 75, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo gp_emlrtRSI = { 77, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo hp_emlrtRSI = { 78, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo ip_emlrtRSI = { 82, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo jp_emlrtRSI = { 87, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo kp_emlrtRSI = { 88, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo lp_emlrtRSI = { 96, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo mp_emlrtRSI = { 97, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo np_emlrtRSI = { 101, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo op_emlrtRSI = { 102, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo pp_emlrtRSI = { 104, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo qp_emlrtRSI = { 107, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo rp_emlrtRSI = { 108, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo sp_emlrtRSI = { 109, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo tp_emlrtRSI = { 110, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo vp_emlrtRSI = { 85, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtMCInfo bb_emlrtMCI = { 22, 5, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtMCInfo cb_emlrtMCI = { 21, 15, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtMCInfo db_emlrtMCI = { 24, 5, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtMCInfo eb_emlrtMCI = { 23, 15, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRTEInfo hb_emlrtRTEI = { 41, 1, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRTEInfo ib_emlrtRTEI = { 65, 1, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRTEInfo jb_emlrtRTEI = { 1, 14, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRTEInfo kb_emlrtRTEI = { 85, 5, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo eq_emlrtRSI = { 24, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

static emlrtRSInfo fq_emlrtRSI = { 22, "bsxfun",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/bsxfun.m" };

/* Function Definitions */
void b_bsxfun(const emlrtStack *sp, const emxArray_real_T *a, const real_T
              b_data[3], const int32_T b_size[2], emxArray_real_T *c)
{
  boolean_T overflow;
  const mxArray *y;
  static const int32_T iv26[2] = { 1, 44 };

  const mxArray *m13;
  char_T cv39[44];
  int32_T i;
  static const char_T cv40[44] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'b', 's', 'x', 'f', 'u', 'n', '_', 'a', 'r', 'r', 'a',
    'y', 'D', 'i', 'm', 'e', 'n', 's', 'i', 'o', 'n', 's', 'M', 'u', 's', 't',
    'M', 'a', 't', 'c', 'h' };

  const mxArray *b_y;
  static const int32_T iv27[2] = { 1, 37 };

  char_T cv41[37];
  static const char_T cv42[37] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'b', 's', 'x', 'f', 'u', 'n', '_', 'd', 'y', 'n',
    'a', 'm', 'i', 'c', 'E', 'x', 'p', 'a', 'n', 's', 'i', 'o', 'n' };

  int32_T a_idx_0;
  int32_T b_a_idx_0;
  emxArray_real_T *av;
  int32_T asub;
  int32_T bsub;
  int32_T ak;
  int32_T bk;
  int32_T b;
  int32_T ck;
  emxArray_real_T *cv;
  int32_T exitg1;
  int32_T exitg2;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = sp;
  b_st.tls = sp->tls;
  c_st.prev = &st;
  c_st.tls = st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  st.site = &dp_emlrtRSI;
  if ((a->size[1] != 1) && (b_size[1] != 1) && (a->size[1] != b_size[1])) {
    overflow = FALSE;
  } else {
    overflow = TRUE;
  }

  if (overflow) {
  } else {
    y = NULL;
    m13 = mxCreateCharArray(2, iv26);
    for (i = 0; i < 44; i++) {
      cv39[i] = cv40[i];
    }

    emlrtInitCharArrayR2013a(sp, 44, m13, cv39);
    emlrtAssign(&y, m13);
    st.site = &dp_emlrtRSI;
    b_st.site = &fq_emlrtRSI;
    b_error(&st, b_message(&b_st, y, &bb_emlrtMCI), &cb_emlrtMCI);
  }

  st.site = &ep_emlrtRSI;
  if (a->size[1] != b_size[1]) {
    overflow = FALSE;
  } else {
    overflow = TRUE;
  }

  if (overflow) {
  } else {
    b_y = NULL;
    m13 = mxCreateCharArray(2, iv27);
    for (i = 0; i < 37; i++) {
      cv41[i] = cv42[i];
    }

    emlrtInitCharArrayR2013a(sp, 37, m13, cv41);
    emlrtAssign(&b_y, m13);
    st.site = &ep_emlrtRSI;
    b_st.site = &eq_emlrtRSI;
    b_error(&st, b_message(&b_st, b_y, &db_emlrtMCI), &eb_emlrtMCI);
  }

  a_idx_0 = c->size[0] * c->size[1];
  c->size[0] = a->size[0];
  c->size[1] = a->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)c, a_idx_0, (int32_T)sizeof(real_T),
                    &hb_emlrtRTEI);
  b_a_idx_0 = a->size[0];
  i = a->size[1];
  if ((b_a_idx_0 == 0) || (i == 0)) {
  } else {
    c_emxInit_real_T(sp, &av, 1, &ib_emlrtRTEI, TRUE);
    i = a->size[0];
    a_idx_0 = av->size[0];
    av->size[0] = i;
    emxEnsureCapacity(sp, (emxArray__common *)av, a_idx_0, (int32_T)sizeof
                      (real_T), &ib_emlrtRTEI);
    asub = 1;
    bsub = 1;
    ak = -1;
    bk = 0;
    st.site = &fp_emlrtRSI;
    b_a_idx_0 = a->size[0];
    i = a->size[1];
    a_idx_0 = a->size[0];
    b = b_a_idx_0 * i - a_idx_0;
    st.site = &fp_emlrtRSI;
    c_st.site = &gb_emlrtRSI;
    b_a_idx_0 = a->size[0];
    if ((b_a_idx_0 == 0) || (0 > b)) {
      overflow = FALSE;
    } else {
      b_a_idx_0 = a->size[0];
      overflow = (b > MAX_int32_T - b_a_idx_0);
    }

    if (overflow) {
      c_st.site = &hb_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }

    ck = 0;
    c_emxInit_real_T(sp, &cv, 1, &kb_emlrtRTEI, TRUE);
    do {
      exitg1 = 0;
      b_a_idx_0 = a->size[0];
      if ((b_a_idx_0 > 0) && (ck <= b)) {
        st.site = &gp_emlrtRSI;
        c_st.site = &gb_emlrtRSI;
        if (1 > a->size[0]) {
          overflow = FALSE;
        } else {
          overflow = (a->size[0] > 2147483646);
        }

        if (overflow) {
          c_st.site = &hb_emlrtRSI;
          check_forloop_overflow_error(&c_st);
        }

        for (i = 1; i <= a->size[0]; i++) {
          st.site = &hp_emlrtRSI;
          av->data[i - 1] = a->data[ak + i];
        }

        st.site = &ip_emlrtRSI;
        st.site = &vp_emlrtRSI;
        c_st.site = &ib_emlrtRSI;
        a_idx_0 = cv->size[0];
        cv->size[0] = av->size[0];
        emxEnsureCapacity(&c_st, (emxArray__common *)cv, a_idx_0, (int32_T)
                          sizeof(real_T), &jb_emlrtRTEI);
        i = av->size[0];
        for (a_idx_0 = 0; a_idx_0 < i; a_idx_0++) {
          cv->data[a_idx_0] = av->data[a_idx_0] / b_data[bk];
        }

        st.site = &jp_emlrtRSI;
        c_st.site = &gb_emlrtRSI;
        b_a_idx_0 = a->size[0];
        if (1 > b_a_idx_0) {
          overflow = FALSE;
        } else {
          b_a_idx_0 = a->size[0];
          overflow = (b_a_idx_0 > 2147483646);
        }

        if (overflow) {
          c_st.site = &hb_emlrtRSI;
          check_forloop_overflow_error(&c_st);
        }

        i = 0;
        do {
          exitg2 = 0;
          b_a_idx_0 = a->size[0];
          if (i + 1 <= b_a_idx_0) {
            st.site = &kp_emlrtRSI;
            c->data[ck + i] = cv->data[i];
            i++;
          } else {
            exitg2 = 1;
          }
        } while (exitg2 == 0);

        st.site = &lp_emlrtRSI;
        c_st.site = &db_emlrtRSI;
        if (asub < a->size[1]) {
          st.site = &mp_emlrtRSI;
          ak += a->size[0];
          st.site = &np_emlrtRSI;
          bk++;
          st.site = &op_emlrtRSI;
          bsub++;
          st.site = &pp_emlrtRSI;
          asub++;
        } else {
          st.site = &qp_emlrtRSI;
          c_st.site = &db_emlrtRSI;
          if (bsub < b_size[1]) {
            st.site = &rp_emlrtRSI;
            st.site = &sp_emlrtRSI;
            bk++;
            st.site = &tp_emlrtRSI;
            bsub++;
          } else {
            asub = 1;
            bsub = 1;
          }
        }

        b_a_idx_0 = a->size[0];
        ck += b_a_idx_0;
      } else {
        exitg1 = 1;
      }
    } while (exitg1 == 0);

    emxFree_real_T(&cv);
    emxFree_real_T(&av);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

void bsxfun(const emlrtStack *sp, const emxArray_real_T *a, const real_T b_data
            [3], const int32_T b_size[2], emxArray_real_T *c)
{
  boolean_T overflow;
  const mxArray *y;
  static const int32_T iv24[2] = { 1, 44 };

  const mxArray *m12;
  char_T cv35[44];
  int32_T i;
  static const char_T cv36[44] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'b', 's', 'x', 'f', 'u', 'n', '_', 'a', 'r', 'r', 'a',
    'y', 'D', 'i', 'm', 'e', 'n', 's', 'i', 'o', 'n', 's', 'M', 'u', 's', 't',
    'M', 'a', 't', 'c', 'h' };

  const mxArray *b_y;
  static const int32_T iv25[2] = { 1, 37 };

  char_T cv37[37];
  static const char_T cv38[37] = { 'C', 'o', 'd', 'e', 'r', ':', 't', 'o', 'o',
    'l', 'b', 'o', 'x', ':', 'b', 's', 'x', 'f', 'u', 'n', '_', 'd', 'y', 'n',
    'a', 'm', 'i', 'c', 'E', 'x', 'p', 'a', 'n', 's', 'i', 'o', 'n' };

  int32_T a_idx_0;
  int32_T b_a_idx_0;
  emxArray_real_T *av;
  int32_T asub;
  int32_T bsub;
  int32_T ak;
  int32_T bk;
  int32_T b;
  int32_T ck;
  emxArray_real_T *cv;
  int32_T exitg1;
  int32_T exitg2;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = sp;
  b_st.tls = sp->tls;
  c_st.prev = &st;
  c_st.tls = st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  st.site = &dp_emlrtRSI;
  if ((a->size[1] != 1) && (b_size[1] != 1) && (a->size[1] != b_size[1])) {
    overflow = FALSE;
  } else {
    overflow = TRUE;
  }

  if (overflow) {
  } else {
    y = NULL;
    m12 = mxCreateCharArray(2, iv24);
    for (i = 0; i < 44; i++) {
      cv35[i] = cv36[i];
    }

    emlrtInitCharArrayR2013a(sp, 44, m12, cv35);
    emlrtAssign(&y, m12);
    st.site = &dp_emlrtRSI;
    b_st.site = &fq_emlrtRSI;
    b_error(&st, b_message(&b_st, y, &bb_emlrtMCI), &cb_emlrtMCI);
  }

  st.site = &ep_emlrtRSI;
  if (a->size[1] != b_size[1]) {
    overflow = FALSE;
  } else {
    overflow = TRUE;
  }

  if (overflow) {
  } else {
    b_y = NULL;
    m12 = mxCreateCharArray(2, iv25);
    for (i = 0; i < 37; i++) {
      cv37[i] = cv38[i];
    }

    emlrtInitCharArrayR2013a(sp, 37, m12, cv37);
    emlrtAssign(&b_y, m12);
    st.site = &ep_emlrtRSI;
    b_st.site = &eq_emlrtRSI;
    b_error(&st, b_message(&b_st, b_y, &db_emlrtMCI), &eb_emlrtMCI);
  }

  a_idx_0 = c->size[0] * c->size[1];
  c->size[0] = a->size[0];
  c->size[1] = a->size[1];
  emxEnsureCapacity(sp, (emxArray__common *)c, a_idx_0, (int32_T)sizeof(real_T),
                    &hb_emlrtRTEI);
  b_a_idx_0 = a->size[0];
  i = a->size[1];
  if ((b_a_idx_0 == 0) || (i == 0)) {
  } else {
    c_emxInit_real_T(sp, &av, 1, &ib_emlrtRTEI, TRUE);
    i = a->size[0];
    a_idx_0 = av->size[0];
    av->size[0] = i;
    emxEnsureCapacity(sp, (emxArray__common *)av, a_idx_0, (int32_T)sizeof
                      (real_T), &ib_emlrtRTEI);
    asub = 1;
    bsub = 1;
    ak = -1;
    bk = 0;
    st.site = &fp_emlrtRSI;
    b_a_idx_0 = a->size[0];
    i = a->size[1];
    a_idx_0 = a->size[0];
    b = b_a_idx_0 * i - a_idx_0;
    st.site = &fp_emlrtRSI;
    b_a_idx_0 = a->size[0];
    if ((b_a_idx_0 == 0) || (0 > b)) {
      overflow = FALSE;
    } else {
      b_a_idx_0 = a->size[0];
      overflow = (b > MAX_int32_T - b_a_idx_0);
    }

    if (overflow) {
      c_st.site = &hb_emlrtRSI;
      check_forloop_overflow_error(&c_st);
    }

    ck = 0;
    c_emxInit_real_T(sp, &cv, 1, &kb_emlrtRTEI, TRUE);
    do {
      exitg1 = 0;
      b_a_idx_0 = a->size[0];
      if ((b_a_idx_0 > 0) && (ck <= b)) {
        st.site = &gp_emlrtRSI;
        if (1 > a->size[0]) {
          overflow = FALSE;
        } else {
          overflow = (a->size[0] > 2147483646);
        }

        if (overflow) {
          c_st.site = &hb_emlrtRSI;
          check_forloop_overflow_error(&c_st);
        }

        for (i = 1; i <= a->size[0]; i++) {
          st.site = &hp_emlrtRSI;
          av->data[i - 1] = a->data[ak + i];
        }

        st.site = &ip_emlrtRSI;
        a_idx_0 = cv->size[0];
        cv->size[0] = av->size[0];
        emxEnsureCapacity(sp, (emxArray__common *)cv, a_idx_0, (int32_T)sizeof
                          (real_T), &jb_emlrtRTEI);
        i = av->size[0];
        for (a_idx_0 = 0; a_idx_0 < i; a_idx_0++) {
          cv->data[a_idx_0] = av->data[a_idx_0] - b_data[bk];
        }

        st.site = &jp_emlrtRSI;
        b_a_idx_0 = a->size[0];
        if (1 > b_a_idx_0) {
          overflow = FALSE;
        } else {
          b_a_idx_0 = a->size[0];
          overflow = (b_a_idx_0 > 2147483646);
        }

        if (overflow) {
          c_st.site = &hb_emlrtRSI;
          check_forloop_overflow_error(&c_st);
        }

        i = 1;
        do {
          exitg2 = 0;
          b_a_idx_0 = a->size[0];
          if (i <= b_a_idx_0) {
            st.site = &kp_emlrtRSI;
            c->data[(ck + i) - 1] = cv->data[i - 1];
            i++;
          } else {
            exitg2 = 1;
          }
        } while (exitg2 == 0);

        st.site = &lp_emlrtRSI;
        if (asub < a->size[1]) {
          st.site = &mp_emlrtRSI;
          ak += a->size[0];
          st.site = &np_emlrtRSI;
          bk++;
          st.site = &op_emlrtRSI;
          bsub++;
          st.site = &pp_emlrtRSI;
          asub++;
        } else {
          st.site = &qp_emlrtRSI;
          if (bsub < b_size[1]) {
            st.site = &rp_emlrtRSI;
            st.site = &sp_emlrtRSI;
            bk++;
            st.site = &tp_emlrtRSI;
            bsub++;
          } else {
            asub = 1;
            bsub = 1;
          }
        }

        b_a_idx_0 = a->size[0];
        ck += b_a_idx_0;
      } else {
        exitg1 = 1;
      }
    } while (exitg1 == 0);

    emxFree_real_T(&cv);
    emxFree_real_T(&av);
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (bsxfun.c) */
