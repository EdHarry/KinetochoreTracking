/*
 * repeatEntries_mex.c
 *
 * Code generation for function 'repeatEntries_mex'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "repeatEntries_mex.h"
#include "round.h"
#include "mmfMex_emxutil.h"
#include "repmat.h"
#include "error.h"
#include "mmfMex_mexutil.h"
#include "mmfMex_data.h"

/* Variable Definitions */
static emlrtRSInfo fb_emlrtRSI = { 20, "eml_index_prod",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/eml/eml_index_prod.m"
};

static emlrtRSInfo wn_emlrtRSI = { 25, "repeatEntries_mex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/repeatEntries_mex.m" };

static emlrtRSInfo xn_emlrtRSI = { 54, "repeatEntries_mex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/repeatEntries_mex.m" };

static emlrtRSInfo yn_emlrtRSI = { 82, "repeatEntries_mex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/repeatEntries_mex.m" };

static emlrtRSInfo bo_emlrtRSI = { 79, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo fo_emlrtRSI = { 282, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo go_emlrtRSI = { 283, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo ho_emlrtRSI = { 289, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo io_emlrtRSI = { 290, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo jo_emlrtRSI = { 291, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo ko_emlrtRSI = { 292, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo lo_emlrtRSI = { 296, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo mo_emlrtRSI = { 299, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo no_emlrtRSI = { 293, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo oo_emlrtRSI = { 294, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo po_emlrtRSI = { 297, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo qo_emlrtRSI = { 300, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo ro_emlrtRSI = { 301, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo so_emlrtRSI = { 348, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo to_emlrtRSI = { 344, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo uo_emlrtRSI = { 337, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo vo_emlrtRSI = { 336, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo wo_emlrtRSI = { 324, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo xo_emlrtRSI = { 314, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo yo_emlrtRSI = { 17, "meshgrid",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/meshgrid.m" };

static emlrtRSInfo ap_emlrtRSI = { 18, "meshgrid",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/meshgrid.m" };

static emlrtRSInfo bp_emlrtRSI = { 54, "repmat",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/repmat.m" };

static emlrtRSInfo cp_emlrtRSI = { 63, "repmat",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/elmat/repmat.m" };

static emlrtMCInfo p_emlrtMCI = { 405, 5, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtMCInfo q_emlrtMCI = { 404, 15, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRTEInfo eb_emlrtRTEI = { 1, 16, "repeatEntries_mex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/repeatEntries_mex.m" };

static emlrtRTEInfo fb_emlrtRTEI = { 284, 1, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRTEInfo gb_emlrtRTEI = { 82, 5, "repeatEntries_mex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/repeatEntries_mex.m" };

static emlrtBCInfo q_emlrtBCI = { -1, -1, 85, 15, "val", "repeatEntries_mex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/repeatEntries_mex.m", 0 };

static emlrtDCInfo emlrtDCI = { 85, 15, "repeatEntries_mex",
  "/Users/edwardharry/Documents/MATLAB/spotFit/mexTest/repeatEntries_mex.m", 1 };

static emlrtRSInfo xp_emlrtRSI = { 404, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

static emlrtRSInfo gq_emlrtRSI = { 405, "colon",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/ops/colon.m" };

/* Function Definitions */
void repeatEntries_mex(const emlrtStack *sp, const emxArray_real_T *val, real_T
  kTimes, emxArray_real_T *out)
{
  real_T valSize[2];
  int32_T i5;
  int32_T cdiff;
  int32_T ndbl;
  int32_T apnd;
  emxArray_real_T *x;
  int32_T k;
  real_T anew;
  real_T b_apnd;
  boolean_T n_too_large;
  real_T b_ndbl;
  real_T b_cdiff;
  const mxArray *y;
  static const int32_T iv22[2] = { 1, 21 };

  const mxArray *m10;
  char_T cv31[21];
  static const char_T cv32[21] = { 'C', 'o', 'd', 'e', 'r', ':', 'M', 'A', 'T',
    'L', 'A', 'B', ':', 'p', 'm', 'a', 'x', 's', 'i', 'z', 'e' };

  emxArray_real_T *b_y;
  emxArray_real_T *idxMat;
  emxArray_real_T *r5;
  emxArray_real_T *b_x;
  real_T c_y[2];
  int32_T d_y[2];
  int32_T outsize[2];
  int32_T mv[2];
  int32_T exitg1;
  int32_T i6;
  emlrtStack st;
  emlrtStack b_st;
  emlrtStack c_st;
  emlrtStack d_st;
  emlrtStack e_st;
  emlrtStack f_st;
  st.prev = sp;
  st.tls = sp->tls;
  b_st.prev = &st;
  b_st.tls = st.tls;
  c_st.prev = &b_st;
  c_st.tls = b_st.tls;
  d_st.prev = &c_st;
  d_st.tls = c_st.tls;
  e_st.prev = &d_st;
  e_st.tls = d_st.tls;
  f_st.prev = &d_st;
  f_st.tls = d_st.tls;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);

  /* REPEATENTRIES fills a matrix with k repeats the rows of the input matrix */
  /*  */
  /*  SYNOPSIS out = repeatEntries(val,kTimes) */
  /*  */
  /*  INPUT    val    : matrix (or vectors) containing the rows to repeat (works for strings, too) */
  /*           kTimes : number of repeats of each row (scalar or vector of size(vlaues,1)) */
  /*  */
  /*  OUTPUT   out    : matrix of size [sum(kTimes) size(values,2)] containing */
  /*                    repeated entries specified with k */
  /*  */
  /*  EXAMPLES     repeatEntries([1;2;3;4],[2;3;1;1]) returns [1;1;2;2;2;3;4] */
  /*  */
  /*               repeatEntries([1;2;3;4],2) returns [1;1;2;2;3;3;4;4] */
  /*  */
  /*  c: jonas, 2/04 */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* =========== */
  /*  test input */
  /* =========== */
  /*  nargin */
  if ((val->size[0] == 0) || (val->size[1] == 0)) {
    st.site = &wn_emlrtRSI;
    error(&st);
  }

  /*  size */
  for (i5 = 0; i5 < 2; i5++) {
    valSize[i5] = val->size[i5];
  }

  /*  decide whether we have scalar k */
  /*  do not care about size of k: we want to make a col vector out of it - and */
  /*  this vector should only contain nonzero positive integers */
  st.site = &xn_emlrtRSI;
  b_round(&kTimes);

  /*  if there are any negative values or zeros, remove the entry */
  if (kTimes < 1.0) {
    i5 = out->size[0] * out->size[1];
    out->size[0] = 0;
    out->size[1] = 0;
    emxEnsureCapacity(sp, (emxArray__common *)out, i5, (int32_T)sizeof(real_T),
                      &eb_emlrtRTEI);
  } else {
    /* kTimes = max(kTimes,ones(size(kTimes))); */
    /* ============ */
    /*  fill in out */
    /* ============ */
    /*  first the elegant case: scalar k */
    /*  build repeat index matrix idxMat */
    st.site = &yn_emlrtRSI;
    cdiff = val->size[0];
    b_st.site = &mf_emlrtRSI;
    c_st.site = &bo_emlrtRSI;
    d_st.site = &fo_emlrtRSI;
    e_st.site = &xo_emlrtRSI;
    e_st.site = &wo_emlrtRSI;
    e_st.site = &vo_emlrtRSI;
    ndbl = (int32_T)muDoubleScalarFloor(((real_T)cdiff - 1.0) + 0.5);
    e_st.site = &uo_emlrtRSI;
    apnd = ndbl + 1;
    cdiff = (ndbl - (int32_T)valSize[0]) + 1;
    e_st.site = &to_emlrtRSI;
    e_st.site = &to_emlrtRSI;
    e_st.site = &to_emlrtRSI;
    if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * valSize[0]) {
      ndbl++;
      apnd = (int32_T)valSize[0];
    } else if (cdiff > 0) {
      e_st.site = &so_emlrtRSI;
      apnd = ndbl;
    } else {
      ndbl++;
    }

    emxInit_real_T(&d_st, &x, 2, &eb_emlrtRTEI, TRUE);
    d_st.site = &go_emlrtRSI;
    i5 = x->size[0] * x->size[1];
    x->size[0] = 1;
    x->size[1] = ndbl;
    emxEnsureCapacity(&c_st, (emxArray__common *)x, i5, (int32_T)sizeof(real_T),
                      &fb_emlrtRTEI);
    if (ndbl > 0) {
      x->data[0] = 1.0;
      if (ndbl > 1) {
        x->data[ndbl - 1] = apnd;
        d_st.site = &ho_emlrtRSI;
        d_st.site = &io_emlrtRSI;
        i5 = ndbl - 1;
        i5 += (i5 < 0);
        if (i5 >= 0) {
          cdiff = (int32_T)((uint32_T)i5 >> 1);
        } else {
          cdiff = ~(int32_T)((uint32_T)~i5 >> 1);
        }

        d_st.site = &jo_emlrtRSI;
        d_st.site = &jo_emlrtRSI;
        e_st.site = &gb_emlrtRSI;
        for (k = 1; k < cdiff; k++) {
          d_st.site = &ko_emlrtRSI;
          d_st.site = &no_emlrtRSI;
          x->data[k] = 1.0 + (real_T)k;
          d_st.site = &oo_emlrtRSI;
          x->data[(ndbl - k) - 1] = apnd - k;
        }

        d_st.site = &lo_emlrtRSI;
        if (cdiff << 1 == ndbl - 1) {
          d_st.site = &po_emlrtRSI;
          x->data[cdiff] = (1.0 + (real_T)apnd) / 2.0;
        } else {
          d_st.site = &mo_emlrtRSI;
          d_st.site = &qo_emlrtRSI;
          x->data[cdiff] = 1.0 + (real_T)cdiff;
          d_st.site = &ro_emlrtRSI;
          x->data[cdiff + 1] = apnd - cdiff;
        }
      }
    }

    st.site = &yn_emlrtRSI;
    b_st.site = &mf_emlrtRSI;
    c_st.site = &bo_emlrtRSI;
    d_st.site = &fo_emlrtRSI;
    e_st.site = &xo_emlrtRSI;
    if (muDoubleScalarIsNaN(kTimes)) {
      ndbl = 0;
      anew = rtNaN;
      b_apnd = kTimes;
      n_too_large = FALSE;
    } else {
      e_st.site = &wo_emlrtRSI;
      if (muDoubleScalarIsInf(kTimes)) {
        ndbl = 0;
        anew = rtNaN;
        b_apnd = kTimes;
        n_too_large = !(1.0 == kTimes);
      } else {
        anew = 1.0;
        e_st.site = &vo_emlrtRSI;
        b_ndbl = muDoubleScalarFloor((kTimes - 1.0) + 0.5);
        e_st.site = &uo_emlrtRSI;
        b_apnd = 1.0 + b_ndbl;
        b_cdiff = (1.0 + b_ndbl) - kTimes;
        e_st.site = &to_emlrtRSI;
        e_st.site = &to_emlrtRSI;
        e_st.site = &to_emlrtRSI;
        if (muDoubleScalarAbs(b_cdiff) < 4.4408920985006262E-16 * kTimes) {
          b_ndbl++;
          b_apnd = kTimes;
        } else if (b_cdiff > 0.0) {
          e_st.site = &so_emlrtRSI;
          b_apnd = 1.0 + (b_ndbl - 1.0);
        } else {
          b_ndbl++;
        }

        n_too_large = FALSE;
        ndbl = (int32_T)b_ndbl - 1;
      }
    }

    d_st.site = &go_emlrtRSI;
    if (!n_too_large) {
    } else {
      y = NULL;
      m10 = mxCreateCharArray(2, iv22);
      for (cdiff = 0; cdiff < 21; cdiff++) {
        cv31[cdiff] = cv32[cdiff];
      }

      emlrtInitCharArrayR2013a(&d_st, 21, m10, cv31);
      emlrtAssign(&y, m10);
      e_st.site = &xp_emlrtRSI;
      f_st.site = &gq_emlrtRSI;
      b_error(&e_st, b_message(&f_st, y, &p_emlrtMCI), &q_emlrtMCI);
    }

    emxInit_real_T(&d_st, &b_y, 2, &eb_emlrtRTEI, TRUE);
    i5 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    b_y->size[1] = ndbl + 1;
    emxEnsureCapacity(&c_st, (emxArray__common *)b_y, i5, (int32_T)sizeof(real_T),
                      &fb_emlrtRTEI);
    if (ndbl + 1 > 0) {
      b_y->data[0] = anew;
      if (ndbl + 1 > 1) {
        b_y->data[ndbl] = b_apnd;
        d_st.site = &ho_emlrtRSI;
        d_st.site = &io_emlrtRSI;
        i5 = ndbl + (ndbl < 0);
        if (i5 >= 0) {
          cdiff = (int32_T)((uint32_T)i5 >> 1);
        } else {
          cdiff = ~(int32_T)((uint32_T)~i5 >> 1);
        }

        d_st.site = &jo_emlrtRSI;
        d_st.site = &jo_emlrtRSI;
        e_st.site = &gb_emlrtRSI;
        for (k = 1; k < cdiff; k++) {
          d_st.site = &ko_emlrtRSI;
          d_st.site = &no_emlrtRSI;
          b_y->data[k] = anew + (real_T)k;
          d_st.site = &oo_emlrtRSI;
          b_y->data[ndbl - k] = b_apnd - (real_T)k;
        }

        d_st.site = &lo_emlrtRSI;
        if (cdiff << 1 == ndbl) {
          d_st.site = &po_emlrtRSI;
          b_y->data[cdiff] = (anew + b_apnd) / 2.0;
        } else {
          d_st.site = &mo_emlrtRSI;
          d_st.site = &qo_emlrtRSI;
          b_y->data[cdiff] = anew + (real_T)cdiff;
          d_st.site = &ro_emlrtRSI;
          b_y->data[cdiff + 1] = b_apnd - (real_T)cdiff;
        }
      }
    }

    st.site = &yn_emlrtRSI;
    emxInit_real_T(&st, &idxMat, 2, &gb_emlrtRTEI, TRUE);
    emxInit_real_T(&st, &r5, 2, &eb_emlrtRTEI, TRUE);
    emxInit_real_T(&st, &b_x, 2, &eb_emlrtRTEI, TRUE);
    if ((x->size[1] == 0) || (b_y->size[1] == 0)) {
      i5 = idxMat->size[0] * idxMat->size[1];
      idxMat->size[0] = 0;
      idxMat->size[1] = 0;
      emxEnsureCapacity(&st, (emxArray__common *)idxMat, i5, (int32_T)sizeof
                        (real_T), &eb_emlrtRTEI);
    } else {
      cdiff = x->size[1];
      i5 = b_x->size[0] * b_x->size[1];
      b_x->size[0] = 1;
      b_x->size[1] = cdiff;
      emxEnsureCapacity(&st, (emxArray__common *)b_x, i5, (int32_T)sizeof(real_T),
                        &eb_emlrtRTEI);
      for (i5 = 0; i5 < cdiff; i5++) {
        b_x->data[b_x->size[0] * i5] = x->data[i5];
      }

      c_y[0] = b_y->size[1];
      c_y[1] = 1.0;
      b_st.site = &yo_emlrtRSI;
      repmat(&b_st, b_x, c_y, r5);
      i5 = idxMat->size[0] * idxMat->size[1];
      idxMat->size[0] = r5->size[0];
      idxMat->size[1] = r5->size[1];
      emxEnsureCapacity(&st, (emxArray__common *)idxMat, i5, (int32_T)sizeof
                        (real_T), &eb_emlrtRTEI);
      cdiff = r5->size[0] * r5->size[1];
      for (i5 = 0; i5 < cdiff; i5++) {
        idxMat->data[i5] = r5->data[i5];
      }

      b_st.site = &ap_emlrtRSI;
      valSize[0] = 1.0;
      valSize[1] = x->size[1];
      c_st.site = &s_emlrtRSI;
      eml_assert_valid_size_arg(&c_st, valSize);
      cdiff = b_y->size[1];
      d_y[0] = cdiff;
      d_y[1] = 1;
      for (i5 = 0; i5 < 2; i5++) {
        outsize[i5] = d_y[i5] * (int32_T)valSize[i5];
        mv[i5] = (int32_T)valSize[i5];
      }

      if (outsize[0] == 0) {
      } else {
        c_st.site = &bp_emlrtRSI;
        d_st.site = &fb_emlrtRSI;
        c_st.site = &t_emlrtRSI;
        d_st.site = &gb_emlrtRSI;
        for (ndbl = 1; ndbl <= mv[1]; ndbl++) {
          c_st.site = &v_emlrtRSI;
          d_st.site = &gb_emlrtRSI;
          c_st.site = &cp_emlrtRSI;
          d_st.site = &gb_emlrtRSI;
          k = 1;
          do {
            exitg1 = 0;
            cdiff = b_y->size[1];
            if (k <= cdiff) {
              c_st.site = &w_emlrtRSI;
              c_st.site = &x_emlrtRSI;
              k++;
            } else {
              exitg1 = 1;
            }
          } while (exitg1 == 0);
        }
      }
    }

    emxFree_real_T(&b_x);
    emxFree_real_T(&r5);
    emxFree_real_T(&b_y);
    emxFree_real_T(&x);

    /*  returns [1;1...2;2;... etc] */
    ndbl = idxMat->size[0] * idxMat->size[1];
    cdiff = val->size[1];
    i5 = out->size[0] * out->size[1];
    out->size[0] = ndbl;
    out->size[1] = cdiff;
    emxEnsureCapacity(sp, (emxArray__common *)out, i5, (int32_T)sizeof(real_T),
                      &eb_emlrtRTEI);
    for (i5 = 0; i5 < cdiff; i5++) {
      for (apnd = 0; apnd < ndbl; apnd++) {
        k = val->size[0];
        anew = idxMat->data[apnd];
        i6 = (int32_T)emlrtIntegerCheckFastR2012b(anew, &emlrtDCI, sp);
        out->data[apnd + out->size[0] * i5] = val->data
          [(emlrtDynamicBoundsCheckFastR2012b(i6, 1, k, &q_emlrtBCI, sp) +
            val->size[0] * i5) - 1];
      }
    }

    emxFree_real_T(&idxMat);

    /*  second: the loop */
  }

  /*  if doScalar */
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (repeatEntries_mex.c) */
