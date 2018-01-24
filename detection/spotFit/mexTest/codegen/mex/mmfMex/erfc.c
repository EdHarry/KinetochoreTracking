/*
 * erfc.c
 *
 * Code generation for function 'erfc'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "erfc.h"
#include "mmfMex_emxutil.h"

/* Variable Definitions */
static emlrtRSInfo ec_emlrtRSI = { 9, "erfc",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/specfun/erfc.m" };

static emlrtRTEInfo m_emlrtRTEI = { 18, 1, "eml_erfcore",
  "/Applications/MATLAB_R2013b.app/toolbox/eml/lib/matlab/specfun/eml_erfcore.m"
};

/* Function Declarations */
static real_T scalar_erf(real_T x);

/* Function Definitions */
static real_T scalar_erf(real_T x)
{
  real_T y;
  real_T absx;
  real_T s;
  real_T R;
  real_T S;
  int32_T eint;

  /* ========================== COPYRIGHT NOTICE ============================ */
  /*  The algorithms for calculating ERF(X) and ERFC(X) are derived           */
  /*  from FDLIBM, which has the following notice:                            */
  /*                                                                          */
  /*  Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.       */
  /*                                                                          */
  /*  Developed at SunSoft, a Sun Microsystems, Inc. business.                */
  /*  Permission to use, copy, modify, and distribute this                    */
  /*  software is freely granted, provided that this notice                   */
  /*  is preserved.                                                           */
  /* =============================    END    ================================ */
  absx = muDoubleScalarAbs(x);
  if (muDoubleScalarIsNaN(x)) {
    y = x;
  } else if (muDoubleScalarIsInf(x)) {
    if (x < 0.0) {
      y = 2.0;
    } else {
      y = 0.0;
    }
  } else if (absx < 0.84375) {
    if (absx < 1.3877787807814457E-17) {
      y = 1.0 - x;
    } else {
      s = x * x;
      y = (0.12837916709551256 + s * (-0.3250421072470015 + s *
            (-0.02848174957559851 + s * (-0.0057702702964894416 + s *
              -2.3763016656650163E-5)))) / (1.0 + s * (0.39791722395915535 + s *
        (0.0650222499887673 + s * (0.0050813062818757656 + s *
        (0.00013249473800432164 + s * -3.9602282787753681E-6)))));
      if (x < 0.25) {
        y = 1.0 - (x + x * y);
      } else {
        y = 0.5 - (x * y + (x - 0.5));
      }
    }
  } else if (absx < 1.25) {
    if (x >= 0.0) {
      y = 0.15493708848953247 - (-0.0023621185607526594 + (absx - 1.0) *
        (0.41485611868374833 + (absx - 1.0) * (-0.37220787603570132 + (absx -
        1.0) * (0.31834661990116175 + (absx - 1.0) * (-0.11089469428239668 +
        (absx - 1.0) * (0.035478304325618236 + (absx - 1.0) *
                        -0.0021663755948687908)))))) / (1.0 + (absx - 1.0) *
        (0.10642088040084423 + (absx - 1.0) * (0.540397917702171 + (absx - 1.0) *
        (0.071828654414196266 + (absx - 1.0) * (0.12617121980876164 + (absx -
        1.0) * (0.013637083912029051 + (absx - 1.0) * 0.011984499846799107))))));
    } else {
      y = 1.0 + (0.84506291151046753 + (-0.0023621185607526594 + (absx - 1.0) *
                  (0.41485611868374833 + (absx - 1.0) * (-0.37220787603570132 +
                    (absx - 1.0) * (0.31834661990116175 + (absx - 1.0) *
        (-0.11089469428239668 + (absx - 1.0) * (0.035478304325618236 + (absx -
        1.0) * -0.0021663755948687908)))))) / (1.0 + (absx - 1.0) *
                  (0.10642088040084423 + (absx - 1.0) * (0.540397917702171 +
        (absx - 1.0) * (0.071828654414196266 + (absx - 1.0) *
                        (0.12617121980876164 + (absx - 1.0) *
                         (0.013637083912029051 + (absx - 1.0) *
                          0.011984499846799107)))))));
    }
  } else if (x < -6.0) {
    y = 2.0;
  } else if (x >= 28.0) {
    y = 0.0;
  } else {
    s = 1.0 / (absx * absx);
    if (absx < 2.8571414947509766) {
      R = -0.0098649440348471482 + s * (-0.69385857270718176 + s *
        (-10.558626225323291 + s * (-62.375332450326006 + s *
        (-162.39666946257347 + s * (-184.60509290671104 + s * (-81.2874355063066
        + s * -9.8143293441691455))))));
      S = 1.0 + s * (19.651271667439257 + s * (137.65775414351904 + s *
        (434.56587747522923 + s * (645.38727173326788 + s * (429.00814002756783
        + s * (108.63500554177944 + s * (6.5702497703192817 + s *
        -0.0604244152148581)))))));
    } else {
      R = -0.0098649429247001 + s * (-0.799283237680523 + s *
        (-17.757954917754752 + s * (-160.63638485582192 + s *
        (-637.56644336838963 + s * (-1025.0951316110772 + s * -483.5191916086514)))));
      S = 1.0 + s * (30.338060743482458 + s * (325.79251299657392 + s *
        (1536.729586084437 + s * (3199.8582195085955 + s * (2553.0504064331644 +
        s * (474.52854120695537 + s * -22.440952446585818))))));
    }

    if ((!muDoubleScalarIsInf(absx)) && (!muDoubleScalarIsNaN(absx))) {
      s = frexp(absx, &eint);
    } else {
      s = absx;
      eint = 0;
    }

    s = muDoubleScalarFloor(s * 2.097152E+6) / 2.097152E+6 * muDoubleScalarPower
      (2.0, eint);
    y = muDoubleScalarExp(-s * s - 0.5625) * muDoubleScalarExp((s - absx) * (s +
      absx) + R / S) / absx;
    if (x < 0.0) {
      y = 2.0 - y;
    }
  }

  return y;
}

void b_erfc(const emlrtStack *sp, const emxArray_real_T *x, emxArray_real_T *y)
{
  int32_T i4;
  int32_T k;
  emlrtStack st;
  st.prev = sp;
  st.tls = sp->tls;
  st.site = &ec_emlrtRSI;
  for (i4 = 0; i4 < 3; i4++) {
    k = y->size[0] * y->size[1] * y->size[2];
    y->size[i4] = x->size[i4];
    emxEnsureCapacity(&st, (emxArray__common *)y, k, (int32_T)sizeof(real_T),
                      &m_emlrtRTEI);
  }

  i4 = x->size[0] * x->size[1] << 1;
  for (k = 0; k < i4; k++) {
    y->data[k] = scalar_erf(x->data[k]);
  }
}

/* End of code generation (erfc.c) */
