/*
 * GaussListND_mexCode.c
 *
 * Code generation for function 'GaussListND_mexCode'
 *
 * C source code generated on: Tue Nov 19 11:16:20 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "GaussListND_mexCode.h"
#include "fitNGaussians3D_mexCode_emxutil.h"
#include "cat.h"
#include "rdivide.h"
#include "repmat.h"

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

void GaussListND_mexCode(const emxArray_real_T *coordList, real_T sigma, real_T
  center, emxArray_real_T *gaussList)
{
  emxArray_real_T *center2;
  real_T b_coordList[2];
  int32_T loop_ub;
  emxArray_real_T *c_coordList;
  emxArray_real_T *b_center2;
  int32_T ixLead;
  int32_T iyLead;
  emxArray_real_T *coordList2;
  emxArray_real_T *b_coordList2;
  real_T tmp2;
  real_T tmp1;
  emxArray_real_T *c_coordList2;
  emxArray_real_T *d_coordList2;
  emxArray_real_T *y;
  int32_T k;
  int32_T ySize[3];
  int32_T stride;
  int32_T ix;
  int32_T iy;
  real_T work;
  int32_T vstride;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  b_emxInit_real_T(&center2, 2, TRUE);

  /* GAUSSLISTND calculates the value of a N-D Gaussian at specific pixel/voxel coordinates */
  /*  */
  /*  SYNOPSIS gaussList = GaussListND(coordList,sigma,center,intNorm,rotation) */
  /*  */
  /*  INPUT    coordList : m-by-n list of coordinates, where m is the number of */
  /*                       coordinates and n the number of dimensions */
  /*           sigma     : 1-by-n (or scalar): sigma of Gaussian */
  /*           center    : (opt) 1-by-n vector of center of Gaussian. */
  /*                       Default: zeros(1,n) */
  /*           intNorm   : (opt) switch for how the Gaussian should be normed */
  /*                       Default: 0 */
  /*                       0 - no norming. Max of Gaussian == 1 */
  /*                       1 - normed so that integral of infinite Gaussian = 1 */
  /*           rotation  : (opt) Equal to the number of degree you want the */
  /*                             coordinate to be rotate for. If rotation is */
  /*                             equal to 1, rotation will be random. */
  /*                             Default: 0;  */
  /*                             Rotation is only supported for 2D and 3D case */
  /*  */
  /*  OUTPUT   gaussList : m-by-1 list of intensities. Intensity is the */
  /*                       integral of the Gaussian over the pixel/voxel */
  /*  */
  /*  REMARKS  The code assumes that a pixel has the edge length 1! */
  /*  */
  /*  c: 2/05 jonas */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  /* ====================== */
  /*  TEST INPUT */
  /* ====================== */
  /*  check number of input arguments */
  /*  the following doesn't work with Matlab 6.5.0 */
  /*  error(nargchk(2,4,nIn,'struct')); */
  /*  % check dimensionality of coordList. */
  /*  if isempty(coordList) */
  /*      error('you have to supply a list of coordinates for GaussList23D') */
  /*  else */
  /*      [nCoords,nDims] = size(coordList); */
  /*  end */
  /*  % sigma */
  /*  ls = length(sigma); */
  /*  switch ls */
  /*      case nDims */
  /*          % make as long as coords */
  /*          sigma = repmat(sigma,[nCoords,1]); */
  /*      case 1 */
  /*          sigma = repmat(sigma,[nCoords,nDims]); */
  /*      otherwise */
  /*          error('sigma has to be a scalar or a 1-by-n vector!') */
  /*  end */
  /*  center */
  /*  if nIn < 3 || isempty(center) */
  /*      center = zeros(nCoords,nDims); */
  /*  else */
  /*      lc = length(center); */
  /*      switch lc */
  /*          case nDims */
  /*              center = repmat(center, [nCoords,1]); */
  /*          case 1 */
  /*              center = repmat(center, [nCoords,3]); */
  /*          otherwise */
  /*              error('center has to be a scalar or a 1-by-n vector!') */
  /*      end */
  /*  end */
  b_coordList[0] = coordList->size[0];
  b_coordList[1] = 1.0;
  repmat(center, b_coordList, center2);

  /* clear center */
  /*  intNorm */
  /* rotation */
  /*  coordDim = size(coordList,2); */
  /*  if nIn < 5 || isempty(rotation) || rotation == 0 */
  /*      rotation = 0; */
  /*      alp = 0; */
  /*      bet = 0; */
  /*      delt = 0; */
  /*  elseif rotation == 1 && coordDim <= 2 */
  /*      rotation = floor(rand(1) * 360); */
  /*  elseif rotation == 1 && coordDim == 3 */
  /*      alp = floor(rand(1) * 180); */
  /*      bet = floor(rand(1) * 180); */
  /*      delt = floor(rand(1) * 360); */
  /*  end */
  /*  if rotation && (nDims < 2 || nDims > 3) */
  /*      error('rotation is only supported for 2-3 dimensions') */
  /*  end */
  /* ====================== */
  /* ====================== */
  /*  CALC GAUSSLIST */
  /* ====================== */
  /*  instead of calculating Gauss-values for very complicated geometries, we */
  /*  make a coordinate transformation so that we can use sigma=1 in all */
  /*  dimensions */
  /*  if rotation ~= 0 */
  /*       */
  /*      %Translate center to origin. */
  /*      coordList = coordList - center2; */
  /*       */
  /*      if coordDim == 2 */
  /*           */
  /*          % 2 Dimension rotation. */
  /*          %Rotation. */
  /*          %Rotation of the coordinate. x' = xcos@ - ysin@. y' = xsin@ + ycos@. */
  /*          tmpX = coordList(:,1) .* cosd(rotation) - coordList(:,2) .* sind(rotation); */
  /*          tmpY = coordList(:,1) .* sind(rotation) + coordList(:,2) .* cosd(rotation); */
  /*           */
  /*          %Translation back to original center. */
  /*          coordList(:,1) = tmpX(:,1) + center2(:,1); */
  /*          coordList(:,2) = tmpY(:,1) + center2(:,2); */
  /*           */
  /*      elseif coordDim == 3 */
  /*           */
  /*          % 3 Dimension rotation. */
  /*          %Rotation of the coordinate. */
  /*          c1 = cos(alp); c2 = cos(bet); c3 = cos(delt); */
  /*          s1 = sin(alp); s2 = sin(bet); s3 = sin(delt); */
  /*           */
  /*          l1 = (c2 * c3) - (c1*s2*s3); l2 = -(c2 * s3) - (c1 * s2 * c3); */
  /*          l3 = s1*s2; */
  /*          m1 = (s2*s3 + c1*c2*s3); m2 = -(s2*s3) + (c1*c2*c3); */
  /*          m3 = -(s1*c2); */
  /*          n1 = s1*s3; n2 = s1*c3; n3 = c1; */
  /*          %Calculation of my new coordinate in function of the rotation. */
  /*          tmpX = coordList(:,1) .* l1 + coordList(:,2) .* l2 + coordList(:,3) .* l3; */
  /*          tmpY = coordList(:,1) .* m1 + coordList(:,2) .* m2 + coordList(:,3) .* m3; */
  /*          tmpZ = coordList(:,1) .* n1 + coordList(:,3) .* n2 + coordList(:,3) .* n3; */
  /*           */
  /*          %Translation back to original center - KJ addition to make */
  /*          %consistent with 2D case */
  /*          %otherwise the code returns nonsense */
  /*          coordList(:,1) = tmpX(:,1) + center2(:,1); */
  /*          coordList(:,2) = tmpY(:,1) + center2(:,2); */
  /*          coordList(:,3) = tmpZ(:,1) + center2(:,3); */
  /*           */
  /*      end */
  /*    */
  /*       */
  /*  end */
  /*  0.5*erfc(-(x+0.5)/sqrt(2))-0.5*erfc(-(x-0.5)/sqrt(2)) gives the integral on the */
  /*  pixel at 1 of a Gaussian with mean 0 and sigma 1 */
  /* clear center2 */
  /*  convert coordList to 0/1 */
  if (1 > coordList->size[0]) {
    loop_ub = 0;
  } else {
    loop_ub = coordList->size[0];
  }

  emxInit_real_T(&c_coordList, 1, TRUE);
  emxInit_real_T(&b_center2, 1, TRUE);
  ixLead = center2->size[0];
  iyLead = b_center2->size[0];
  b_center2->size[0] = ixLead;
  emxEnsureCapacity((emxArray__common *)b_center2, iyLead, (int32_T)sizeof
                    (real_T));
  for (iyLead = 0; iyLead < ixLead; iyLead++) {
    b_center2->data[iyLead] = center2->data[iyLead];
  }

  emxFree_real_T(&center2);
  iyLead = c_coordList->size[0];
  c_coordList->size[0] = loop_ub;
  emxEnsureCapacity((emxArray__common *)c_coordList, iyLead, (int32_T)sizeof
                    (real_T));
  for (iyLead = 0; iyLead < loop_ub; iyLead++) {
    c_coordList->data[iyLead] = coordList->data[iyLead] - b_center2->data[iyLead];
  }

  emxFree_real_T(&b_center2);
  emxInit_real_T(&coordList2, 1, TRUE);
  emxInit_real_T(&b_coordList2, 1, TRUE);
  rdivide(c_coordList, sigma, coordList2);

  /* clear coordList center3 */
  /*  double coordList as preparation for erfc */
  /* fixed bug: must divide the 0.5 by sigma - KJ */
  tmp2 = 0.5 / sigma;
  tmp1 = 0.5 / sigma;
  iyLead = b_coordList2->size[0];
  b_coordList2->size[0] = coordList2->size[0];
  emxEnsureCapacity((emxArray__common *)b_coordList2, iyLead, (int32_T)sizeof
                    (real_T));
  loop_ub = coordList2->size[0];
  emxFree_real_T(&c_coordList);
  for (iyLead = 0; iyLead < loop_ub; iyLead++) {
    b_coordList2->data[iyLead] = coordList2->data[iyLead] - tmp2;
  }

  emxInit_real_T(&c_coordList2, 1, TRUE);
  iyLead = c_coordList2->size[0];
  c_coordList2->size[0] = coordList2->size[0];
  emxEnsureCapacity((emxArray__common *)c_coordList2, iyLead, (int32_T)sizeof
                    (real_T));
  loop_ub = coordList2->size[0];
  for (iyLead = 0; iyLead < loop_ub; iyLead++) {
    c_coordList2->data[iyLead] = coordList2->data[iyLead] + tmp1;
  }

  emxFree_real_T(&coordList2);
  c_emxInit_real_T(&d_coordList2, 3, TRUE);
  cat(b_coordList2, c_coordList2, d_coordList2);

  /*  calculate gaussList */
  /* Jonas was missing the minus sign in erfc. I corrected that - KJ */
  iyLead = d_coordList2->size[0] * d_coordList2->size[1] * d_coordList2->size[2];
  emxEnsureCapacity((emxArray__common *)d_coordList2, iyLead, (int32_T)sizeof
                    (real_T));
  loop_ub = d_coordList2->size[0];
  ixLead = d_coordList2->size[1];
  iyLead = d_coordList2->size[2];
  loop_ub = loop_ub * ixLead * iyLead;
  emxFree_real_T(&c_coordList2);
  emxFree_real_T(&b_coordList2);
  for (iyLead = 0; iyLead < loop_ub; iyLead++) {
    d_coordList2->data[iyLead] = -d_coordList2->data[iyLead] /
      1.4142135623730951;
  }

  c_emxInit_real_T(&y, 3, TRUE);
  for (iyLead = 0; iyLead < 3; iyLead++) {
    loop_ub = y->size[0] * y->size[1] * y->size[2];
    y->size[iyLead] = d_coordList2->size[iyLead];
    emxEnsureCapacity((emxArray__common *)y, loop_ub, (int32_T)sizeof(real_T));
  }

  iyLead = d_coordList2->size[0] * d_coordList2->size[1] * d_coordList2->size[2];
  for (k = 0; k < iyLead; k++) {
    y->data[k] = scalar_erf(d_coordList2->data[k]);
  }

  iyLead = y->size[0] * y->size[1] * y->size[2];
  emxEnsureCapacity((emxArray__common *)y, iyLead, (int32_T)sizeof(real_T));
  loop_ub = y->size[0];
  ixLead = y->size[1];
  iyLead = y->size[2];
  loop_ub = loop_ub * ixLead * iyLead;
  for (iyLead = 0; iyLead < loop_ub; iyLead++) {
    y->data[iyLead] *= 0.5;
  }

  if (y->size[2] <= 1) {
    for (iyLead = 0; iyLead < 3; iyLead++) {
      ySize[iyLead] = y->size[iyLead];
    }

    iyLead = d_coordList2->size[0] * d_coordList2->size[1] * d_coordList2->size
      [2];
    d_coordList2->size[0] = ySize[0];
    emxEnsureCapacity((emxArray__common *)d_coordList2, iyLead, (int32_T)sizeof
                      (real_T));
    iyLead = d_coordList2->size[0] * d_coordList2->size[1] * d_coordList2->size
      [2];
    d_coordList2->size[1] = ySize[1];
    d_coordList2->size[2] = 0;
    emxEnsureCapacity((emxArray__common *)d_coordList2, iyLead, (int32_T)sizeof
                      (real_T));
  } else {
    for (iyLead = 0; iyLead < 3; iyLead++) {
      ySize[iyLead] = y->size[iyLead];
    }

    iyLead = d_coordList2->size[0] * d_coordList2->size[1] * d_coordList2->size
      [2];
    d_coordList2->size[0] = ySize[0];
    d_coordList2->size[1] = ySize[1];
    d_coordList2->size[2] = y->size[2] - 1;
    emxEnsureCapacity((emxArray__common *)d_coordList2, iyLead, (int32_T)sizeof
                      (real_T));
    loop_ub = y->size[2] - 1;
    if (!((ySize[0] == 0) || (ySize[1] == 0) || (loop_ub == 0))) {
      k = 3;
      while ((k > 2) && (y->size[2] == 1)) {
        k = 2;
      }

      if (3 > k) {
        stride = y->size[0] * y->size[1] * y->size[2];
      } else {
        stride = 1;
        for (k = 0; k < 2; k++) {
          iyLead = y->size[k];
          stride *= iyLead;
        }
      }

      ix = 0;
      iy = 0;
      for (loop_ub = 1; loop_ub <= stride; loop_ub++) {
        ixLead = ix + stride;
        iyLead = iy;
        work = y->data[ix];
        for (vstride = 2; vstride <= y->size[2]; vstride++) {
          tmp2 = work;
          work = y->data[ixLead];
          tmp1 = y->data[ixLead] - tmp2;
          ixLead += stride;
          d_coordList2->data[iyLead] = tmp1;
          iyLead += stride;
        }

        ix++;
        iy++;
      }
    }
  }

  emxFree_real_T(&y);
  for (iyLead = 0; iyLead < 3; iyLead++) {
    ySize[iyLead] = d_coordList2->size[iyLead];
  }

  iyLead = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
  gaussList->size[0] = ySize[0];
  gaussList->size[1] = 1;
  gaussList->size[2] = ySize[2];
  emxEnsureCapacity((emxArray__common *)gaussList, iyLead, (int32_T)sizeof
                    (real_T));
  if ((d_coordList2->size[0] == 0) || (d_coordList2->size[1] == 0) ||
      (d_coordList2->size[2] == 0)) {
    iyLead = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
    gaussList->size[0] = ySize[0];
    gaussList->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)gaussList, iyLead, (int32_T)sizeof
                      (real_T));
    iyLead = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
    gaussList->size[2] = ySize[2];
    emxEnsureCapacity((emxArray__common *)gaussList, iyLead, (int32_T)sizeof
                      (real_T));
    loop_ub = ySize[0] * ySize[2];
    for (iyLead = 0; iyLead < loop_ub; iyLead++) {
      gaussList->data[iyLead] = 1.0;
    }
  } else {
    vstride = d_coordList2->size[0];
    stride = 1;
    k = 3;
    while ((k > 2) && (d_coordList2->size[2] == 1)) {
      k = 2;
    }

    loop_ub = 3;
    while (loop_ub <= k) {
      iyLead = d_coordList2->size[2];
      stride *= iyLead;
      loop_ub = 4;
    }

    ix = 0;
    iy = -1;
    for (loop_ub = 1; loop_ub <= stride; loop_ub++) {
      ixLead = ix;
      for (iyLead = 1; iyLead <= vstride; iyLead++) {
        ixLead++;
        ix = ixLead;
        tmp2 = d_coordList2->data[ixLead - 1];
        for (k = 2; k <= d_coordList2->size[1]; k++) {
          ix += vstride;
          tmp2 *= d_coordList2->data[ix - 1];
        }

        iy++;
        gaussList->data[iy] = tmp2;
      }
    }
  }

  emxFree_real_T(&d_coordList2);

  /*  norm gaussList */
  /*  "un-norm" Gaussian */
  tmp2 = 2.5066282746310002 * sigma;
  iyLead = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
  gaussList->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)gaussList, iyLead, (int32_T)sizeof
                    (real_T));
  loop_ub = gaussList->size[0];
  ixLead = gaussList->size[1];
  iyLead = gaussList->size[2];
  loop_ub = loop_ub * ixLead * iyLead;
  for (iyLead = 0; iyLead < loop_ub; iyLead++) {
    gaussList->data[iyLead] *= tmp2;
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (GaussListND_mexCode.c) */
