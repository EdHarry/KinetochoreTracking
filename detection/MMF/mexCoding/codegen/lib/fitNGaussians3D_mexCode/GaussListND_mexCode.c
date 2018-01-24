/*
 * GaussListND_mexCode.c
 *
 * Code generation for function 'GaussListND_mexCode'
 *
 * C source code generated on: Thu May  3 12:56:20 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "GaussListND_mexCode.h"
#include "fitNGaussians3D_mexCode_emxutil.h"
#include "rdivide.h"
#include "fitNGaussians3D_mexCode_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void GaussListND_mexCode(const emxArray_real_T *coordList, real_T sigma, real_T
  center, emxArray_real_T *gaussList)
{
  int32_T nCoords;
  int32_T b_nCoords[2];
  int32_T outsize[2];
  int32_T vstride;
  emxArray_real_T *b_coordList;
  emxArray_real_T *b_center;
  int32_T stride;
  emxArray_real_T *coordList2;
  emxArray_real_T *varargin_1;
  real_T z;
  real_T p;
  int32_T ysize[3];
  int32_T sz1[2];
  int32_T j;
  emxArray_real_T *b_coordList2;
  int32_T iy;
  emxArray_real_T *y;
  int32_T k;
  real_T absx;
  real_T R;
  real_T S;
  int32_T ix;
  emxArray_real_T *b_y1;

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
  /*  check dimensionality of coordList. */
  if (coordList->size[0] == 0) {
  } else {
    nCoords = coordList->size[0];
  }

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
  b_nCoords[0] = nCoords;
  b_nCoords[1] = 1;
  for (vstride = 0; vstride < 2; vstride++) {
    outsize[vstride] = b_nCoords[vstride];
  }

  c_emxInit_real_T(&b_coordList, 1);
  c_emxInit_real_T(&b_center, 1);

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
  vstride = b_center->size[0];
  b_center->size[0] = outsize[0];
  emxEnsureCapacity((emxArray__common *)b_center, vstride, (int32_T)sizeof
                    (real_T));
  stride = outsize[0] - 1;
  for (vstride = 0; vstride <= stride; vstride++) {
    b_center->data[vstride] = center;
  }

  vstride = b_coordList->size[0];
  b_coordList->size[0] = nCoords;
  emxEnsureCapacity((emxArray__common *)b_coordList, vstride, (int32_T)sizeof
                    (real_T));
  stride = nCoords - 1;
  for (vstride = 0; vstride <= stride; vstride++) {
    b_coordList->data[vstride] = coordList->data[vstride] - b_center->
      data[vstride];
  }

  emxFree_real_T(&b_center);
  c_emxInit_real_T(&coordList2, 1);
  c_emxInit_real_T(&varargin_1, 1);
  rdivide(b_coordList, sigma, coordList2);

  /* clear coordList center3 */
  /*  double coordList as preparation for erfc */
  /* fixed bug: must divide the 0.5 by sigma - KJ */
  z = 0.5 / sigma;
  p = 0.5 / sigma;
  vstride = varargin_1->size[0];
  varargin_1->size[0] = coordList2->size[0];
  emxEnsureCapacity((emxArray__common *)varargin_1, vstride, (int32_T)sizeof
                    (real_T));
  emxFree_real_T(&b_coordList);
  stride = coordList2->size[0] - 1;
  for (vstride = 0; vstride <= stride; vstride++) {
    varargin_1->data[vstride] = coordList2->data[vstride] - z;
  }

  vstride = coordList2->size[0];
  coordList2->size[0] = coordList2->size[0];
  emxEnsureCapacity((emxArray__common *)coordList2, vstride, (int32_T)sizeof
                    (real_T));
  stride = coordList2->size[0] - 1;
  for (vstride = 0; vstride <= stride; vstride++) {
    coordList2->data[vstride] += p;
  }

  for (vstride = 0; vstride < 3; vstride++) {
    ysize[vstride] = 1;
  }

  sz1[0] = varargin_1->size[0];
  sz1[1] = 1;
  for (j = 0; j < 2; j++) {
    ysize[j] = sz1[j];
  }

  b_emxInit_real_T(&b_coordList2, 3);
  ysize[2] = 2;
  vstride = b_coordList2->size[0] * b_coordList2->size[1] * b_coordList2->size[2];
  b_coordList2->size[0] = ysize[0];
  b_coordList2->size[1] = ysize[1];
  b_coordList2->size[2] = 2;
  emxEnsureCapacity((emxArray__common *)b_coordList2, vstride, (int32_T)sizeof
                    (real_T));
  iy = -1;
  for (j = 1; j <= varargin_1->size[0]; j++) {
    iy++;
    b_coordList2->data[iy] = varargin_1->data[j - 1];
  }

  emxFree_real_T(&varargin_1);
  for (j = 1; j <= coordList2->size[0]; j++) {
    iy++;
    b_coordList2->data[iy] = coordList2->data[j - 1];
  }

  emxFree_real_T(&coordList2);

  /*  calculate gaussList */
  /* Jonas was missing the minus sign in erfc. I corrected that - KJ */
  vstride = b_coordList2->size[0] * b_coordList2->size[1] * b_coordList2->size[2];
  b_coordList2->size[0] = b_coordList2->size[0];
  b_coordList2->size[1] = b_coordList2->size[1];
  b_coordList2->size[2] = 2;
  emxEnsureCapacity((emxArray__common *)b_coordList2, vstride, (int32_T)sizeof
                    (real_T));
  vstride = b_coordList2->size[0];
  nCoords = b_coordList2->size[1];
  stride = b_coordList2->size[2];
  stride = vstride * nCoords * stride - 1;
  for (vstride = 0; vstride <= stride; vstride++) {
    b_coordList2->data[vstride] = -b_coordList2->data[vstride] /
      1.4142135623730951;
  }

  b_emxInit_real_T(&y, 3);
  for (vstride = 0; vstride < 3; vstride++) {
    nCoords = y->size[0] * y->size[1] * y->size[2];
    y->size[vstride] = b_coordList2->size[vstride];
    emxEnsureCapacity((emxArray__common *)y, nCoords, (int32_T)sizeof(real_T));
  }

  vstride = b_coordList2->size[0] * b_coordList2->size[1] << 1;
  for (k = 0; k <= vstride - 1; k++) {
    /* <LEGAL> The algorithms for calculating ERF(X) and ERFC(X) are derived */
    /* <LEGAL> from FDLIBM, which has the following notice: */
    /* <LEGAL> */
    /* <LEGAL> Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved. */
    /* <LEGAL> */
    /* <LEGAL> Developed at SunSoft, a Sun Microsystems, Inc. business. */
    /* <LEGAL> Permission to use, copy, modify, and distribute this */
    /* <LEGAL> software is freely granted, provided that this notice */
    /* <LEGAL> is preserved. */
    absx = fabs(b_coordList2->data[k]);
    if (rtIsNaN(b_coordList2->data[k])) {
      p = b_coordList2->data[k];
    } else if (rtIsInf(b_coordList2->data[k])) {
      if (b_coordList2->data[k] < 0.0) {
        p = 2.0;
      } else {
        p = 0.0;
      }
    } else if (absx < 0.84375) {
      if (absx < 1.3877787807814457E-17) {
        p = 1.0 - b_coordList2->data[k];
      } else {
        z = b_coordList2->data[k] * b_coordList2->data[k];
        p = (0.12837916709551256 + z * (-0.3250421072470015 + z *
              (-0.02848174957559851 + z * (-0.0057702702964894416 + z *
                -2.3763016656650163E-5)))) / (1.0 + z * (0.39791722395915535 + z
          * (0.0650222499887673 + z * (0.0050813062818757656 + z *
          (0.00013249473800432164 + z * -3.9602282787753681E-6)))));
        if (b_coordList2->data[k] < 0.25) {
          p = 1.0 - (b_coordList2->data[k] + b_coordList2->data[k] * p);
        } else {
          p = 0.5 - (b_coordList2->data[k] * p + (b_coordList2->data[k] - 0.5));
        }
      }
    } else if (absx < 1.25) {
      if (b_coordList2->data[k] >= 0.0) {
        p = 0.15493708848953247 - (-0.0023621185607526594 + (absx - 1.0) *
          (0.41485611868374833 + (absx - 1.0) * (-0.37220787603570132 + (absx -
          1.0) * (0.31834661990116175 + (absx - 1.0) * (-0.11089469428239668 +
          (absx - 1.0) * (0.035478304325618236 + (absx - 1.0) *
                          -0.0021663755948687908)))))) / (1.0 + (absx - 1.0) *
          (0.10642088040084423 + (absx - 1.0) * (0.540397917702171 + (absx - 1.0)
          * (0.071828654414196266 + (absx - 1.0) * (0.12617121980876164 + (absx
          - 1.0) * (0.013637083912029051 + (absx - 1.0) * 0.011984499846799107))))));
      } else {
        p = 1.0 + (0.84506291151046753 + (-0.0023621185607526594 + (absx - 1.0) *
                    (0.41485611868374833 + (absx - 1.0) * (-0.37220787603570132
          + (absx - 1.0) * (0.31834661990116175 + (absx - 1.0) *
                            (-0.11089469428239668 + (absx - 1.0) *
                             (0.035478304325618236 + (absx - 1.0) *
                              -0.0021663755948687908)))))) / (1.0 + (absx - 1.0)
                    * (0.10642088040084423 + (absx - 1.0) * (0.540397917702171 +
                      (absx - 1.0) * (0.071828654414196266 + (absx - 1.0) *
          (0.12617121980876164 + (absx - 1.0) * (0.013637083912029051 + (absx -
          1.0) * 0.011984499846799107)))))));
      }
    } else if (b_coordList2->data[k] < -6.0) {
      p = 2.0;
    } else if (b_coordList2->data[k] >= 28.0) {
      p = 0.0;
    } else {
      p = 1.0 / (absx * absx);
      if (absx < 2.8571414947509766) {
        R = -0.0098649440348471482 + p * (-0.69385857270718176 + p *
          (-10.558626225323291 + p * (-62.375332450326006 + p *
          (-162.39666946257347 + p * (-184.60509290671104 + p *
          (-81.2874355063066 + p * -9.8143293441691455))))));
        S = 1.0 + p * (19.651271667439257 + p * (137.65775414351904 + p *
          (434.56587747522923 + p * (645.38727173326788 + p *
          (429.00814002756783 + p * (108.63500554177944 + p *
          (6.5702497703192817 + p * -0.0604244152148581)))))));
      } else {
        R = -0.0098649429247001 + p * (-0.799283237680523 + p *
          (-17.757954917754752 + p * (-160.63638485582192 + p *
          (-637.56644336838963 + p * (-1025.0951316110772 + p *
          -483.5191916086514)))));
        S = 1.0 + p * (30.338060743482458 + p * (325.79251299657392 + p *
          (1536.729586084437 + p * (3199.8582195085955 + p * (2553.0504064331644
          + p * (474.52854120695537 + p * -22.440952446585818))))));
      }

      if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
        p = frexp(absx, &ix);
        nCoords = ix;
      } else {
        p = absx;
        nCoords = 0;
      }

      z = floor(p * 2.097152E+6) / 2.097152E+6 * rt_powd_snf(2.0, (real_T)
        nCoords);
      p = exp(-z * z - 0.5625) * exp((z - absx) * (z + absx) + R / S) / absx;
      if (b_coordList2->data[k] < 0.0) {
        p = 2.0 - p;
      }
    }

    y->data[k] = p;
  }

  emxFree_real_T(&b_coordList2);
  vstride = y->size[0] * y->size[1] * y->size[2];
  y->size[0] = y->size[0];
  y->size[1] = y->size[1];
  y->size[2] = 2;
  emxEnsureCapacity((emxArray__common *)y, vstride, (int32_T)sizeof(real_T));
  nCoords = y->size[0];
  stride = y->size[1];
  vstride = y->size[2];
  stride = nCoords * stride * vstride - 1;
  for (vstride = 0; vstride <= stride; vstride++) {
    y->data[vstride] *= 0.5;
  }

  for (vstride = 0; vstride < 3; vstride++) {
    ysize[vstride] = y->size[vstride];
  }

  b_emxInit_real_T(&b_y1, 3);
  ysize[2] = 1;
  vstride = b_y1->size[0] * b_y1->size[1] * b_y1->size[2];
  b_y1->size[0] = ysize[0];
  b_y1->size[1] = ysize[1];
  b_y1->size[2] = 1;
  emxEnsureCapacity((emxArray__common *)b_y1, vstride, (int32_T)sizeof(real_T));
  stride = 1;
  for (k = 0; k < 2; k++) {
    stride *= y->size[k];
  }

  ix = 0;
  iy = 0;
  for (nCoords = 1; nCoords <= stride; nCoords++) {
    p = y->data[ix];
    p = y->data[ix + stride] - p;
    b_y1->data[iy] = p;
    ix++;
    iy++;
  }

  emxFree_real_T(&y);
  for (vstride = 0; vstride < 3; vstride++) {
    ysize[vstride] = b_y1->size[vstride];
  }

  ysize[1] = 1;
  vstride = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
  gaussList->size[0] = ysize[0];
  gaussList->size[1] = 1;
  gaussList->size[2] = 1;
  emxEnsureCapacity((emxArray__common *)gaussList, vstride, (int32_T)sizeof
                    (real_T));
  nCoords = b_y1->size[1];
  vstride = b_y1->size[0];
  iy = -1;
  stride = -1;
  for (j = 1; j <= vstride; j++) {
    stride++;
    ix = stride;
    p = b_y1->data[stride];
    for (k = 2; k <= nCoords; k++) {
      ix += vstride;
      p *= b_y1->data[ix];
    }

    iy++;
    gaussList->data[iy] = p;
  }

  emxFree_real_T(&b_y1);

  /*  norm gaussList */
  /*  "un-norm" Gaussian */
  p = 2.5066282746310002 * sigma;
  vstride = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
  gaussList->size[0] = gaussList->size[0];
  gaussList->size[1] = 1;
  gaussList->size[2] = 1;
  emxEnsureCapacity((emxArray__common *)gaussList, vstride, (int32_T)sizeof
                    (real_T));
  nCoords = gaussList->size[0];
  stride = gaussList->size[1];
  vstride = gaussList->size[2];
  stride = nCoords * stride * vstride - 1;
  for (vstride = 0; vstride <= stride; vstride++) {
    gaussList->data[vstride] *= p;
  }
}

/* End of code generation (GaussListND_mexCode.c) */
