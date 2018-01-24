/*
 * GaussListND_mexCode.c
 *
 * Code generation for function 'GaussListND_mexCode'
 *
 * C source code generated on: Sun Dec  2 23:59:22 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "multiGaussListND_mexCode.h"
#include "GaussListND_mexCode.h"
#include "multiGaussListND_mexCode_emxutil.h"
#include "erfc.h"
#include "repmat.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void GaussListND_mexCode(const emxArray_real_T *coordList, const real_T
  sigma_data[3], const int32_T sigma_size[2], const real_T center_data[3], const
  int32_T center_size[2], emxArray_real_T *gaussList)
{
  emxArray_real_T *sigma2;
  emxArray_real_T *center2;
  real_T b_coordList[2];
  real_T c_coordList[2];
  int32_T stride;
  int32_T loop_ub;
  emxArray_real_T *coordList2;
  int32_T ix;
  int32_T ixstart;
  emxArray_real_T *varargin_2;
  int32_T ysize[3];
  int32_T sz1[2];
  emxArray_real_T *b_coordList2;
  int32_T iy;
  emxArray_real_T *c_coordList2;
  emxArray_real_T *b_y1;
  int32_T k;
  real_T p;
  int32_T exitg1;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  b_emxInit_real_T(&sigma2, 2, TRUE);
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
  /*                             Default: 0; */
  /*                             Rotation is only supported for 2D and 3D case */
  /*  */
  /*  OUTPUT   gaussList : m-by-1 list of intensities. Intensity is the */
  /*                       integral of the Gaussian over the pixel/voxel */
  /*  */
  /*  REMARKS  The code assumes that a pixel has the edge length 1! */
  /*  */
  /*  c: 2/05 jonas */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
  b_coordList[0] = (real_T)coordList->size[0];
  b_coordList[1] = 1.0;
  repmat(sigma_data, sigma_size, b_coordList, sigma2);
  c_coordList[0] = (real_T)coordList->size[0];
  c_coordList[1] = 1.0;
  repmat(center_data, center_size, c_coordList, center2);

  /* ====================== */
  /*  CALC GAUSSLIST */
  /* ====================== */
  /*  instead of calculating Gauss-values for very complicated geometries, we */
  /*  make a coordinate transformation so that we can use sigma=1 in all */
  /*  dimensions */
  /*  0.5*erfc(-(x+0.5)/sqrt(2))-0.5*erfc(-(x-0.5)/sqrt(2)) gives the integral on the */
  /*  pixel at 1 of a Gaussian with mean 0 and sigma 1 */
  /* center3 = center2(:,1); */
  /* clear center2 */
  /*  convert coordList to 0/1 */
  if (1 > coordList->size[0]) {
    stride = 0;
  } else {
    stride = coordList->size[0];
  }

  if (1 > coordList->size[1]) {
    loop_ub = 0;
  } else {
    loop_ub = coordList->size[1];
  }

  b_emxInit_real_T(&coordList2, 2, TRUE);
  ix = coordList2->size[0] * coordList2->size[1];
  coordList2->size[0] = stride;
  coordList2->size[1] = loop_ub;
  emxEnsureCapacity((emxArray__common *)coordList2, ix, (int32_T)sizeof(real_T));
  for (ix = 0; ix < loop_ub; ix++) {
    for (ixstart = 0; ixstart < stride; ixstart++) {
      coordList2->data[ixstart + coordList2->size[0] * ix] = (coordList->
        data[ixstart + coordList->size[0] * ix] - center2->data[ixstart +
        center2->size[0] * ix]) / sigma2->data[ixstart + sigma2->size[0] * ix];
    }
  }

  /* clear coordList center3 */
  /*  double coordList as preparation for erfc */
  /* fixed bug: must divide the 0.5 by sigma - KJ */
  ix = center2->size[0] * center2->size[1];
  center2->size[0] = coordList2->size[0];
  center2->size[1] = coordList2->size[1];
  emxEnsureCapacity((emxArray__common *)center2, ix, (int32_T)sizeof(real_T));
  stride = coordList2->size[1];
  for (ix = 0; ix < stride; ix++) {
    loop_ub = coordList2->size[0];
    for (ixstart = 0; ixstart < loop_ub; ixstart++) {
      center2->data[ixstart + center2->size[0] * ix] = coordList2->data[ixstart
        + coordList2->size[0] * ix] - 0.5 / sigma2->data[ixstart + sigma2->size
        [0] * ix];
    }
  }

  b_emxInit_real_T(&varargin_2, 2, TRUE);
  ix = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = coordList2->size[0];
  varargin_2->size[1] = coordList2->size[1];
  emxEnsureCapacity((emxArray__common *)varargin_2, ix, (int32_T)sizeof(real_T));
  stride = coordList2->size[1];
  for (ix = 0; ix < stride; ix++) {
    loop_ub = coordList2->size[0];
    for (ixstart = 0; ixstart < loop_ub; ixstart++) {
      varargin_2->data[ixstart + varargin_2->size[0] * ix] = coordList2->
        data[ixstart + coordList2->size[0] * ix] + 0.5 / sigma2->data[ixstart +
        sigma2->size[0] * ix];
    }
  }

  emxFree_real_T(&coordList2);
  for (ix = 0; ix < 3; ix++) {
    ysize[ix] = 1;
  }

  for (ix = 0; ix < 2; ix++) {
    sz1[ix] = center2->size[ix];
  }

  for (loop_ub = 0; loop_ub < 2; loop_ub++) {
    ysize[loop_ub] = sz1[loop_ub];
  }

  emxInit_real_T(&b_coordList2, 3, TRUE);
  ix = b_coordList2->size[0] * b_coordList2->size[1] * b_coordList2->size[2];
  b_coordList2->size[0] = ysize[0];
  b_coordList2->size[1] = ysize[1];
  b_coordList2->size[2] = 2;
  emxEnsureCapacity((emxArray__common *)b_coordList2, ix, (int32_T)sizeof(real_T));
  iy = -1;
  ix = center2->size[0] * center2->size[1];
  for (loop_ub = 1; loop_ub <= ix; loop_ub++) {
    iy++;
    b_coordList2->data[iy] = center2->data[loop_ub - 1];
  }

  emxFree_real_T(&center2);
  ix = varargin_2->size[0] * varargin_2->size[1];
  for (loop_ub = 1; loop_ub <= ix; loop_ub++) {
    iy++;
    b_coordList2->data[iy] = varargin_2->data[loop_ub - 1];
  }

  emxFree_real_T(&varargin_2);
  emxInit_real_T(&c_coordList2, 3, TRUE);

  /*  calculate gaussList */
  /* Jonas was missing the minus sign in erfc. I corrected that - KJ */
  ix = c_coordList2->size[0] * c_coordList2->size[1] * c_coordList2->size[2];
  c_coordList2->size[0] = b_coordList2->size[0];
  c_coordList2->size[1] = b_coordList2->size[1];
  c_coordList2->size[2] = 2;
  emxEnsureCapacity((emxArray__common *)c_coordList2, ix, (int32_T)sizeof(real_T));
  stride = b_coordList2->size[0] * b_coordList2->size[1] * b_coordList2->size[2];
  for (ix = 0; ix < stride; ix++) {
    c_coordList2->data[ix] = -b_coordList2->data[ix] / 1.4142135623730951;
  }

  b_erfc(c_coordList2, b_coordList2);
  ix = b_coordList2->size[0] * b_coordList2->size[1] * b_coordList2->size[2];
  b_coordList2->size[0] = b_coordList2->size[0];
  b_coordList2->size[1] = b_coordList2->size[1];
  b_coordList2->size[2] = 2;
  emxEnsureCapacity((emxArray__common *)b_coordList2, ix, (int32_T)sizeof(real_T));
  stride = b_coordList2->size[0];
  loop_ub = b_coordList2->size[1];
  ixstart = b_coordList2->size[2];
  stride = stride * loop_ub * ixstart;
  emxFree_real_T(&c_coordList2);
  for (ix = 0; ix < stride; ix++) {
    b_coordList2->data[ix] *= 0.5;
  }

  for (ix = 0; ix < 3; ix++) {
    ysize[ix] = b_coordList2->size[ix];
  }

  emxInit_real_T(&b_y1, 3, TRUE);
  ix = b_y1->size[0] * b_y1->size[1] * b_y1->size[2];
  b_y1->size[0] = ysize[0];
  b_y1->size[1] = ysize[1];
  b_y1->size[2] = 1;
  emxEnsureCapacity((emxArray__common *)b_y1, ix, (int32_T)sizeof(real_T));
  stride = 1;
  for (k = 0; k < 2; k++) {
    stride *= b_coordList2->size[k];
  }

  ix = 0;
  iy = 0;
  for (ixstart = 1; ixstart <= stride; ixstart++) {
    p = b_coordList2->data[ix];
    p = b_coordList2->data[ix + stride] - p;
    b_y1->data[iy] = p;
    ix++;
    iy++;
  }

  emxFree_real_T(&b_coordList2);
  for (ix = 0; ix < 3; ix++) {
    ysize[ix] = b_y1->size[ix];
  }

  ix = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
  gaussList->size[0] = ysize[0];
  gaussList->size[1] = 1;
  gaussList->size[2] = 1;
  emxEnsureCapacity((emxArray__common *)gaussList, ix, (int32_T)sizeof(real_T));
  if ((b_y1->size[0] == 0) || (b_y1->size[1] == 0)) {
    ix = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
    gaussList->size[0] = ysize[0];
    gaussList->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)gaussList, ix, (int32_T)sizeof(real_T));
    ix = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
    gaussList->size[2] = 1;
    emxEnsureCapacity((emxArray__common *)gaussList, ix, (int32_T)sizeof(real_T));
    stride = ysize[0];
    for (ix = 0; ix < stride; ix++) {
      gaussList->data[ix] = 1.0;
    }
  } else {
    iy = -1;
    ixstart = 0;
    for (loop_ub = 1; loop_ub <= b_y1->size[0]; loop_ub++) {
      ixstart++;
      ix = ixstart - 1;
      p = b_y1->data[ixstart - 1];
      for (k = 2; k <= b_y1->size[1]; k++) {
        ix += b_y1->size[0];
        p *= b_y1->data[ix];
      }

      iy++;
      gaussList->data[iy] = p;
    }
  }

  emxFree_real_T(&b_y1);

  /*  norm gaussList */
  ix = sigma2->size[1];
  if (ix == 0) {
    p = 1.0;
  } else {
    p = sigma2->data[0];
    k = 2;
    do {
      exitg1 = 0;
      ix = sigma2->size[1];
      if (k <= ix) {
        p *= sigma2->data[sigma2->size[0] * (k - 1)];
        k++;
      } else {
        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  emxFree_real_T(&sigma2);
  p *= muDoubleScalarPower(6.2831853071795862, 0.5 * (real_T)coordList->size[1]);
  ix = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
  gaussList->size[0] = gaussList->size[0];
  gaussList->size[1] = 1;
  gaussList->size[2] = gaussList->size[2];
  emxEnsureCapacity((emxArray__common *)gaussList, ix, (int32_T)sizeof(real_T));
  stride = gaussList->size[0];
  loop_ub = gaussList->size[1];
  ixstart = gaussList->size[2];
  stride = stride * loop_ub * ixstart;
  for (ix = 0; ix < stride; ix++) {
    gaussList->data[ix] *= p;
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (GaussListND_mexCode.c) */
