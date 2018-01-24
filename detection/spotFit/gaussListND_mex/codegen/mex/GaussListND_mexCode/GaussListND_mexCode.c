/*
 * GaussListND_mexCode.c
 *
 * Code generation for function 'GaussListND_mexCode'
 *
 * C source code generated on: Sun Dec  2 22:57:02 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GaussListND_mexCode.h"
#include "GaussListND_mexCode_emxutil.h"
#include "erfc.h"
#include "rdivide.h"
#include "repmat.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void GaussListND_mexCode(const emxArray_real_T *coordList, const emxArray_real_T
  *sigma, const emxArray_real_T *center, emxArray_real_T *gaussList)
{
  emxArray_real_T *sigma2;
  emxArray_real_T *center2;
  real_T b_coordList[2];
  real_T c_coordList[2];
  int32_T loop_ub;
  int32_T ix;
  emxArray_real_T *coordList2;
  int32_T i0;
  int32_T iy;
  int32_T stride;
  int32_T ixstart;
  emxArray_real_T *b_sigma2;
  emxArray_real_T *c_sigma2;
  emxArray_real_T *r0;
  int32_T ysize[3];
  uint32_T sz1[2];
  emxArray_real_T *b_coordList2;
  emxArray_real_T *c_coordList2;
  emxArray_real_T *b_y1;
  int32_T k;
  real_T p;
  emxArray_int32_T *r1;
  emxArray_int32_T *r2;
  int32_T exitg1;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  emxInit_real_T(&sigma2, 2, TRUE);
  emxInit_real_T(&center2, 2, TRUE);

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
  repmat(sigma, b_coordList, sigma2);
  c_coordList[0] = (real_T)coordList->size[0];
  c_coordList[1] = 1.0;
  b_repmat(center, c_coordList, center2);

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
    loop_ub = 0;
  } else {
    loop_ub = coordList->size[0];
  }

  if (1 > coordList->size[1]) {
    ix = 0;
  } else {
    ix = coordList->size[1];
  }

  emxInit_real_T(&coordList2, 2, TRUE);
  i0 = coordList2->size[0] * coordList2->size[1];
  coordList2->size[0] = loop_ub;
  coordList2->size[1] = ix;
  emxEnsureCapacity((emxArray__common *)coordList2, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < ix; i0++) {
    for (iy = 0; iy < loop_ub; iy++) {
      coordList2->data[iy + coordList2->size[0] * i0] = (coordList->data[iy +
        coordList->size[0] * i0] - center2->data[iy + center2->size[0] * i0]) /
        sigma2->data[iy + sigma2->size[0] * i0];
    }
  }

  /* clear coordList center3 */
  /*  double coordList as preparation for erfc */
  /* fixed bug: must divide the 0.5 by sigma - KJ */
  if (1 > coordList->size[0]) {
    loop_ub = 0;
  } else {
    loop_ub = coordList->size[0];
  }

  if (1 > coordList->size[1]) {
    ix = 0;
  } else {
    ix = coordList->size[1];
  }

  if (1 > coordList->size[0]) {
    stride = 0;
  } else {
    stride = coordList->size[0];
  }

  if (1 > coordList->size[1]) {
    ixstart = 0;
  } else {
    ixstart = coordList->size[1];
  }

  emxInit_real_T(&b_sigma2, 2, TRUE);
  i0 = b_sigma2->size[0] * b_sigma2->size[1];
  b_sigma2->size[0] = loop_ub;
  b_sigma2->size[1] = ix;
  emxEnsureCapacity((emxArray__common *)b_sigma2, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < ix; i0++) {
    for (iy = 0; iy < loop_ub; iy++) {
      b_sigma2->data[iy + b_sigma2->size[0] * i0] = sigma2->data[iy +
        sigma2->size[0] * i0];
    }
  }

  rdivide(0.5, b_sigma2, center2);
  i0 = center2->size[0] * center2->size[1];
  center2->size[0] = coordList2->size[0];
  center2->size[1] = coordList2->size[1];
  emxEnsureCapacity((emxArray__common *)center2, i0, (int32_T)sizeof(real_T));
  loop_ub = coordList2->size[0] * coordList2->size[1];
  emxFree_real_T(&b_sigma2);
  for (i0 = 0; i0 < loop_ub; i0++) {
    center2->data[i0] = coordList2->data[i0] - center2->data[i0];
  }

  emxInit_real_T(&c_sigma2, 2, TRUE);
  i0 = c_sigma2->size[0] * c_sigma2->size[1];
  c_sigma2->size[0] = stride;
  c_sigma2->size[1] = ixstart;
  emxEnsureCapacity((emxArray__common *)c_sigma2, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < ixstart; i0++) {
    for (iy = 0; iy < stride; iy++) {
      c_sigma2->data[iy + c_sigma2->size[0] * i0] = sigma2->data[iy +
        sigma2->size[0] * i0];
    }
  }

  emxInit_real_T(&r0, 2, TRUE);
  rdivide(0.5, c_sigma2, r0);
  i0 = coordList2->size[0] * coordList2->size[1];
  coordList2->size[0] = coordList2->size[0];
  coordList2->size[1] = coordList2->size[1];
  emxEnsureCapacity((emxArray__common *)coordList2, i0, (int32_T)sizeof(real_T));
  ixstart = coordList2->size[0];
  stride = coordList2->size[1];
  loop_ub = ixstart * stride;
  emxFree_real_T(&c_sigma2);
  for (i0 = 0; i0 < loop_ub; i0++) {
    coordList2->data[i0] += r0->data[i0];
  }

  emxFree_real_T(&r0);
  for (i0 = 0; i0 < 3; i0++) {
    ysize[i0] = 1;
  }

  for (i0 = 0; i0 < 2; i0++) {
    sz1[i0] = (uint32_T)center2->size[i0];
  }

  for (stride = 0; stride < 2; stride++) {
    ysize[stride] = (int32_T)sz1[stride];
  }

  b_emxInit_real_T(&b_coordList2, 3, TRUE);
  i0 = b_coordList2->size[0] * b_coordList2->size[1] * b_coordList2->size[2];
  b_coordList2->size[0] = ysize[0];
  b_coordList2->size[1] = ysize[1];
  b_coordList2->size[2] = 2;
  emxEnsureCapacity((emxArray__common *)b_coordList2, i0, (int32_T)sizeof(real_T));
  iy = -1;
  i0 = center2->size[0] * center2->size[1];
  for (stride = 1; stride <= i0; stride++) {
    iy++;
    b_coordList2->data[iy] = center2->data[stride - 1];
  }

  emxFree_real_T(&center2);
  i0 = coordList2->size[0] * coordList2->size[1];
  for (stride = 1; stride <= i0; stride++) {
    iy++;
    b_coordList2->data[iy] = coordList2->data[stride - 1];
  }

  emxFree_real_T(&coordList2);
  b_emxInit_real_T(&c_coordList2, 3, TRUE);

  /*  calculate gaussList */
  /* Jonas was missing the minus sign in erfc. I corrected that - KJ */
  i0 = c_coordList2->size[0] * c_coordList2->size[1] * c_coordList2->size[2];
  c_coordList2->size[0] = b_coordList2->size[0];
  c_coordList2->size[1] = b_coordList2->size[1];
  c_coordList2->size[2] = 2;
  emxEnsureCapacity((emxArray__common *)c_coordList2, i0, (int32_T)sizeof(real_T));
  loop_ub = b_coordList2->size[0] * b_coordList2->size[1] * b_coordList2->size[2];
  for (i0 = 0; i0 < loop_ub; i0++) {
    c_coordList2->data[i0] = -b_coordList2->data[i0] / 1.4142135623730951;
  }

  b_erfc(c_coordList2, b_coordList2);
  i0 = b_coordList2->size[0] * b_coordList2->size[1] * b_coordList2->size[2];
  b_coordList2->size[0] = b_coordList2->size[0];
  b_coordList2->size[1] = b_coordList2->size[1];
  b_coordList2->size[2] = 2;
  emxEnsureCapacity((emxArray__common *)b_coordList2, i0, (int32_T)sizeof(real_T));
  ixstart = b_coordList2->size[0];
  stride = b_coordList2->size[1];
  ix = b_coordList2->size[2];
  loop_ub = ixstart * stride * ix;
  emxFree_real_T(&c_coordList2);
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_coordList2->data[i0] *= 0.5;
  }

  for (i0 = 0; i0 < 3; i0++) {
    ysize[i0] = b_coordList2->size[i0];
  }

  b_emxInit_real_T(&b_y1, 3, TRUE);
  i0 = b_y1->size[0] * b_y1->size[1] * b_y1->size[2];
  b_y1->size[0] = ysize[0];
  b_y1->size[1] = ysize[1];
  b_y1->size[2] = 1;
  emxEnsureCapacity((emxArray__common *)b_y1, i0, (int32_T)sizeof(real_T));
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
  for (i0 = 0; i0 < 3; i0++) {
    ysize[i0] = b_y1->size[i0];
  }

  i0 = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
  gaussList->size[0] = ysize[0];
  gaussList->size[1] = 1;
  gaussList->size[2] = 1;
  emxEnsureCapacity((emxArray__common *)gaussList, i0, (int32_T)sizeof(real_T));
  if ((b_y1->size[0] == 0) || (b_y1->size[1] == 0)) {
    i0 = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
    gaussList->size[0] = ysize[0];
    gaussList->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)gaussList, i0, (int32_T)sizeof(real_T));
    i0 = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
    gaussList->size[2] = 1;
    emxEnsureCapacity((emxArray__common *)gaussList, i0, (int32_T)sizeof(real_T));
    loop_ub = ysize[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      gaussList->data[i0] = 1.0;
    }
  } else {
    iy = -1;
    ixstart = 0;
    for (stride = 1; stride <= b_y1->size[0]; stride++) {
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
  emxInit_int32_T(&r1, 1, TRUE);

  /*  norm gaussList */
  loop_ub = sigma2->size[1];
  i0 = r1->size[0];
  r1->size[0] = loop_ub;
  emxEnsureCapacity((emxArray__common *)r1, i0, (int32_T)sizeof(int32_T));
  for (i0 = 0; i0 < loop_ub; i0++) {
    r1->data[i0] = 1 + i0;
  }

  if (r1->size[0] == 0) {
    p = 1.0;
  } else {
    p = sigma2->data[0];
    k = 2;
    emxInit_int32_T(&r2, 1, TRUE);
    do {
      exitg1 = 0;
      loop_ub = sigma2->size[1];
      i0 = r2->size[0];
      r2->size[0] = loop_ub;
      emxEnsureCapacity((emxArray__common *)r2, i0, (int32_T)sizeof(int32_T));
      for (i0 = 0; i0 < loop_ub; i0++) {
        r2->data[i0] = 1 + i0;
      }

      if (k <= r2->size[0]) {
        p *= sigma2->data[sigma2->size[0] * (k - 1)];
        k++;
      } else {
        exitg1 = 1;
      }
    } while (exitg1 == 0);

    emxFree_int32_T(&r2);
  }

  emxFree_int32_T(&r1);
  emxFree_real_T(&sigma2);
  p *= muDoubleScalarPower(6.2831853071795862, 0.5 * (real_T)coordList->size[1]);
  i0 = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
  gaussList->size[0] = gaussList->size[0];
  gaussList->size[1] = 1;
  gaussList->size[2] = gaussList->size[2];
  emxEnsureCapacity((emxArray__common *)gaussList, i0, (int32_T)sizeof(real_T));
  ixstart = gaussList->size[0];
  stride = gaussList->size[1];
  ix = gaussList->size[2];
  loop_ub = ixstart * stride * ix;
  for (i0 = 0; i0 < loop_ub; i0++) {
    gaussList->data[i0] *= p;
  }

  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (GaussListND_mexCode.c) */
