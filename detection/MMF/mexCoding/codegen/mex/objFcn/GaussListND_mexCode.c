/*
 * GaussListND_mexCode.c
 *
 * Code generation for function 'GaussListND_mexCode'
 *
 * C source code generated on: Fri May 25 21:48:51 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "objFcn.h"
#include "GaussListND_mexCode.h"
#include "objFcn_emxutil.h"
#include "erfc.h"
#include "rdivide.h"
#include "repmat.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

/*
 * function gaussList = GaussListND_mexCode(coordList,sigma,center,intNorm)
 */
void GaussListND_mexCode(const emxArray_real_T *coordList, real_T sigma, real_T center, emxArray_real_T *gaussList)
{
    emxArray_real_T *center2;
    int32_T vlen;
    real_T nCoords[2];
    emxArray_real_T *b_coordList;
    emxArray_real_T *b_center2;
    int32_T iyLead;
    int32_T vstride;
    emxArray_real_T *coordList2;
    emxArray_real_T *varargin_1;
    real_T tmp2;
    real_T tmp1;
    int32_T ysize[3];
    int32_T sz1[2];
    int32_T dimSize;
    emxArray_real_T *b_coordList2;
    int32_T iy;
    emxArray_real_T *c_coordList2;
    int32_T ixLead;
    emxArray_real_T *b_y1;
    int32_T stride;
    int32_T k;
    int32_T ix;
    real_T work;
    emlrtHeapReferenceStackEnterFcn();
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
    /* 'GaussListND_mexCode:34' nIn = nargin; */
    /*  the following doesn't work with Matlab 6.5.0 */
    /*  error(nargchk(2,4,nIn,'struct')); */
    /* 'GaussListND_mexCode:37' if nIn < 2 || nIn > 5 */
    /*  % check dimensionality of coordList. */
    /*  if isempty(coordList) */
    /*      error('you have to supply a list of coordinates for GaussList23D') */
    /*  else */
    /*      [nCoords,nDims] = size(coordList); */
    /*  end */
    /* 'GaussListND_mexCode:48' [nCoords,nDims] = size(coordList); */
    vlen = coordList->size[0];
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
    /* 'GaussListND_mexCode:77' center2 = repmat(center, [nCoords,1]); */
    nCoords[0] = (real_T)vlen;
    nCoords[1] = 1.0;
    repmat(center, nCoords, center2);
    /* clear center */
    /*  intNorm */
    /* 'GaussListND_mexCode:81' if nIn < 4 || isempty(intNorm) */
    /* 'GaussListND_mexCode:82' intNorm = 0; */
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
    /* 'GaussListND_mexCode:163' center3 = center2(:,1); */
    /* clear center2 */
    /*  convert coordList to 0/1 */
    /* 'GaussListND_mexCode:167' coordList2 = (coordList(1:nCoords) - center3(1:nCoords))./sigma; */
    if (1 > vlen) {
        vlen = 0;
    }
    emxInit_real_T(&b_coordList, 1, TRUE);
    emxInit_real_T(&b_center2, 1, TRUE);
    iyLead = b_center2->size[0];
    b_center2->size[0] = center2->size[0];
    emxEnsureCapacity((emxArray__common *)b_center2, iyLead, (int32_T)sizeof(real_T));
    vstride = center2->size[0] - 1;
    for (iyLead = 0; iyLead <= vstride; iyLead++) {
        b_center2->data[iyLead] = center2->data[iyLead];
    }
    emxFree_real_T(&center2);
    iyLead = b_coordList->size[0];
    b_coordList->size[0] = vlen;
    emxEnsureCapacity((emxArray__common *)b_coordList, iyLead, (int32_T)sizeof(real_T));
    vstride = vlen - 1;
    for (iyLead = 0; iyLead <= vstride; iyLead++) {
        b_coordList->data[iyLead] = coordList->data[iyLead] - b_center2->data[iyLead];
    }
    emxFree_real_T(&b_center2);
    emxInit_real_T(&coordList2, 1, TRUE);
    emxInit_real_T(&varargin_1, 1, TRUE);
    rdivide(b_coordList, sigma, coordList2);
    /* clear coordList center3 */
    /*  double coordList as preparation for erfc */
    /* fixed bug: must divide the 0.5 by sigma - KJ */
    /* 'GaussListND_mexCode:172' coordList2 = cat(3,coordList2-0.5./sigma, coordList2+0.5./sigma); */
    tmp2 = 0.5 / sigma;
    tmp1 = 0.5 / sigma;
    iyLead = varargin_1->size[0];
    varargin_1->size[0] = coordList2->size[0];
    emxEnsureCapacity((emxArray__common *)varargin_1, iyLead, (int32_T)sizeof(real_T));
    emxFree_real_T(&b_coordList);
    vstride = coordList2->size[0] - 1;
    for (iyLead = 0; iyLead <= vstride; iyLead++) {
        varargin_1->data[iyLead] = coordList2->data[iyLead] - tmp2;
    }
    iyLead = coordList2->size[0];
    emxEnsureCapacity((emxArray__common *)coordList2, iyLead, (int32_T)sizeof(real_T));
    vstride = coordList2->size[0] - 1;
    for (iyLead = 0; iyLead <= vstride; iyLead++) {
        coordList2->data[iyLead] += tmp1;
    }
    for (iyLead = 0; iyLead < 3; iyLead++) {
        ysize[iyLead] = 1;
    }
    sz1[0] = varargin_1->size[0];
    sz1[1] = 1;
    for (dimSize = 0; dimSize < 2; dimSize++) {
        ysize[dimSize] = sz1[dimSize];
    }
    c_emxInit_real_T(&b_coordList2, 3, TRUE);
    ysize[2] = 2;
    iyLead = b_coordList2->size[0] * b_coordList2->size[1] * b_coordList2->size[2];
    b_coordList2->size[0] = ysize[0];
    b_coordList2->size[1] = ysize[1];
    b_coordList2->size[2] = 2;
    emxEnsureCapacity((emxArray__common *)b_coordList2, iyLead, (int32_T)sizeof(real_T));
    iy = -1;
    vstride = varargin_1->size[0];
    for (dimSize = 1; dimSize <= vstride; dimSize++) {
        iy++;
        b_coordList2->data[iy] = varargin_1->data[dimSize - 1];
    }
    emxFree_real_T(&varargin_1);
    vstride = coordList2->size[0];
    for (dimSize = 1; dimSize <= vstride; dimSize++) {
        iy++;
        b_coordList2->data[iy] = coordList2->data[dimSize - 1];
    }
    emxFree_real_T(&coordList2);
    c_emxInit_real_T(&c_coordList2, 3, TRUE);
    /*  calculate gaussList */
    /* Jonas was missing the minus sign in erfc. I corrected that - KJ */
    /* 'GaussListND_mexCode:176' gaussList = diff(0.5 * erfc(-coordList2/sqrt(2)),1,3); */
    iyLead = c_coordList2->size[0] * c_coordList2->size[1] * c_coordList2->size[2];
    c_coordList2->size[0] = b_coordList2->size[0];
    c_coordList2->size[1] = b_coordList2->size[1];
    c_coordList2->size[2] = b_coordList2->size[2];
    emxEnsureCapacity((emxArray__common *)c_coordList2, iyLead, (int32_T)sizeof(real_T));
    vstride = b_coordList2->size[0] * b_coordList2->size[1] * b_coordList2->size[2] - 1;
    for (iyLead = 0; iyLead <= vstride; iyLead++) {
        c_coordList2->data[iyLead] = -b_coordList2->data[iyLead] / 1.4142135623730951;
    }
    b_erfc(c_coordList2, b_coordList2);
    iyLead = b_coordList2->size[0] * b_coordList2->size[1] * b_coordList2->size[2];
    emxEnsureCapacity((emxArray__common *)b_coordList2, iyLead, (int32_T)sizeof(real_T));
    vlen = b_coordList2->size[0];
    ixLead = b_coordList2->size[1];
    iyLead = b_coordList2->size[2];
    emxFree_real_T(&c_coordList2);
    vstride = vlen * ixLead * iyLead - 1;
    for (iyLead = 0; iyLead <= vstride; iyLead++) {
        b_coordList2->data[iyLead] *= 0.5;
    }
    dimSize = b_coordList2->size[2];
    c_emxInit_real_T(&b_y1, 3, TRUE);
    if (dimSize <= 1) {
        for (iyLead = 0; iyLead < 3; iyLead++) {
            ysize[iyLead] = b_coordList2->size[iyLead];
        }
        ysize[2] = 0;
        iyLead = b_y1->size[0] * b_y1->size[1] * b_y1->size[2];
        b_y1->size[0] = ysize[0];
        emxEnsureCapacity((emxArray__common *)b_y1, iyLead, (int32_T)sizeof(real_T));
        iyLead = b_y1->size[0] * b_y1->size[1] * b_y1->size[2];
        b_y1->size[1] = ysize[1];
        b_y1->size[2] = 0;
        emxEnsureCapacity((emxArray__common *)b_y1, iyLead, (int32_T)sizeof(real_T));
    } else {
        for (iyLead = 0; iyLead < 3; iyLead++) {
            ysize[iyLead] = b_coordList2->size[iyLead];
        }
        ysize[2] = dimSize - 1;
        iyLead = b_y1->size[0] * b_y1->size[1] * b_y1->size[2];
        b_y1->size[0] = ysize[0];
        b_y1->size[1] = ysize[1];
        b_y1->size[2] = ysize[2];
        emxEnsureCapacity((emxArray__common *)b_y1, iyLead, (int32_T)sizeof(real_T));
        stride = 1;
        for (k = 0; k < 2; k++) {
            vlen = b_coordList2->size[k];
            stride *= vlen;
        }
        ix = 0;
        iy = 0;
        for (vlen = 1; vlen <= stride; vlen++) {
            ixLead = ix + stride;
            iyLead = iy;
            work = b_coordList2->data[ix];
            for (vstride = 2; vstride <= dimSize; vstride++) {
                tmp2 = work;
                work = b_coordList2->data[ixLead];
                tmp1 = b_coordList2->data[ixLead] - tmp2;
                ixLead += stride;
                b_y1->data[iyLead] = tmp1;
                iyLead += stride;
            }
            ix++;
            iy++;
        }
    }
    emxFree_real_T(&b_coordList2);
    /* 'GaussListND_mexCode:177' gaussList = prod(gaussList,2); */
    for (iyLead = 0; iyLead < 3; iyLead++) {
        ysize[iyLead] = b_y1->size[iyLead];
    }
    ysize[1] = 1;
    iyLead = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
    gaussList->size[0] = ysize[0];
    gaussList->size[1] = 1;
    gaussList->size[2] = ysize[2];
    emxEnsureCapacity((emxArray__common *)gaussList, iyLead, (int32_T)sizeof(real_T));
    if ((b_y1->size[0] == 0) || (b_y1->size[1] == 0) || (b_y1->size[2] == 0)) {
        iyLead = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
        gaussList->size[1] = 1;
        emxEnsureCapacity((emxArray__common *)gaussList, iyLead, (int32_T)sizeof(real_T));
        vstride = gaussList->size[2] - 1;
        for (iyLead = 0; iyLead <= vstride; iyLead++) {
            vlen = gaussList->size[0] - 1;
            for (ixLead = 0; ixLead <= vlen; ixLead++) {
                gaussList->data[ixLead + gaussList->size[0] * gaussList->size[1] * iyLead] = 1.0;
            }
        }
    } else {
        vlen = b_y1->size[1];
        vstride = b_y1->size[0];
        ixLead = b_y1->size[2];
        ix = 0;
        iy = -1;
        for (iyLead = 1; iyLead <= ixLead; iyLead++) {
            stride = ix;
            for (dimSize = 1; dimSize <= vstride; dimSize++) {
                stride++;
                ix = stride;
                tmp2 = b_y1->data[stride - 1];
                for (k = 2; k <= vlen; k++) {
                    ix += vstride;
                    tmp2 *= b_y1->data[ix - 1];
                }
                iy++;
                gaussList->data[iy] = tmp2;
            }
        }
    }
    emxFree_real_T(&b_y1);
    /*  norm gaussList */
    /* 'GaussListND_mexCode:180' switch intNorm */
    /* 'GaussListND_mexCode:181' case 0 */
    /*  "un-norm" Gaussian */
    /* 'GaussListND_mexCode:183' gaussList = gaussList*((2*pi)^(0.5*nDims)*prod(sigma(1,:))); */
    tmp2 = 2.5066282746310002 * sigma;
    iyLead = gaussList->size[0] * gaussList->size[1] * gaussList->size[2];
    gaussList->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)gaussList, iyLead, (int32_T)sizeof(real_T));
    vlen = gaussList->size[0];
    ixLead = gaussList->size[1];
    iyLead = gaussList->size[2];
    vstride = vlen * ixLead * iyLead - 1;
    for (iyLead = 0; iyLead <= vstride; iyLead++) {
        gaussList->data[iyLead] *= tmp2;
    }
    emlrtHeapReferenceStackLeaveFcn();
}
/* End of code generation (GaussListND_mexCode.c) */
