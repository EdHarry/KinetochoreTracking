/*
 * fitNGaussians3D_mexCode_F.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode_F'
 *
 * C source code generated on: Mon May 28 14:50:33 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_F.h"
#include "squeeze.h"
#include "GaussListND_mexCode.h"
#include "fitNGaussians3D_mexCode_F_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static real_T rt_roundd_snf(real_T u);

/* Function Definitions */
static real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

/*
 * function F = fitNGaussians3D_mexCode_F(x0,image,index,psfSigma)
 *  Edit of fitNGaussians2D to work in 3D
 *    EHarry March 2012
 */
void fitNGaussians3D_mexCode_F(emxArray_real_T *x0, const emxArray_real_T *image,
  const emxArray_real_T *b_index, const real_T psfSigma[2], emxArray_real_T *F)
{
  real_T bgAmp;
  int32_T i0;
  emxArray_real_T *b_x0;
  int32_T ix;
  int32_T loop_ub;
  real_T numPSF;
  int32_T nx;
  int32_T sz[2];
  emxArray_real_T *c_x0;
  int32_T ixstart;
  emxArray_real_T *d_x0;
  emxArray_int32_T *r0;
  int32_T n;
  real_T mtmp;
  boolean_T exitg7;
  emxArray_int32_T *r1;
  real_T b_mtmp;
  boolean_T exitg6;
  emxArray_int32_T *r2;
  real_T c_mtmp;
  boolean_T exitg5;
  emxArray_int32_T *r3;
  real_T d_mtmp;
  boolean_T exitg4;
  emxArray_int32_T *r4;
  real_T e_mtmp;
  boolean_T exitg3;
  emxArray_int32_T *r5;
  real_T f_mtmp;
  boolean_T exitg2;
  emxArray_real_T *psfIntegX;
  int32_T i;
  emxArray_real_T *temp;
  emxArray_real_T *temp2;
  emxArray_real_T *y;
  emxArray_real_T *b_y;
  real_T anew;
  real_T apnd;
  real_T ndbl;
  real_T cdiff;
  real_T absa;
  real_T absb;
  int32_T npages;
  emxArray_real_T *psfIntegY;
  emxArray_real_T *c_y;
  emxArray_real_T *psfIntegZ;
  emxArray_real_T *d_y;
  emxArray_int32_T *r6;
  int32_T iv0[2];
  int32_T outsize[2];
  emxArray_real_T *x;
  emxArray_int32_T *r7;
  int32_T exitg1;
  emxArray_real_T *r8;
  emxArray_real_T *r9;
  emxArray_real_T *r10;
  uint32_T b_sz[2];
  emxArray_real_T *r11;
  emxArray_boolean_T *indxPixel;
  emxArray_int32_T *r12;
  emxArray_real_T *b_F;
  emxArray_real_T *c_F;

  /* % ORIGINAL HEADER */
  /*  % %FITNGAUSSIANS2D yields F, the difference between an image and a theoretical image produced by N Gaussians, and J, the Jacobian of F. */
  /*  % % */
  /*  % %SYNOPSIS [F,J] = fitNGaussians2D(x0,image,index,psfSigma) */
  /*  % % */
  /*  % %INPUT  x0      : initial guess of PSF positions and amplitudes and */
  /*  % %                 background noise. */
  /*  % %       image   : Image part being analyzed. */
  /*  % %       index   : x,y-indices of pixels considered. */
  /*  % %       psfSigma: Standard deviation of point spread function (in pixels). */
  /*  % */
  /*  % %OUTPUT F       : Residuals from fitting an image with supplied */
  /*  % %                 Gaussians. */
  /*  % %       J       : The Jacobian matrix of F. */
  /*  % %       errFlag : 0 if function executes normally, 1 otherwise. */
  /*  % % */
  /*  % %REMARKS F = model image - real image, important to know if the sign of the */
  /*  % %residuals matters. */
  /*  % % */
  /*  % %Khuloud Jaqaman, August 2005 */
  /* % Output */
  /* 'fitNGaussians3D_mexCode_F:28' F = []; */
  /*  J = []; */
  /* % Input */
  /*  %check whether correct number of input arguments was used */
  /*  if nargin ~= 4 */
  /*      disp('--fitNGaussians2D: Incorrect number of input arguments!'); */
  /*      return */
  /*  end */
  /* check whether correct number of input arguments was used */
  /* 'fitNGaussians3D_mexCode_F:39' if nargin ~= 4 */
  /* % Calculating F & J */
  /* extract background intensity from x0 and remove from vector */
  /* 'fitNGaussians3D_mexCode_F:47' bgAmp = x0(end); */
  bgAmp = x0->data[x0->size[0] - 1];

  /* 'fitNGaussians3D_mexCode_F:48' x0 = x0(1:end-1); */
  if (1 > (x0->size[0] - 1)) {
    i0 = -1;
  } else {
    i0 = x0->size[0] - 2;
  }

  c_emxInit_real_T(&b_x0, 1);
  ix = b_x0->size[0];
  b_x0->size[0] = i0 + 1;
  emxEnsureCapacity((emxArray__common *)b_x0, ix, (int32_T)(sizeof(real_T)));
  for (ix = 0; ix <= i0; ix++) {
    b_x0->data[ix] = x0->data[ix];
  }

  i0 = x0->size[0];
  x0->size[0] = b_x0->size[0];
  emxEnsureCapacity((emxArray__common *)x0, i0, (int32_T)(sizeof(real_T)));
  loop_ub = b_x0->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    x0->data[i0] = b_x0->data[i0];
  }

  emxFree_real_T(&b_x0);

  /* get number of PSFs considered */
  /*  numPSF = length(x0)/3; */
  /* 'fitNGaussians3D_mexCode_F:52' numPSF = length(x0)/4; */
  numPSF = ((real_T)x0->size[0]) / 4.0;

  /*  %reshape 3nx1 vector x0 into nx3 matrix */
  /*  x0 = reshape(x0,3,numPSF); */
  /*  x0 = x0'; */
  /* reshape 4nx1 vector x0 into nx4 matrix */
  /* 'fitNGaussians3D_mexCode_F:59' x0 = reshape(x0,4,numPSF); */
  nx = x0->size[0];
  for (i0 = 0; i0 < 2; i0++) {
    sz[i0] = 0;
  }

  emxInit_real_T(&c_x0, 2);
  sz[0] = 4;
  i0 = (int32_T)rt_roundd_snf(numPSF);
  sz[1] = i0;
  i0 = c_x0->size[0] * c_x0->size[1];
  c_x0->size[0] = 4;
  c_x0->size[1] = sz[1];
  emxEnsureCapacity((emxArray__common *)c_x0, i0, (int32_T)(sizeof(real_T)));
  for (ixstart = 0; (ixstart + 1) <= nx; ixstart++) {
    c_x0->data[ixstart] = x0->data[ixstart];
  }

  emxInit_real_T(&d_x0, 2);

  /* 'fitNGaussians3D_mexCode_F:60' x0 = x0'; */
  i0 = d_x0->size[0] * d_x0->size[1];
  d_x0->size[0] = c_x0->size[1];
  d_x0->size[1] = 4;
  emxEnsureCapacity((emxArray__common *)d_x0, i0, (int32_T)(sizeof(real_T)));
  for (i0 = 0; i0 < 4; i0++) {
    loop_ub = c_x0->size[1] - 1;
    for (ix = 0; ix <= loop_ub; ix++) {
      d_x0->data[ix + (d_x0->size[0] * i0)] = c_x0->data[i0 + (c_x0->size[0] *
        ix)];
    }
  }

  emxFree_real_T(&c_x0);
  emxInit_int32_T(&r0, 1);

  /* extract PSF center positions and amplitudes */
  /*  psfPos = x0(:,1:2); */
  /*  psfAmp = x0(:,3); */
  /* 'fitNGaussians3D_mexCode_F:66' psfPos = x0(:,1:3); */
  /* 'fitNGaussians3D_mexCode_F:67' psfAmp = x0(:,4); */
  /* find minimum and maximum pixel indices */
  /* 'fitNGaussians3D_mexCode_F:70' minIndxX = min(index(:,1)); */
  ixstart = 1;
  i0 = b_index->size[0];
  ix = r0->size[0];
  r0->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r0, ix, (int32_T)(sizeof(int32_T)));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r0->data[i0] = 1 + i0;
  }

  n = r0->size[0];
  mtmp = b_index->data[0];
  emxFree_int32_T(&r0);
  if (n > 1) {
    if (rtIsNaN(b_index->data[0])) {
      ix = 2;
      exitg7 = FALSE;
      while ((exitg7 == 0U) && (ix <= n)) {
        ixstart = ix;
        if (!rtIsNaN(b_index->data[ix - 1])) {
          mtmp = b_index->data[ix - 1];
          exitg7 = TRUE;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < n) {
      while ((ixstart + 1) <= n) {
        if (b_index->data[ixstart] < mtmp) {
          mtmp = b_index->data[ixstart];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r1, 1);

  /* 'fitNGaussians3D_mexCode_F:71' maxIndxX = max(index(:,1)); */
  ixstart = 1;
  i0 = b_index->size[0];
  ix = r1->size[0];
  r1->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r1, ix, (int32_T)(sizeof(int32_T)));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r1->data[i0] = 1 + i0;
  }

  n = r1->size[0];
  b_mtmp = b_index->data[0];
  emxFree_int32_T(&r1);
  if (n > 1) {
    if (rtIsNaN(b_index->data[0])) {
      ix = 2;
      exitg6 = FALSE;
      while ((exitg6 == 0U) && (ix <= n)) {
        ixstart = ix;
        if (!rtIsNaN(b_index->data[ix - 1])) {
          b_mtmp = b_index->data[ix - 1];
          exitg6 = TRUE;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < n) {
      while ((ixstart + 1) <= n) {
        if (b_index->data[ixstart] > b_mtmp) {
          b_mtmp = b_index->data[ixstart];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r2, 1);

  /* 'fitNGaussians3D_mexCode_F:72' minIndxY = min(index(:,2)); */
  ixstart = 1;
  i0 = b_index->size[0];
  ix = r2->size[0];
  r2->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r2, ix, (int32_T)(sizeof(int32_T)));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r2->data[i0] = 1 + i0;
  }

  n = r2->size[0];
  c_mtmp = b_index->data[b_index->size[0]];
  emxFree_int32_T(&r2);
  if (n > 1) {
    if (rtIsNaN(b_index->data[b_index->size[0]])) {
      ix = 2;
      exitg5 = FALSE;
      while ((exitg5 == 0U) && (ix <= n)) {
        ixstart = ix;
        if (!rtIsNaN(b_index->data[(ix + b_index->size[0]) - 1])) {
          c_mtmp = b_index->data[(ix + b_index->size[0]) - 1];
          exitg5 = TRUE;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < n) {
      while ((ixstart + 1) <= n) {
        if (b_index->data[ixstart + b_index->size[0]] < c_mtmp) {
          c_mtmp = b_index->data[ixstart + b_index->size[0]];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r3, 1);

  /* 'fitNGaussians3D_mexCode_F:73' maxIndxY = max(index(:,2)); */
  ixstart = 1;
  i0 = b_index->size[0];
  ix = r3->size[0];
  r3->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r3, ix, (int32_T)(sizeof(int32_T)));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r3->data[i0] = 1 + i0;
  }

  n = r3->size[0];
  d_mtmp = b_index->data[b_index->size[0]];
  emxFree_int32_T(&r3);
  if (n > 1) {
    if (rtIsNaN(b_index->data[b_index->size[0]])) {
      ix = 2;
      exitg4 = FALSE;
      while ((exitg4 == 0U) && (ix <= n)) {
        ixstart = ix;
        if (!rtIsNaN(b_index->data[(ix + b_index->size[0]) - 1])) {
          d_mtmp = b_index->data[(ix + b_index->size[0]) - 1];
          exitg4 = TRUE;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < n) {
      while ((ixstart + 1) <= n) {
        if (b_index->data[ixstart + b_index->size[0]] > d_mtmp) {
          d_mtmp = b_index->data[ixstart + b_index->size[0]];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r4, 1);

  /* 'fitNGaussians3D_mexCode_F:74' minIndxZ = min(index(:,3)); */
  ixstart = 1;
  i0 = b_index->size[0];
  ix = r4->size[0];
  r4->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r4, ix, (int32_T)(sizeof(int32_T)));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r4->data[i0] = 1 + i0;
  }

  n = r4->size[0];
  e_mtmp = b_index->data[b_index->size[0] << 1];
  emxFree_int32_T(&r4);
  if (n > 1) {
    if (rtIsNaN(b_index->data[b_index->size[0] << 1])) {
      ix = 2;
      exitg3 = FALSE;
      while ((exitg3 == 0U) && (ix <= n)) {
        ixstart = ix;
        if (!rtIsNaN(b_index->data[(ix + (b_index->size[0] << 1)) - 1])) {
          e_mtmp = b_index->data[(ix + (b_index->size[0] << 1)) - 1];
          exitg3 = TRUE;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < n) {
      while ((ixstart + 1) <= n) {
        if (b_index->data[ixstart + (b_index->size[0] << 1)] < e_mtmp) {
          e_mtmp = b_index->data[ixstart + (b_index->size[0] << 1)];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r5, 1);

  /* 'fitNGaussians3D_mexCode_F:75' maxIndxZ = max(index(:,3)); */
  ixstart = 1;
  i0 = b_index->size[0];
  ix = r5->size[0];
  r5->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r5, ix, (int32_T)(sizeof(int32_T)));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r5->data[i0] = 1 + i0;
  }

  n = r5->size[0];
  f_mtmp = b_index->data[b_index->size[0] << 1];
  emxFree_int32_T(&r5);
  if (n > 1) {
    if (rtIsNaN(b_index->data[b_index->size[0] << 1])) {
      ix = 2;
      exitg2 = FALSE;
      while ((exitg2 == 0U) && (ix <= n)) {
        ixstart = ix;
        if (!rtIsNaN(b_index->data[(ix + (b_index->size[0] << 1)) - 1])) {
          f_mtmp = b_index->data[(ix + (b_index->size[0] << 1)) - 1];
          exitg2 = TRUE;
        } else {
          ix++;
        }
      }
    }

    if (ixstart < n) {
      while ((ixstart + 1) <= n) {
        if (b_index->data[ixstart + (b_index->size[0] << 1)] > f_mtmp) {
          f_mtmp = b_index->data[ixstart + (b_index->size[0] << 1)];
        }

        ixstart++;
      }
    }
  }

  emxInit_real_T(&psfIntegX, 2);

  /* determine the contribution of each PSF (assuming amplitude 1) to a */
  /* pixel based on its x-coordinate (needed to calculate F & J) */
  /* 'fitNGaussians3D_mexCode_F:79' psfIntegX = zeros(maxIndxX-minIndxX+1,numPSF); */
  i0 = psfIntegX->size[0] * psfIntegX->size[1];
  psfIntegX->size[0] = (int32_T)((b_mtmp - mtmp) + 1.0);
  psfIntegX->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfIntegX, i0, (int32_T)(sizeof(real_T)));
  loop_ub = (((int32_T)((b_mtmp - mtmp) + 1.0)) * ((int32_T)numPSF)) - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfIntegX->data[i0] = 0.0;
  }

  /* 'fitNGaussians3D_mexCode_F:80' for i=1:numPSF */
  i = 0;
  b_emxInit_real_T(&temp, 3);
  emxInit_real_T(&temp2, 2);
  emxInit_real_T(&y, 2);
  c_emxInit_real_T(&b_y, 1);
  while (i <= (((int32_T)numPSF) - 1)) {
    /* 'fitNGaussians3D_mexCode_F:81' temp = GaussListND_mexCode((minIndxX:maxIndxX)',... */
    /* 'fitNGaussians3D_mexCode_F:82'         psfSigma(1),psfPos(i,1)); */
    if ((rtIsNaN(mtmp)) || (rtIsNaN(b_mtmp))) {
      n = 1;
      anew = rtNaN;
      apnd = b_mtmp;
    } else if (b_mtmp < mtmp) {
      n = 0;
      anew = mtmp;
      apnd = b_mtmp;
    } else if ((rtIsInf(mtmp)) || (rtIsInf(b_mtmp))) {
      n = 1;
      anew = rtNaN;
      apnd = b_mtmp;
    } else {
      anew = mtmp;
      ndbl = floor((b_mtmp - mtmp) + 0.5);
      apnd = mtmp + ndbl;
      cdiff = apnd - b_mtmp;
      absa = fabs(mtmp);
      absb = fabs(b_mtmp);
      if (absa > absb) {
        absb = absa;
      }

      if (fabs(cdiff) < (4.4408920985006262E-16 * absb)) {
        ndbl++;
        apnd = b_mtmp;
      } else if (cdiff > 0.0) {
        apnd = mtmp + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int32_T)ndbl;
      } else {
        n = 0;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)(sizeof(real_T)));
    if (n > 0) {
      y->data[0] = anew;
      if (n > 1) {
        y->data[n - 1] = apnd;
        nx = n - 1;
        npages = nx / 2;
        for (ixstart = 1; ixstart <= (npages - 1); ixstart++) {
          y->data[ixstart] = anew + ((real_T)ixstart);
          y->data[(n - ixstart) - 1] = apnd - ((real_T)ixstart);
        }

        if ((npages << 1) == nx) {
          y->data[npages] = (anew + apnd) / 2.0;
        } else {
          y->data[npages] = anew + ((real_T)npages);
          y->data[npages + 1] = apnd - ((real_T)npages);
        }
      }
    }

    i0 = b_y->size[0];
    b_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)b_y, i0, (int32_T)(sizeof(real_T)));
    loop_ub = y->size[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      b_y->data[i0] = y->data[i0];
    }

    GaussListND_mexCode(b_y, psfSigma[0], d_x0->data[i], temp);

    /* 'fitNGaussians3D_mexCode_F:84' temp2 = squeeze(temp); */
    squeeze(temp, temp2);

    /* clear temp */
    /* 'fitNGaussians3D_mexCode_F:87' psfIntegX(:,i) = temp2(:,1); */
    i0 = temp2->size[0] - 1;
    for (ix = 0; ix <= i0; ix++) {
      psfIntegX->data[ix + (psfIntegX->size[0] * i)] = temp2->data[ix];
    }

    /* clear temp2 */
    i++;
  }

  emxFree_real_T(&b_y);
  emxInit_real_T(&psfIntegY, 2);

  /* determine the contribution of each PSF (assuming amplitude 1) to a */
  /* pixel based on its y-coordinate (needed to calculate F & J) */
  /* 'fitNGaussians3D_mexCode_F:93' psfIntegY = zeros(maxIndxY-minIndxY+1,numPSF); */
  i0 = psfIntegY->size[0] * psfIntegY->size[1];
  psfIntegY->size[0] = (int32_T)((d_mtmp - c_mtmp) + 1.0);
  psfIntegY->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfIntegY, i0, (int32_T)(sizeof(real_T)));
  loop_ub = (((int32_T)((d_mtmp - c_mtmp) + 1.0)) * ((int32_T)numPSF)) - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfIntegY->data[i0] = 0.0;
  }

  /* 'fitNGaussians3D_mexCode_F:94' for i=1:numPSF */
  i = 0;
  c_emxInit_real_T(&c_y, 1);
  while (i <= (((int32_T)numPSF) - 1)) {
    /* 'fitNGaussians3D_mexCode_F:95' temp = GaussListND_mexCode((minIndxY:maxIndxY)',... */
    /* 'fitNGaussians3D_mexCode_F:96'         psfSigma(1),psfPos(i,2)); */
    if ((rtIsNaN(c_mtmp)) || (rtIsNaN(d_mtmp))) {
      n = 1;
      anew = rtNaN;
      apnd = d_mtmp;
    } else if (d_mtmp < c_mtmp) {
      n = 0;
      anew = c_mtmp;
      apnd = d_mtmp;
    } else if ((rtIsInf(c_mtmp)) || (rtIsInf(d_mtmp))) {
      n = 1;
      anew = rtNaN;
      apnd = d_mtmp;
    } else {
      anew = c_mtmp;
      ndbl = floor((d_mtmp - c_mtmp) + 0.5);
      apnd = c_mtmp + ndbl;
      cdiff = apnd - d_mtmp;
      absa = fabs(c_mtmp);
      absb = fabs(d_mtmp);
      if (absa > absb) {
        absb = absa;
      }

      if (fabs(cdiff) < (4.4408920985006262E-16 * absb)) {
        ndbl++;
        apnd = d_mtmp;
      } else if (cdiff > 0.0) {
        apnd = c_mtmp + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int32_T)ndbl;
      } else {
        n = 0;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)(sizeof(real_T)));
    if (n > 0) {
      y->data[0] = anew;
      if (n > 1) {
        y->data[n - 1] = apnd;
        nx = n - 1;
        npages = nx / 2;
        for (ixstart = 1; ixstart <= (npages - 1); ixstart++) {
          y->data[ixstart] = anew + ((real_T)ixstart);
          y->data[(n - ixstart) - 1] = apnd - ((real_T)ixstart);
        }

        if ((npages << 1) == nx) {
          y->data[npages] = (anew + apnd) / 2.0;
        } else {
          y->data[npages] = anew + ((real_T)npages);
          y->data[npages + 1] = apnd - ((real_T)npages);
        }
      }
    }

    i0 = c_y->size[0];
    c_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)c_y, i0, (int32_T)(sizeof(real_T)));
    loop_ub = y->size[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      c_y->data[i0] = y->data[i0];
    }

    GaussListND_mexCode(c_y, psfSigma[0], d_x0->data[i + d_x0->size[0]], temp);

    /* 'fitNGaussians3D_mexCode_F:98' temp2 = squeeze(temp); */
    squeeze(temp, temp2);

    /* clear temp */
    /* 'fitNGaussians3D_mexCode_F:101' psfIntegY(:,i) = temp2(:,1); */
    i0 = temp2->size[0] - 1;
    for (ix = 0; ix <= i0; ix++) {
      psfIntegY->data[ix + (psfIntegY->size[0] * i)] = temp2->data[ix];
    }

    /* clear temp2 */
    i++;
  }

  emxFree_real_T(&c_y);
  emxInit_real_T(&psfIntegZ, 2);

  /* determine the contribution of each PSF (assuming amplitude 1) to a */
  /* pixel based on its z-coordinate (needed to calculate F & J) */
  /* 'fitNGaussians3D_mexCode_F:107' psfIntegZ = zeros(maxIndxZ-minIndxZ+1,numPSF); */
  i0 = psfIntegZ->size[0] * psfIntegZ->size[1];
  psfIntegZ->size[0] = (int32_T)((f_mtmp - e_mtmp) + 1.0);
  psfIntegZ->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfIntegZ, i0, (int32_T)(sizeof(real_T)));
  loop_ub = (((int32_T)((f_mtmp - e_mtmp) + 1.0)) * ((int32_T)numPSF)) - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfIntegZ->data[i0] = 0.0;
  }

  /* 'fitNGaussians3D_mexCode_F:108' for i=1:numPSF */
  i = 0;
  c_emxInit_real_T(&d_y, 1);
  while (i <= (((int32_T)numPSF) - 1)) {
    /* 'fitNGaussians3D_mexCode_F:109' temp = GaussListND_mexCode((minIndxZ:maxIndxZ)',... */
    /* 'fitNGaussians3D_mexCode_F:110'         psfSigma(2),psfPos(i,3)); */
    if ((rtIsNaN(e_mtmp)) || (rtIsNaN(f_mtmp))) {
      n = 1;
      anew = rtNaN;
      apnd = f_mtmp;
    } else if (f_mtmp < e_mtmp) {
      n = 0;
      anew = e_mtmp;
      apnd = f_mtmp;
    } else if ((rtIsInf(e_mtmp)) || (rtIsInf(f_mtmp))) {
      n = 1;
      anew = rtNaN;
      apnd = f_mtmp;
    } else {
      anew = e_mtmp;
      ndbl = floor((f_mtmp - e_mtmp) + 0.5);
      apnd = e_mtmp + ndbl;
      cdiff = apnd - f_mtmp;
      absa = fabs(e_mtmp);
      absb = fabs(f_mtmp);
      if (absa > absb) {
        absb = absa;
      }

      if (fabs(cdiff) < (4.4408920985006262E-16 * absb)) {
        ndbl++;
        apnd = f_mtmp;
      } else if (cdiff > 0.0) {
        apnd = e_mtmp + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int32_T)ndbl;
      } else {
        n = 0;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)(sizeof(real_T)));
    if (n > 0) {
      y->data[0] = anew;
      if (n > 1) {
        y->data[n - 1] = apnd;
        nx = n - 1;
        npages = nx / 2;
        for (ixstart = 1; ixstart <= (npages - 1); ixstart++) {
          y->data[ixstart] = anew + ((real_T)ixstart);
          y->data[(n - ixstart) - 1] = apnd - ((real_T)ixstart);
        }

        if ((npages << 1) == nx) {
          y->data[npages] = (anew + apnd) / 2.0;
        } else {
          y->data[npages] = anew + ((real_T)npages);
          y->data[npages + 1] = apnd - ((real_T)npages);
        }
      }
    }

    i0 = d_y->size[0];
    d_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)d_y, i0, (int32_T)(sizeof(real_T)));
    loop_ub = y->size[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      d_y->data[i0] = y->data[i0];
    }

    GaussListND_mexCode(d_y, psfSigma[1], d_x0->data[i + (d_x0->size[0] << 1)],
                        temp);

    /* 'fitNGaussians3D_mexCode_F:113' temp2 = squeeze(temp); */
    squeeze(temp, temp2);

    /* clear temp */
    /* 'fitNGaussians3D_mexCode_F:116' psfIntegZ(:,i) = temp2(:,1); */
    i0 = temp2->size[0] - 1;
    for (ix = 0; ix <= i0; ix++) {
      psfIntegZ->data[ix + (psfIntegZ->size[0] * i)] = temp2->data[ix];
    }

    /* clear temp2 */
    i++;
  }

  emxFree_real_T(&d_y);
  emxFree_real_T(&temp2);
  emxFree_real_T(&temp);
  emxInit_int32_T(&r6, 1);

  /*  %calculate the value of each PSF (assuming amplitude 1) at the */
  /*  %x-coordinates of the corners of all pixels (needed to calculate J) */
  /*  psfValueX = zeros(maxIndxX-minIndxX+2,numPSF); */
  /*  for i=1:numPSF */
  /*      psfValueX(:,i) = exp(-((minIndxX-0.5:maxIndxX+0.5)'... */
  /*          -psfPos(i,1)).^2/2/psfSigma(1)^2); */
  /*  end */
  /*   */
  /*  %calculate the value of each PSF (assuming amplitude 1) at the */
  /*  %y-coordinates of the corners of all pixels (needed to calculate J) */
  /*  psfValueY = zeros(maxIndxY-minIndxY+2,numPSF); */
  /*  for i=1:numPSF */
  /*      psfValueY(:,i) = exp(-((minIndxY-0.5:maxIndxY+0.5)'... */
  /*          -psfPos(i,2)).^2/2/psfSigma(1)^2); */
  /*  end */
  /*   */
  /*  %calculate the value of each PSF (assuming amplitude 1) at the */
  /*  %z-coordinates of the corners of all pixels (needed to calculate J) */
  /*  psfValueZ = zeros(maxIndxZ-minIndxZ+2,numPSF); */
  /*  for i=1:numPSF */
  /*      psfValueZ(:,i) = exp(-((minIndxZ-0.5:maxIndxZ+0.5)'... */
  /*          -psfPos(i,3)).^2/2/psfSigma(2)^2); */
  /*  end */
  /* get number of pixels in image */
  /* 'fitNGaussians3D_mexCode_F:145' numPixel = length(image); */
  /*  %get xy-indices relative to minimum */
  /*  relIndxX = index(:,1) - minIndxX + 1; */
  /*  relIndxY = index(:,2) - minIndxY + 1; */
  /* get xy-indices relative to minimum */
  /* 'fitNGaussians3D_mexCode_F:151' relIndxX = index(:,1) - minIndxX + 1; */
  /* 'fitNGaussians3D_mexCode_F:152' relIndxY = index(:,2) - minIndxY + 1; */
  /* 'fitNGaussians3D_mexCode_F:153' relIndxZ = index(:,3) - minIndxZ + 1; */
  /*  %calculate the value of F at all pixels */
  /*  F = (sum(repmat(psfAmp,1,numPixel).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)',1))' ... */
  /*      + repmat(bgAmp,numPixel,1) - image; */
  /* calculate the value of F at all pixels */
  /* 'fitNGaussians3D_mexCode_F:160' F = (sum(repmat(psfAmp,1,numPixel).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)'.*psfIntegZ(relIndxZ,:)',1))' ... */
  /* 'fitNGaussians3D_mexCode_F:161'     + repmat(bgAmp,numPixel,1) - image; */
  sz[0] = 1;
  sz[1] = image->size[0];
  i0 = d_x0->size[0];
  ix = r6->size[0];
  r6->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r6, ix, (int32_T)(sizeof(int32_T)));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r6->data[i0] = 1 + i0;
  }

  iv0[0] = r6->size[0];
  iv0[1] = 1;
  emxFree_int32_T(&r6);
  for (i0 = 0; i0 < 2; i0++) {
    outsize[i0] = iv0[i0] * sz[i0];
  }

  emxInit_real_T(&x, 2);
  i0 = x->size[0] * x->size[1];
  x->size[0] = outsize[0];
  x->size[1] = outsize[1];
  emxEnsureCapacity((emxArray__common *)x, i0, (int32_T)(sizeof(real_T)));
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
  } else {
    nx = 0;
    npages = 1;
    emxInit_int32_T(&r7, 1);
    while (npages <= sz[1]) {
      n = 0;
      ixstart = 1;
      do {
        exitg1 = 0;
        i0 = d_x0->size[0];
        ix = r7->size[0];
        r7->size[0] = i0;
        emxEnsureCapacity((emxArray__common *)r7, ix, (int32_T)(sizeof(int32_T)));
        loop_ub = i0 - 1;
        for (i0 = 0; i0 <= loop_ub; i0++) {
          r7->data[i0] = 1 + i0;
        }

        if (ixstart <= r7->size[0]) {
          x->data[nx] = d_x0->data[n + (d_x0->size[0] * 3)];
          n++;
          nx++;
          ixstart++;
        } else {
          exitg1 = 1;
        }
      } while (exitg1 == 0U);

      npages++;
    }

    emxFree_int32_T(&r7);
  }

  emxFree_real_T(&d_x0);
  emxInit_real_T(&r8, 2);
  i0 = psfIntegX->size[1];
  ix = b_index->size[0];
  nx = r8->size[0] * r8->size[1];
  r8->size[0] = i0;
  r8->size[1] = ix;
  emxEnsureCapacity((emxArray__common *)r8, nx, (int32_T)(sizeof(real_T)));
  loop_ub = ix - 1;
  for (ix = 0; ix <= loop_ub; ix++) {
    npages = i0 - 1;
    for (nx = 0; nx <= npages; nx++) {
      r8->data[nx + (r8->size[0] * ix)] = psfIntegX->data[(((int32_T)
        ((b_index->data[ix] - mtmp) + 1.0)) + (psfIntegX->size[0] * nx)) - 1];
    }
  }

  emxFree_real_T(&psfIntegX);
  emxInit_real_T(&r9, 2);
  i0 = psfIntegY->size[1];
  ix = b_index->size[0];
  nx = r9->size[0] * r9->size[1];
  r9->size[0] = i0;
  r9->size[1] = ix;
  emxEnsureCapacity((emxArray__common *)r9, nx, (int32_T)(sizeof(real_T)));
  loop_ub = ix - 1;
  for (ix = 0; ix <= loop_ub; ix++) {
    npages = i0 - 1;
    for (nx = 0; nx <= npages; nx++) {
      r9->data[nx + (r9->size[0] * ix)] = psfIntegY->data[(((int32_T)
        ((b_index->data[ix + b_index->size[0]] - c_mtmp) + 1.0)) +
        (psfIntegY->size[0] * nx)) - 1];
    }
  }

  emxFree_real_T(&psfIntegY);
  emxInit_real_T(&r10, 2);
  i0 = psfIntegZ->size[1];
  ix = b_index->size[0];
  nx = r10->size[0] * r10->size[1];
  r10->size[0] = i0;
  r10->size[1] = ix;
  emxEnsureCapacity((emxArray__common *)r10, nx, (int32_T)(sizeof(real_T)));
  loop_ub = ix - 1;
  for (ix = 0; ix <= loop_ub; ix++) {
    npages = i0 - 1;
    for (nx = 0; nx <= npages; nx++) {
      r10->data[nx + (r10->size[0] * ix)] = psfIntegZ->data[(((int32_T)
        ((b_index->data[ix + (b_index->size[0] << 1)] - e_mtmp) + 1.0)) +
        (psfIntegZ->size[0] * nx)) - 1];
    }
  }

  emxFree_real_T(&psfIntegZ);
  i0 = x->size[0] * x->size[1];
  x->size[0] = x->size[0];
  x->size[1] = x->size[1];
  emxEnsureCapacity((emxArray__common *)x, i0, (int32_T)(sizeof(real_T)));
  nx = x->size[0];
  npages = x->size[1];
  loop_ub = (nx * npages) - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    x->data[i0] = ((x->data[i0] * r8->data[i0]) * r9->data[i0]) * r10->data[i0];
  }

  emxFree_real_T(&r10);
  emxFree_real_T(&r9);
  emxFree_real_T(&r8);
  for (i0 = 0; i0 < 2; i0++) {
    b_sz[i0] = (uint32_T)x->size[i0];
  }

  b_sz[0] = 1U;
  i0 = y->size[0] * y->size[1];
  y->size[0] = 1;
  y->size[1] = (int32_T)b_sz[1];
  emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)(sizeof(real_T)));
  if ((x->size[0] == 0) || (x->size[1] == 0)) {
    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = y->size[1];
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)(sizeof(real_T)));
    loop_ub = y->size[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      y->data[y->size[0] * i0] = 0.0;
    }
  } else {
    nx = x->size[0];
    npages = x->size[1];
    ix = -1;
    n = -1;
    for (i = 1; i <= npages; i++) {
      ixstart = ix + 1;
      ix++;
      numPSF = x->data[ixstart];
      for (ixstart = 2; ixstart <= nx; ixstart++) {
        ix++;
        numPSF += x->data[ix];
      }

      n++;
      y->data[n] = numPSF;
    }
  }

  emxFree_real_T(&x);
  c_emxInit_real_T(&r11, 1);
  i0 = r11->size[0];
  r11->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)r11, i0, (int32_T)(sizeof(real_T)));
  loop_ub = y->size[1] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r11->data[i0] = y->data[i0];
  }

  emxFree_real_T(&y);
  nx = r11->size[0];
  i0 = F->size[0] * F->size[1];
  F->size[0] = nx;
  emxEnsureCapacity((emxArray__common *)F, i0, (int32_T)(sizeof(real_T)));
  i0 = F->size[0] * F->size[1];
  F->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)F, i0, (int32_T)(sizeof(real_T)));
  loop_ub = r11->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    F->data[i0] = (r11->data[i0] + bgAmp) - image->data[i0];
  }

  emxFree_real_T(&r11);
  emxInit_boolean_T(&indxPixel, 1);

  /* remove pixels with NaN (which means they are out of the cropped image */
  /* area) */
  /* 'fitNGaussians3D_mexCode_F:165' indxPixel = ~isnan(image); */
  i0 = indxPixel->size[0];
  indxPixel->size[0] = image->size[0];
  emxEnsureCapacity((emxArray__common *)indxPixel, i0, (int32_T)(sizeof
    (boolean_T)));
  loop_ub = image->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    indxPixel->data[i0] = rtIsNaN(image->data[i0]);
  }

  i0 = indxPixel->size[0];
  indxPixel->size[0] = indxPixel->size[0];
  emxEnsureCapacity((emxArray__common *)indxPixel, i0, (int32_T)(sizeof
    (boolean_T)));
  loop_ub = indxPixel->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    indxPixel->data[i0] = !indxPixel->data[i0];
  }

  /* 'fitNGaussians3D_mexCode_F:167' F = F(indxPixel); */
  n = indxPixel->size[0];
  ixstart = 0;
  for (i = 1; i <= n; i++) {
    if (indxPixel->data[i - 1]) {
      ixstart++;
    }
  }

  emxInit_int32_T(&r12, 1);
  i0 = r12->size[0];
  r12->size[0] = ixstart;
  emxEnsureCapacity((emxArray__common *)r12, i0, (int32_T)(sizeof(int32_T)));
  npages = 0;
  for (i = 1; i <= n; i++) {
    if (indxPixel->data[i - 1]) {
      r12->data[npages] = i;
      npages++;
    }
  }

  emxFree_boolean_T(&indxPixel);
  c_emxInit_real_T(&b_F, 1);
  i0 = b_F->size[0];
  b_F->size[0] = r12->size[0];
  emxEnsureCapacity((emxArray__common *)b_F, i0, (int32_T)(sizeof(real_T)));
  loop_ub = r12->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    b_F->data[i0] = F->data[r12->data[i0] - 1];
  }

  c_emxInit_real_T(&c_F, 1);
  i0 = c_F->size[0];
  c_F->size[0] = r12->size[0];
  emxEnsureCapacity((emxArray__common *)c_F, i0, (int32_T)(sizeof(real_T)));
  loop_ub = r12->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    c_F->data[i0] = F->data[r12->data[i0] - 1];
  }

  emxFree_int32_T(&r12);
  nx = b_F->size[0];
  i0 = F->size[0] * F->size[1];
  F->size[0] = nx;
  F->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)F, i0, (int32_T)(sizeof(real_T)));
  emxFree_real_T(&b_F);
  i0 = 0;
  while (i0 <= 0) {
    loop_ub = nx - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      F->data[i0] = c_F->data[i0];
    }

    i0 = 1;
  }

  emxFree_real_T(&c_F);

  /*  if nargout > 1 */
  /*      %calculate the derivative at all pixels */
  /*      %     J = ones(numPixel,3*numPSF+1); %(last column for background amplitude) */
  /*      %     J(:,1:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueX(relIndxX,:)-... */
  /*      %         psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:); %w.r.t. x */
  /*      %     J(:,2:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueY(relIndxY,:)-... */
  /*      %         psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:); %w.r.t. y */
  /*      %     J(:,3:3:3*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); %w.r.t. amp */
  /*       */
  /*      J = ones(numPixel,4*numPSF+1); %(last column for background amplitude) */
  /*      J(:,1:4:4*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueX(relIndxX,:)-... */
  /*          psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:).*psfIntegZ(relIndxZ,:); %w.r.t. x */
  /*      J(:,2:4:4*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueY(relIndxY,:)-... */
  /*          psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:).*psfIntegZ(relIndxZ,:); %w.r.t. y */
  /*      J(:,3:4:4*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueZ(relIndxZ,:)-... */
  /*          psfValueZ(relIndxZ+1,:)).*psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); %w.r.t. z */
  /*      J(:,4:4:4*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:).*psfIntegZ(relIndxZ,:); %w.r.t. amp */
  /*       */
  /*      %remove pixels with NaN (which means they are out of the cropped image */
  /*      %area) */
  /*      J = J(indxPixel,:); */
  /*  end */
  /* % ~~ the end ~~ */
  /* % OLD CODE */
  /*  % J = ones(numPixel,3*numPSF+1); */
  /*  % F = ones(numPixel,1); */
  /*  % */
  /*  % for i=1:numPixel %for each pixel */
  /*  % */
  /*  %     %get xy-indices relative to minimum */
  /*  %     relIndxX = index(i,1) - minIndxX + 1; */
  /*  %     relIndxY = index(i,2) - minIndxY + 1; */
  /*  % */
  /*  %     %calculate the value of F */
  /*  %     F(i) = sum(psfAmp.*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)') ... */
  /*  %         + bgAmp - image(i); */
  /*  % */
  /*  %     %calculate the derivative wrt x-coordinate */
  /*  %     J(i,1:3:3*numPSF) = psfAmp'.*(psfValueX(relIndxX,:)-... */
  /*  %         psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:)/psfSigma^2; */
  /*  % */
  /*  %     %calculate the derivative wrt y-coordinate */
  /*  %     J(i,2:3:3*numPSF) = psfAmp'.*(psfValueY(relIndxY,:)-... */
  /*  %         psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:)/psfSigma^2; */
  /*  % */
  /*  %     %calculate the derivative wrt amplitude */
  /*  %     J(i,3:3:3*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); */
  /*  % */
  /*  %     %since the derivative wrt background intensity = 1, this is already */
  /*  %     %accounted for in the initial assignment of J. */
  /*  % */
  /*  % end */
}

/* End of code generation (fitNGaussians3D_mexCode_F.c) */
