/*
 * fitNGaussians3D_mexCode.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode'
 *
 * C source code generated on: Thu May  3 12:56:20 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "squeeze.h"
#include "GaussListND_mexCode.h"
#include "fitNGaussians3D_mexCode_emxutil.h"
#include "exp.h"
#include "rdivide.h"
#include "power.h"
#include "repmat.h"
#include "sum.h"
#include "fitNGaussians3D_mexCode_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void fitNGaussians3D_mexCode(emxArray_real_T *x0, const emxArray_real_T *image,
  const emxArray_real_T *b_index, const real_T psfSigma[2], emxArray_real_T *F,
  emxArray_real_T *J)
{
  real_T bgAmp;
  int32_T i0;
  emxArray_real_T *b_x0;
  int32_T i1;
  int32_T loop_ub;
  real_T numPSF;
  int32_T nx;
  int32_T sz[2];
  emxArray_real_T *c_x0;
  int32_T k;
  emxArray_real_T *d_x0;
  emxArray_int32_T *r0;
  int32_T ixstart;
  real_T mtmp;
  int32_T nm1;
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
  emxArray_real_T *psfIntegY;
  emxArray_real_T *c_y;
  emxArray_real_T *psfIntegZ;
  emxArray_real_T *d_y;
  emxArray_real_T *psfValueX;
  emxArray_real_T *indxPixel;
  emxArray_real_T *b_indxPixel;
  emxArray_real_T *c_indxPixel;
  emxArray_real_T *e_y;
  emxArray_real_T *psfValueY;
  emxArray_real_T *d_indxPixel;
  emxArray_real_T *e_indxPixel;
  emxArray_real_T *f_y;
  emxArray_real_T *psfValueZ;
  emxArray_real_T *f_indxPixel;
  emxArray_real_T *g_indxPixel;
  emxArray_real_T *g_y;
  emxArray_real_T *relIndxX;
  emxArray_real_T *relIndxY;
  emxArray_real_T *relIndxZ;
  emxArray_real_T *e_x0;
  emxArray_real_T *r6;
  emxArray_real_T *r7;
  emxArray_real_T *r8;
  emxArray_real_T *r9;
  emxArray_real_T *r10;
  emxArray_real_T *r11;
  emxArray_boolean_T *x;
  emxArray_int32_T *ii;
  boolean_T exitg1;
  boolean_T guard1 = FALSE;
  emxArray_int32_T *b_ii;
  emxArray_real_T *b_F;
  emxArray_real_T *c_F;
  uint32_T h_y;
  emxArray_real_T *f_x0;
  emxArray_real_T *r12;
  emxArray_real_T *g_x0;
  emxArray_real_T *h_x0;
  emxArray_real_T *b_J;

  /*  Edit of fitNGaussians2D to work in 3D */
  /*    EHarry March 2012 */
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
  /* % Input */
  /*  %check whether correct number of input arguments was used */
  /*  if nargin ~= 4 */
  /*      disp('--fitNGaussians2D: Incorrect number of input arguments!'); */
  /*      return */
  /*  end */
  /* check whether correct number of input arguments was used */
  /* % Calculating F & J */
  /* extract background intensity from x0 and remove from vector */
  bgAmp = x0->data[x0->size[0] - 1];
  if (1 > x0->size[0] - 1) {
    i0 = -1;
  } else {
    i0 = x0->size[0] - 2;
  }

  c_emxInit_real_T(&b_x0, 1);
  i1 = b_x0->size[0];
  b_x0->size[0] = i0 + 1;
  emxEnsureCapacity((emxArray__common *)b_x0, i1, (int32_T)sizeof(real_T));
  for (i1 = 0; i1 <= i0; i1++) {
    b_x0->data[i1] = x0->data[i1];
  }

  i0 = x0->size[0];
  x0->size[0] = b_x0->size[0];
  emxEnsureCapacity((emxArray__common *)x0, i0, (int32_T)sizeof(real_T));
  loop_ub = b_x0->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    x0->data[i0] = b_x0->data[i0];
  }

  emxFree_real_T(&b_x0);

  /* get number of PSFs considered */
  /*  numPSF = length(x0)/3; */
  numPSF = (real_T)x0->size[0] / 4.0;

  /*  %reshape 3nx1 vector x0 into nx3 matrix */
  /*  x0 = reshape(x0,3,numPSF); */
  /*  x0 = x0'; */
  /* reshape 4nx1 vector x0 into nx4 matrix */
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
  emxEnsureCapacity((emxArray__common *)c_x0, i0, (int32_T)sizeof(real_T));
  for (k = 0; k + 1 <= nx; k++) {
    c_x0->data[k] = x0->data[k];
  }

  emxInit_real_T(&d_x0, 2);
  i0 = d_x0->size[0] * d_x0->size[1];
  d_x0->size[0] = c_x0->size[1];
  d_x0->size[1] = 4;
  emxEnsureCapacity((emxArray__common *)d_x0, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < 4; i0++) {
    loop_ub = c_x0->size[1] - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      d_x0->data[i1 + d_x0->size[0] * i0] = c_x0->data[i0 + c_x0->size[0] * i1];
    }
  }

  emxFree_real_T(&c_x0);
  emxInit_int32_T(&r0, 1);

  /* extract PSF center positions and amplitudes */
  /*  psfPos = x0(:,1:2); */
  /*  psfAmp = x0(:,3); */
  /* find minimum and maximum pixel indices */
  ixstart = 1;
  i0 = b_index->size[0];
  i1 = r0->size[0];
  r0->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r0, i1, (int32_T)sizeof(int32_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r0->data[i0] = 1 + i0;
  }

  nx = r0->size[0];
  mtmp = b_index->data[0];
  emxFree_int32_T(&r0);
  if (nx > 1) {
    if (rtIsNaN(b_index->data[0])) {
      nm1 = 2;
      exitg7 = FALSE;
      while ((exitg7 == 0U) && (nm1 <= nx)) {
        ixstart = nm1;
        if (!rtIsNaN(b_index->data[nm1 - 1])) {
          mtmp = b_index->data[nm1 - 1];
          exitg7 = TRUE;
        } else {
          nm1++;
        }
      }
    }

    if (ixstart < nx) {
      while (ixstart + 1 <= nx) {
        if (b_index->data[ixstart] < mtmp) {
          mtmp = b_index->data[ixstart];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r1, 1);
  ixstart = 1;
  i0 = b_index->size[0];
  i1 = r1->size[0];
  r1->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r1, i1, (int32_T)sizeof(int32_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r1->data[i0] = 1 + i0;
  }

  nx = r1->size[0];
  b_mtmp = b_index->data[0];
  emxFree_int32_T(&r1);
  if (nx > 1) {
    if (rtIsNaN(b_index->data[0])) {
      nm1 = 2;
      exitg6 = FALSE;
      while ((exitg6 == 0U) && (nm1 <= nx)) {
        ixstart = nm1;
        if (!rtIsNaN(b_index->data[nm1 - 1])) {
          b_mtmp = b_index->data[nm1 - 1];
          exitg6 = TRUE;
        } else {
          nm1++;
        }
      }
    }

    if (ixstart < nx) {
      while (ixstart + 1 <= nx) {
        if (b_index->data[ixstart] > b_mtmp) {
          b_mtmp = b_index->data[ixstart];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r2, 1);
  ixstart = 1;
  i0 = b_index->size[0];
  i1 = r2->size[0];
  r2->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r2, i1, (int32_T)sizeof(int32_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r2->data[i0] = 1 + i0;
  }

  nx = r2->size[0];
  c_mtmp = b_index->data[b_index->size[0]];
  emxFree_int32_T(&r2);
  if (nx > 1) {
    if (rtIsNaN(b_index->data[b_index->size[0]])) {
      nm1 = 2;
      exitg5 = FALSE;
      while ((exitg5 == 0U) && (nm1 <= nx)) {
        ixstart = nm1;
        if (!rtIsNaN(b_index->data[(nm1 + b_index->size[0]) - 1])) {
          c_mtmp = b_index->data[(nm1 + b_index->size[0]) - 1];
          exitg5 = TRUE;
        } else {
          nm1++;
        }
      }
    }

    if (ixstart < nx) {
      while (ixstart + 1 <= nx) {
        if (b_index->data[ixstart + b_index->size[0]] < c_mtmp) {
          c_mtmp = b_index->data[ixstart + b_index->size[0]];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r3, 1);
  ixstart = 1;
  i0 = b_index->size[0];
  i1 = r3->size[0];
  r3->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r3, i1, (int32_T)sizeof(int32_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r3->data[i0] = 1 + i0;
  }

  nx = r3->size[0];
  d_mtmp = b_index->data[b_index->size[0]];
  emxFree_int32_T(&r3);
  if (nx > 1) {
    if (rtIsNaN(b_index->data[b_index->size[0]])) {
      nm1 = 2;
      exitg4 = FALSE;
      while ((exitg4 == 0U) && (nm1 <= nx)) {
        ixstart = nm1;
        if (!rtIsNaN(b_index->data[(nm1 + b_index->size[0]) - 1])) {
          d_mtmp = b_index->data[(nm1 + b_index->size[0]) - 1];
          exitg4 = TRUE;
        } else {
          nm1++;
        }
      }
    }

    if (ixstart < nx) {
      while (ixstart + 1 <= nx) {
        if (b_index->data[ixstart + b_index->size[0]] > d_mtmp) {
          d_mtmp = b_index->data[ixstart + b_index->size[0]];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r4, 1);
  ixstart = 1;
  i0 = b_index->size[0];
  i1 = r4->size[0];
  r4->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r4, i1, (int32_T)sizeof(int32_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r4->data[i0] = 1 + i0;
  }

  nx = r4->size[0];
  e_mtmp = b_index->data[b_index->size[0] << 1];
  emxFree_int32_T(&r4);
  if (nx > 1) {
    if (rtIsNaN(b_index->data[b_index->size[0] << 1])) {
      nm1 = 2;
      exitg3 = FALSE;
      while ((exitg3 == 0U) && (nm1 <= nx)) {
        ixstart = nm1;
        if (!rtIsNaN(b_index->data[(nm1 + (b_index->size[0] << 1)) - 1])) {
          e_mtmp = b_index->data[(nm1 + (b_index->size[0] << 1)) - 1];
          exitg3 = TRUE;
        } else {
          nm1++;
        }
      }
    }

    if (ixstart < nx) {
      while (ixstart + 1 <= nx) {
        if (b_index->data[ixstart + (b_index->size[0] << 1)] < e_mtmp) {
          e_mtmp = b_index->data[ixstart + (b_index->size[0] << 1)];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r5, 1);
  ixstart = 1;
  i0 = b_index->size[0];
  i1 = r5->size[0];
  r5->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r5, i1, (int32_T)sizeof(int32_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r5->data[i0] = 1 + i0;
  }

  nx = r5->size[0];
  f_mtmp = b_index->data[b_index->size[0] << 1];
  emxFree_int32_T(&r5);
  if (nx > 1) {
    if (rtIsNaN(b_index->data[b_index->size[0] << 1])) {
      nm1 = 2;
      exitg2 = FALSE;
      while ((exitg2 == 0U) && (nm1 <= nx)) {
        ixstart = nm1;
        if (!rtIsNaN(b_index->data[(nm1 + (b_index->size[0] << 1)) - 1])) {
          f_mtmp = b_index->data[(nm1 + (b_index->size[0] << 1)) - 1];
          exitg2 = TRUE;
        } else {
          nm1++;
        }
      }
    }

    if (ixstart < nx) {
      while (ixstart + 1 <= nx) {
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
  i0 = psfIntegX->size[0] * psfIntegX->size[1];
  psfIntegX->size[0] = (int32_T)((b_mtmp - mtmp) + 1.0);
  psfIntegX->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfIntegX, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((b_mtmp - mtmp) + 1.0) * (int32_T)numPSF - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfIntegX->data[i0] = 0.0;
  }

  i = 0;
  b_emxInit_real_T(&temp, 3);
  emxInit_real_T(&temp2, 2);
  emxInit_real_T(&y, 2);
  c_emxInit_real_T(&b_y, 1);
  while (i <= (int32_T)numPSF - 1) {
    if (rtIsNaN(mtmp) || rtIsNaN(b_mtmp)) {
      nx = 1;
      anew = rtNaN;
      apnd = b_mtmp;
    } else if (b_mtmp < mtmp) {
      nx = 0;
      anew = mtmp;
      apnd = b_mtmp;
    } else if (rtIsInf(mtmp) || rtIsInf(b_mtmp)) {
      nx = 1;
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

      if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
        ndbl++;
        apnd = b_mtmp;
      } else if (cdiff > 0.0) {
        apnd = mtmp + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        nx = (int32_T)ndbl;
      } else {
        nx = 0;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = nx;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (nx > 0) {
      y->data[0] = anew;
      if (nx > 1) {
        y->data[nx - 1] = apnd;
        nm1 = nx - 1;
        ixstart = nm1 / 2;
        for (k = 1; k <= ixstart - 1; k++) {
          y->data[k] = anew + (real_T)k;
          y->data[(nx - k) - 1] = apnd - (real_T)k;
        }

        if (ixstart << 1 == nm1) {
          y->data[ixstart] = (anew + apnd) / 2.0;
        } else {
          y->data[ixstart] = anew + (real_T)ixstart;
          y->data[ixstart + 1] = apnd - (real_T)ixstart;
        }
      }
    }

    i0 = b_y->size[0];
    b_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)b_y, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      b_y->data[i0] = y->data[i0];
    }

    GaussListND_mexCode(b_y, psfSigma[0], d_x0->data[i], temp);
    squeeze(temp, temp2);

    /* clear temp */
    i0 = temp2->size[0] - 1;
    for (i1 = 0; i1 <= i0; i1++) {
      psfIntegX->data[i1 + psfIntegX->size[0] * i] = temp2->data[i1];
    }

    /* clear temp2 */
    i++;
  }

  emxFree_real_T(&b_y);
  emxInit_real_T(&psfIntegY, 2);

  /* determine the contribution of each PSF (assuming amplitude 1) to a */
  /* pixel based on its y-coordinate (needed to calculate F & J) */
  i0 = psfIntegY->size[0] * psfIntegY->size[1];
  psfIntegY->size[0] = (int32_T)((d_mtmp - c_mtmp) + 1.0);
  psfIntegY->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfIntegY, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((d_mtmp - c_mtmp) + 1.0) * (int32_T)numPSF - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfIntegY->data[i0] = 0.0;
  }

  i = 0;
  c_emxInit_real_T(&c_y, 1);
  while (i <= (int32_T)numPSF - 1) {
    if (rtIsNaN(c_mtmp) || rtIsNaN(d_mtmp)) {
      nx = 1;
      anew = rtNaN;
      apnd = d_mtmp;
    } else if (d_mtmp < c_mtmp) {
      nx = 0;
      anew = c_mtmp;
      apnd = d_mtmp;
    } else if (rtIsInf(c_mtmp) || rtIsInf(d_mtmp)) {
      nx = 1;
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

      if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
        ndbl++;
        apnd = d_mtmp;
      } else if (cdiff > 0.0) {
        apnd = c_mtmp + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        nx = (int32_T)ndbl;
      } else {
        nx = 0;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = nx;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (nx > 0) {
      y->data[0] = anew;
      if (nx > 1) {
        y->data[nx - 1] = apnd;
        nm1 = nx - 1;
        ixstart = nm1 / 2;
        for (k = 1; k <= ixstart - 1; k++) {
          y->data[k] = anew + (real_T)k;
          y->data[(nx - k) - 1] = apnd - (real_T)k;
        }

        if (ixstart << 1 == nm1) {
          y->data[ixstart] = (anew + apnd) / 2.0;
        } else {
          y->data[ixstart] = anew + (real_T)ixstart;
          y->data[ixstart + 1] = apnd - (real_T)ixstart;
        }
      }
    }

    i0 = c_y->size[0];
    c_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)c_y, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      c_y->data[i0] = y->data[i0];
    }

    GaussListND_mexCode(c_y, psfSigma[0], d_x0->data[i + d_x0->size[0]], temp);
    squeeze(temp, temp2);

    /* clear temp */
    i0 = temp2->size[0] - 1;
    for (i1 = 0; i1 <= i0; i1++) {
      psfIntegY->data[i1 + psfIntegY->size[0] * i] = temp2->data[i1];
    }

    /* clear temp2 */
    i++;
  }

  emxFree_real_T(&c_y);
  emxInit_real_T(&psfIntegZ, 2);

  /* determine the contribution of each PSF (assuming amplitude 1) to a */
  /* pixel based on its z-coordinate (needed to calculate F & J) */
  i0 = psfIntegZ->size[0] * psfIntegZ->size[1];
  psfIntegZ->size[0] = (int32_T)((f_mtmp - e_mtmp) + 1.0);
  psfIntegZ->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfIntegZ, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((f_mtmp - e_mtmp) + 1.0) * (int32_T)numPSF - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfIntegZ->data[i0] = 0.0;
  }

  i = 0;
  c_emxInit_real_T(&d_y, 1);
  while (i <= (int32_T)numPSF - 1) {
    if (rtIsNaN(e_mtmp) || rtIsNaN(f_mtmp)) {
      nx = 1;
      anew = rtNaN;
      apnd = f_mtmp;
    } else if (f_mtmp < e_mtmp) {
      nx = 0;
      anew = e_mtmp;
      apnd = f_mtmp;
    } else if (rtIsInf(e_mtmp) || rtIsInf(f_mtmp)) {
      nx = 1;
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

      if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
        ndbl++;
        apnd = f_mtmp;
      } else if (cdiff > 0.0) {
        apnd = e_mtmp + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        nx = (int32_T)ndbl;
      } else {
        nx = 0;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = nx;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (nx > 0) {
      y->data[0] = anew;
      if (nx > 1) {
        y->data[nx - 1] = apnd;
        nm1 = nx - 1;
        ixstart = nm1 / 2;
        for (k = 1; k <= ixstart - 1; k++) {
          y->data[k] = anew + (real_T)k;
          y->data[(nx - k) - 1] = apnd - (real_T)k;
        }

        if (ixstart << 1 == nm1) {
          y->data[ixstart] = (anew + apnd) / 2.0;
        } else {
          y->data[ixstart] = anew + (real_T)ixstart;
          y->data[ixstart + 1] = apnd - (real_T)ixstart;
        }
      }
    }

    i0 = d_y->size[0];
    d_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)d_y, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      d_y->data[i0] = y->data[i0];
    }

    GaussListND_mexCode(d_y, psfSigma[1], d_x0->data[i + (d_x0->size[0] << 1)],
                        temp);
    squeeze(temp, temp2);

    /* clear temp */
    i0 = temp2->size[0] - 1;
    for (i1 = 0; i1 <= i0; i1++) {
      psfIntegZ->data[i1 + psfIntegZ->size[0] * i] = temp2->data[i1];
    }

    /* clear temp2 */
    i++;
  }

  emxFree_real_T(&d_y);
  emxFree_real_T(&temp2);
  emxFree_real_T(&temp);
  emxInit_real_T(&psfValueX, 2);

  /* calculate the value of each PSF (assuming amplitude 1) at the */
  /* x-coordinates of the corners of all pixels (needed to calculate J) */
  i0 = psfValueX->size[0] * psfValueX->size[1];
  psfValueX->size[0] = (int32_T)((b_mtmp - mtmp) + 2.0);
  psfValueX->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfValueX, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((b_mtmp - mtmp) + 2.0) * (int32_T)numPSF - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfValueX->data[i0] = 0.0;
  }

  i = 0;
  c_emxInit_real_T(&indxPixel, 1);
  c_emxInit_real_T(&b_indxPixel, 1);
  c_emxInit_real_T(&c_indxPixel, 1);
  c_emxInit_real_T(&e_y, 1);
  while (i <= (int32_T)numPSF - 1) {
    if (rtIsNaN(mtmp - 0.5) || rtIsNaN(b_mtmp + 0.5)) {
      nx = 1;
      anew = rtNaN;
      apnd = b_mtmp + 0.5;
    } else if (b_mtmp + 0.5 < mtmp - 0.5) {
      nx = 0;
      anew = mtmp - 0.5;
      apnd = b_mtmp + 0.5;
    } else if (rtIsInf(mtmp - 0.5) || rtIsInf(b_mtmp + 0.5)) {
      nx = 1;
      anew = rtNaN;
      apnd = b_mtmp + 0.5;
    } else {
      anew = mtmp - 0.5;
      ndbl = floor(((b_mtmp + 0.5) - (mtmp - 0.5)) + 0.5);
      apnd = (mtmp - 0.5) + ndbl;
      cdiff = apnd - (b_mtmp + 0.5);
      absa = fabs(mtmp - 0.5);
      absb = fabs(b_mtmp + 0.5);
      if (absa > absb) {
        absb = absa;
      }

      if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
        ndbl++;
        apnd = b_mtmp + 0.5;
      } else if (cdiff > 0.0) {
        apnd = (mtmp - 0.5) + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        nx = (int32_T)ndbl;
      } else {
        nx = 0;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = nx;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (nx > 0) {
      y->data[0] = anew;
      if (nx > 1) {
        y->data[nx - 1] = apnd;
        nm1 = nx - 1;
        ixstart = nm1 / 2;
        for (k = 1; k <= ixstart - 1; k++) {
          y->data[k] = anew + (real_T)k;
          y->data[(nx - k) - 1] = apnd - (real_T)k;
        }

        if (ixstart << 1 == nm1) {
          y->data[ixstart] = (anew + apnd) / 2.0;
        } else {
          y->data[ixstart] = anew + (real_T)ixstart;
          y->data[ixstart + 1] = apnd - (real_T)ixstart;
        }
      }
    }

    anew = d_x0->data[i];
    i0 = e_y->size[0];
    e_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)e_y, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      e_y->data[i0] = y->data[i0] - anew;
    }

    power(e_y, indxPixel);
    i0 = indxPixel->size[0];
    indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)indxPixel, i0, (int32_T)sizeof(real_T));
    loop_ub = indxPixel->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      indxPixel->data[i0] = -indxPixel->data[i0];
    }

    i0 = c_indxPixel->size[0];
    c_indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)c_indxPixel, i0, (int32_T)sizeof
                      (real_T));
    loop_ub = indxPixel->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      c_indxPixel->data[i0] = indxPixel->data[i0];
    }

    rdivide(c_indxPixel, 2.0, indxPixel);
    i0 = b_indxPixel->size[0];
    b_indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)b_indxPixel, i0, (int32_T)sizeof
                      (real_T));
    loop_ub = indxPixel->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      b_indxPixel->data[i0] = indxPixel->data[i0];
    }

    rdivide(b_indxPixel, rt_powd_snf(psfSigma[0], 2.0), indxPixel);
    b_exp(indxPixel);
    loop_ub = indxPixel->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      psfValueX->data[i0 + psfValueX->size[0] * i] = indxPixel->data[i0];
    }

    i++;
  }

  emxFree_real_T(&e_y);
  emxFree_real_T(&c_indxPixel);
  emxFree_real_T(&b_indxPixel);
  emxInit_real_T(&psfValueY, 2);

  /* calculate the value of each PSF (assuming amplitude 1) at the */
  /* y-coordinates of the corners of all pixels (needed to calculate J) */
  i0 = psfValueY->size[0] * psfValueY->size[1];
  psfValueY->size[0] = (int32_T)((d_mtmp - c_mtmp) + 2.0);
  psfValueY->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfValueY, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((d_mtmp - c_mtmp) + 2.0) * (int32_T)numPSF - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfValueY->data[i0] = 0.0;
  }

  i = 0;
  c_emxInit_real_T(&d_indxPixel, 1);
  c_emxInit_real_T(&e_indxPixel, 1);
  c_emxInit_real_T(&f_y, 1);
  while (i <= (int32_T)numPSF - 1) {
    if (rtIsNaN(c_mtmp - 0.5) || rtIsNaN(d_mtmp + 0.5)) {
      nx = 1;
      anew = rtNaN;
      apnd = d_mtmp + 0.5;
    } else if (d_mtmp + 0.5 < c_mtmp - 0.5) {
      nx = 0;
      anew = c_mtmp - 0.5;
      apnd = d_mtmp + 0.5;
    } else if (rtIsInf(c_mtmp - 0.5) || rtIsInf(d_mtmp + 0.5)) {
      nx = 1;
      anew = rtNaN;
      apnd = d_mtmp + 0.5;
    } else {
      anew = c_mtmp - 0.5;
      ndbl = floor(((d_mtmp + 0.5) - (c_mtmp - 0.5)) + 0.5);
      apnd = (c_mtmp - 0.5) + ndbl;
      cdiff = apnd - (d_mtmp + 0.5);
      absa = fabs(c_mtmp - 0.5);
      absb = fabs(d_mtmp + 0.5);
      if (absa > absb) {
        absb = absa;
      }

      if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
        ndbl++;
        apnd = d_mtmp + 0.5;
      } else if (cdiff > 0.0) {
        apnd = (c_mtmp - 0.5) + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        nx = (int32_T)ndbl;
      } else {
        nx = 0;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = nx;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (nx > 0) {
      y->data[0] = anew;
      if (nx > 1) {
        y->data[nx - 1] = apnd;
        nm1 = nx - 1;
        ixstart = nm1 / 2;
        for (k = 1; k <= ixstart - 1; k++) {
          y->data[k] = anew + (real_T)k;
          y->data[(nx - k) - 1] = apnd - (real_T)k;
        }

        if (ixstart << 1 == nm1) {
          y->data[ixstart] = (anew + apnd) / 2.0;
        } else {
          y->data[ixstart] = anew + (real_T)ixstart;
          y->data[ixstart + 1] = apnd - (real_T)ixstart;
        }
      }
    }

    anew = d_x0->data[i + d_x0->size[0]];
    i0 = f_y->size[0];
    f_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)f_y, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      f_y->data[i0] = y->data[i0] - anew;
    }

    power(f_y, indxPixel);
    i0 = indxPixel->size[0];
    indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)indxPixel, i0, (int32_T)sizeof(real_T));
    loop_ub = indxPixel->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      indxPixel->data[i0] = -indxPixel->data[i0];
    }

    i0 = e_indxPixel->size[0];
    e_indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)e_indxPixel, i0, (int32_T)sizeof
                      (real_T));
    loop_ub = indxPixel->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      e_indxPixel->data[i0] = indxPixel->data[i0];
    }

    rdivide(e_indxPixel, 2.0, indxPixel);
    i0 = d_indxPixel->size[0];
    d_indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)d_indxPixel, i0, (int32_T)sizeof
                      (real_T));
    loop_ub = indxPixel->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      d_indxPixel->data[i0] = indxPixel->data[i0];
    }

    rdivide(d_indxPixel, rt_powd_snf(psfSigma[0], 2.0), indxPixel);
    b_exp(indxPixel);
    loop_ub = indxPixel->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      psfValueY->data[i0 + psfValueY->size[0] * i] = indxPixel->data[i0];
    }

    i++;
  }

  emxFree_real_T(&f_y);
  emxFree_real_T(&e_indxPixel);
  emxFree_real_T(&d_indxPixel);
  emxInit_real_T(&psfValueZ, 2);

  /* calculate the value of each PSF (assuming amplitude 1) at the */
  /* z-coordinates of the corners of all pixels (needed to calculate J) */
  i0 = psfValueZ->size[0] * psfValueZ->size[1];
  psfValueZ->size[0] = (int32_T)((f_mtmp - e_mtmp) + 2.0);
  psfValueZ->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfValueZ, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((f_mtmp - e_mtmp) + 2.0) * (int32_T)numPSF - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfValueZ->data[i0] = 0.0;
  }

  i = 0;
  c_emxInit_real_T(&f_indxPixel, 1);
  c_emxInit_real_T(&g_indxPixel, 1);
  c_emxInit_real_T(&g_y, 1);
  while (i <= (int32_T)numPSF - 1) {
    if (rtIsNaN(e_mtmp - 0.5) || rtIsNaN(f_mtmp + 0.5)) {
      nx = 1;
      anew = rtNaN;
      apnd = f_mtmp + 0.5;
    } else if (f_mtmp + 0.5 < e_mtmp - 0.5) {
      nx = 0;
      anew = e_mtmp - 0.5;
      apnd = f_mtmp + 0.5;
    } else if (rtIsInf(e_mtmp - 0.5) || rtIsInf(f_mtmp + 0.5)) {
      nx = 1;
      anew = rtNaN;
      apnd = f_mtmp + 0.5;
    } else {
      anew = e_mtmp - 0.5;
      ndbl = floor(((f_mtmp + 0.5) - (e_mtmp - 0.5)) + 0.5);
      apnd = (e_mtmp - 0.5) + ndbl;
      cdiff = apnd - (f_mtmp + 0.5);
      absa = fabs(e_mtmp - 0.5);
      absb = fabs(f_mtmp + 0.5);
      if (absa > absb) {
        absb = absa;
      }

      if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
        ndbl++;
        apnd = f_mtmp + 0.5;
      } else if (cdiff > 0.0) {
        apnd = (e_mtmp - 0.5) + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        nx = (int32_T)ndbl;
      } else {
        nx = 0;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = nx;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (nx > 0) {
      y->data[0] = anew;
      if (nx > 1) {
        y->data[nx - 1] = apnd;
        nm1 = nx - 1;
        ixstart = nm1 / 2;
        for (k = 1; k <= ixstart - 1; k++) {
          y->data[k] = anew + (real_T)k;
          y->data[(nx - k) - 1] = apnd - (real_T)k;
        }

        if (ixstart << 1 == nm1) {
          y->data[ixstart] = (anew + apnd) / 2.0;
        } else {
          y->data[ixstart] = anew + (real_T)ixstart;
          y->data[ixstart + 1] = apnd - (real_T)ixstart;
        }
      }
    }

    anew = d_x0->data[i + (d_x0->size[0] << 1)];
    i0 = g_y->size[0];
    g_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)g_y, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      g_y->data[i0] = y->data[i0] - anew;
    }

    power(g_y, indxPixel);
    i0 = indxPixel->size[0];
    indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)indxPixel, i0, (int32_T)sizeof(real_T));
    loop_ub = indxPixel->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      indxPixel->data[i0] = -indxPixel->data[i0];
    }

    i0 = g_indxPixel->size[0];
    g_indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)g_indxPixel, i0, (int32_T)sizeof
                      (real_T));
    loop_ub = indxPixel->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      g_indxPixel->data[i0] = indxPixel->data[i0];
    }

    rdivide(g_indxPixel, 2.0, indxPixel);
    i0 = f_indxPixel->size[0];
    f_indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)f_indxPixel, i0, (int32_T)sizeof
                      (real_T));
    loop_ub = indxPixel->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      f_indxPixel->data[i0] = indxPixel->data[i0];
    }

    rdivide(f_indxPixel, rt_powd_snf(psfSigma[1], 2.0), indxPixel);
    b_exp(indxPixel);
    loop_ub = indxPixel->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      psfValueZ->data[i0 + psfValueZ->size[0] * i] = indxPixel->data[i0];
    }

    i++;
  }

  emxFree_real_T(&g_y);
  emxFree_real_T(&g_indxPixel);
  emxFree_real_T(&f_indxPixel);
  c_emxInit_real_T(&relIndxX, 1);

  /* get number of pixels in image */
  /*  %get xy-indices relative to minimum */
  /*  relIndxX = index(:,1) - minIndxX + 1; */
  /*  relIndxY = index(:,2) - minIndxY + 1; */
  /* get xy-indices relative to minimum */
  i0 = b_index->size[0];
  i1 = relIndxX->size[0];
  relIndxX->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)relIndxX, i1, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    relIndxX->data[i0] = (b_index->data[i0] - mtmp) + 1.0;
  }

  c_emxInit_real_T(&relIndxY, 1);
  i0 = b_index->size[0];
  i1 = relIndxY->size[0];
  relIndxY->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)relIndxY, i1, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    relIndxY->data[i0] = (b_index->data[i0 + b_index->size[0]] - c_mtmp) + 1.0;
  }

  c_emxInit_real_T(&relIndxZ, 1);
  i0 = b_index->size[0];
  i1 = relIndxZ->size[0];
  relIndxZ->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)relIndxZ, i1, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    relIndxZ->data[i0] = (b_index->data[i0 + (b_index->size[0] << 1)] - e_mtmp)
      + 1.0;
  }

  c_emxInit_real_T(&e_x0, 1);

  /*  %calculate the value of F at all pixels */
  /*  F = (sum(repmat(psfAmp,1,numPixel).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)',1))' ... */
  /*      + repmat(bgAmp,numPixel,1) - image; */
  /* calculate the value of F at all pixels */
  i0 = d_x0->size[0];
  i1 = e_x0->size[0];
  e_x0->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)e_x0, i1, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    e_x0->data[i0] = d_x0->data[i0 + d_x0->size[0] * 3];
  }

  emxInit_real_T(&r6, 2);
  emxInit_real_T(&r7, 2);
  repmat(e_x0, (real_T)image->size[0], r6);
  i0 = psfIntegX->size[1];
  i1 = r7->size[0] * r7->size[1];
  r7->size[0] = i0;
  r7->size[1] = relIndxX->size[0];
  emxEnsureCapacity((emxArray__common *)r7, i1, (int32_T)sizeof(real_T));
  emxFree_real_T(&e_x0);
  loop_ub = relIndxX->size[0] - 1;
  for (i1 = 0; i1 <= loop_ub; i1++) {
    nx = i0 - 1;
    for (nm1 = 0; nm1 <= nx; nm1++) {
      r7->data[nm1 + r7->size[0] * i1] = psfIntegX->data[((int32_T)
        relIndxX->data[i1] + psfIntegX->size[0] * nm1) - 1];
    }
  }

  emxInit_real_T(&r8, 2);
  i0 = psfIntegY->size[1];
  i1 = r8->size[0] * r8->size[1];
  r8->size[0] = i0;
  r8->size[1] = relIndxY->size[0];
  emxEnsureCapacity((emxArray__common *)r8, i1, (int32_T)sizeof(real_T));
  loop_ub = relIndxY->size[0] - 1;
  for (i1 = 0; i1 <= loop_ub; i1++) {
    nx = i0 - 1;
    for (nm1 = 0; nm1 <= nx; nm1++) {
      r8->data[nm1 + r8->size[0] * i1] = psfIntegY->data[((int32_T)
        relIndxY->data[i1] + psfIntegY->size[0] * nm1) - 1];
    }
  }

  emxInit_real_T(&r9, 2);
  i0 = psfIntegZ->size[1];
  i1 = r9->size[0] * r9->size[1];
  r9->size[0] = i0;
  r9->size[1] = relIndxZ->size[0];
  emxEnsureCapacity((emxArray__common *)r9, i1, (int32_T)sizeof(real_T));
  loop_ub = relIndxZ->size[0] - 1;
  for (i1 = 0; i1 <= loop_ub; i1++) {
    nx = i0 - 1;
    for (nm1 = 0; nm1 <= nx; nm1++) {
      r9->data[nm1 + r9->size[0] * i1] = psfIntegZ->data[((int32_T)
        relIndxZ->data[i1] + psfIntegZ->size[0] * nm1) - 1];
    }
  }

  emxInit_real_T(&r10, 2);
  i0 = r10->size[0] * r10->size[1];
  r10->size[0] = r6->size[0];
  r10->size[1] = r6->size[1];
  emxEnsureCapacity((emxArray__common *)r10, i0, (int32_T)sizeof(real_T));
  loop_ub = r6->size[0] * r6->size[1] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r10->data[i0] = r6->data[i0] * r7->data[i0] * r8->data[i0] * r9->data[i0];
  }

  emxFree_real_T(&r9);
  emxFree_real_T(&r8);
  emxFree_real_T(&r7);
  emxFree_real_T(&r6);
  sum(r10, y);
  i0 = indxPixel->size[0];
  indxPixel->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)indxPixel, i0, (int32_T)sizeof(real_T));
  emxFree_real_T(&r10);
  loop_ub = y->size[1] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    indxPixel->data[i0] = y->data[i0];
  }

  emxFree_real_T(&y);
  c_emxInit_real_T(&r11, 1);
  b_repmat(bgAmp, (real_T)image->size[0], r11);
  nm1 = indxPixel->size[0];
  i0 = F->size[0] * F->size[1];
  F->size[0] = nm1;
  emxEnsureCapacity((emxArray__common *)F, i0, (int32_T)sizeof(real_T));
  i0 = F->size[0] * F->size[1];
  F->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)F, i0, (int32_T)sizeof(real_T));
  loop_ub = indxPixel->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    F->data[i0] = (indxPixel->data[i0] + r11->data[i0]) - image->data[i0];
  }

  emxFree_real_T(&r11);
  emxInit_boolean_T(&x, 1);

  /* remove pixels with NaN (which means they are out of the cropped image */
  /* area) */
  i0 = x->size[0];
  x->size[0] = image->size[0];
  emxEnsureCapacity((emxArray__common *)x, i0, (int32_T)sizeof(boolean_T));
  loop_ub = image->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    x->data[i0] = rtIsNaN(image->data[i0]);
  }

  i0 = x->size[0];
  x->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)x, i0, (int32_T)sizeof(boolean_T));
  loop_ub = x->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    x->data[i0] = !x->data[i0];
  }

  emxInit_int32_T(&ii, 1);
  nx = x->size[0];
  ixstart = 0;
  i0 = ii->size[0];
  ii->size[0] = nx;
  emxEnsureCapacity((emxArray__common *)ii, i0, (int32_T)sizeof(int32_T));
  nm1 = 1;
  exitg1 = FALSE;
  while ((exitg1 == 0U) && (nm1 <= nx)) {
    guard1 = FALSE;
    if (x->data[nm1 - 1]) {
      ixstart++;
      ii->data[ixstart - 1] = nm1;
      if (ixstart >= nx) {
        exitg1 = TRUE;
      } else {
        guard1 = TRUE;
      }
    } else {
      guard1 = TRUE;
    }

    if (guard1 == TRUE) {
      nm1++;
    }
  }

  emxFree_boolean_T(&x);
  if (nx == 1) {
    if (ixstart == 0) {
      i0 = ii->size[0];
      ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)ii, i0, (int32_T)sizeof(int32_T));
    }
  } else {
    if (1 > ixstart) {
      ixstart = 0;
    }

    emxInit_int32_T(&b_ii, 1);
    i0 = b_ii->size[0];
    b_ii->size[0] = ixstart;
    emxEnsureCapacity((emxArray__common *)b_ii, i0, (int32_T)sizeof(int32_T));
    loop_ub = ixstart - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      b_ii->data[i0] = ii->data[i0];
    }

    i0 = ii->size[0];
    ii->size[0] = b_ii->size[0];
    emxEnsureCapacity((emxArray__common *)ii, i0, (int32_T)sizeof(int32_T));
    loop_ub = b_ii->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      ii->data[i0] = b_ii->data[i0];
    }

    emxFree_int32_T(&b_ii);
  }

  i0 = indxPixel->size[0];
  indxPixel->size[0] = ii->size[0];
  emxEnsureCapacity((emxArray__common *)indxPixel, i0, (int32_T)sizeof(real_T));
  loop_ub = ii->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    indxPixel->data[i0] = (real_T)ii->data[i0];
  }

  emxFree_int32_T(&ii);
  c_emxInit_real_T(&b_F, 1);
  i0 = b_F->size[0];
  b_F->size[0] = indxPixel->size[0];
  emxEnsureCapacity((emxArray__common *)b_F, i0, (int32_T)sizeof(real_T));
  loop_ub = indxPixel->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    b_F->data[i0] = F->data[(int32_T)indxPixel->data[i0] - 1];
  }

  c_emxInit_real_T(&c_F, 1);
  i0 = c_F->size[0];
  c_F->size[0] = indxPixel->size[0];
  emxEnsureCapacity((emxArray__common *)c_F, i0, (int32_T)sizeof(real_T));
  loop_ub = indxPixel->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    c_F->data[i0] = F->data[(int32_T)indxPixel->data[i0] - 1];
  }

  nm1 = b_F->size[0];
  i0 = F->size[0] * F->size[1];
  F->size[0] = nm1;
  F->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)F, i0, (int32_T)sizeof(real_T));
  emxFree_real_T(&b_F);
  i0 = 0;
  while (i0 <= 0) {
    loop_ub = nm1 - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      F->data[i0] = c_F->data[i0];
    }

    i0 = 1;
  }

  emxFree_real_T(&c_F);

  /* calculate the derivative at all pixels */
  /*      J = ones(numPixel,3*numPSF+1); %(last column for background amplitude) */
  /*      J(:,1:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueX(relIndxX,:)-... */
  /*          psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:); %w.r.t. x */
  /*      J(:,2:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueY(relIndxY,:)-... */
  /*          psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:); %w.r.t. y */
  /*      J(:,3:3:3*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); %w.r.t. amp */
  h_y = (uint32_T)numPSF << 2;
  nm1 = image->size[0];
  i0 = J->size[0] * J->size[1];
  J->size[0] = nm1;
  emxEnsureCapacity((emxArray__common *)J, i0, (int32_T)sizeof(real_T));
  i0 = J->size[0] * J->size[1];
  J->size[1] = (int32_T)((real_T)h_y + 1.0);
  emxEnsureCapacity((emxArray__common *)J, i0, (int32_T)sizeof(real_T));
  loop_ub = image->size[0] * (int32_T)((real_T)h_y + 1.0) - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    J->data[i0] = 1.0;
  }

  emxInit_real_T(&f_x0, 2);

  /* (last column for background amplitude) */
  i0 = d_x0->size[0];
  i1 = f_x0->size[0] * f_x0->size[1];
  f_x0->size[0] = 1;
  f_x0->size[1] = i0;
  emxEnsureCapacity((emxArray__common *)f_x0, i1, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    f_x0->data[f_x0->size[0] * i0] = d_x0->data[i0 + d_x0->size[0] * 3];
  }

  emxInit_real_T(&r12, 2);
  c_repmat(f_x0, (real_T)image->size[0], r12);
  emxFree_real_T(&f_x0);
  if (1U > ((uint32_T)numPSF << 2)) {
    i0 = 1;
  } else {
    i0 = 4;
  }

  loop_ub = r12->size[1] - 1;
  for (i1 = 0; i1 <= loop_ub; i1++) {
    nx = r12->size[0] - 1;
    for (nm1 = 0; nm1 <= nx; nm1++) {
      J->data[nm1 + J->size[0] * (i0 * i1)] = r12->data[nm1 + r12->size[0] * i1]
        * (psfValueX->data[((int32_T)relIndxX->data[nm1] + psfValueX->size[0] *
                            i1) - 1] - psfValueX->data[((int32_T)(relIndxX->
             data[nm1] + 1.0) + psfValueX->size[0] * i1) - 1]) * psfIntegY->
        data[((int32_T)relIndxY->data[nm1] + psfIntegY->size[0] * i1) - 1] *
        psfIntegZ->data[((int32_T)relIndxZ->data[nm1] + psfIntegZ->size[0] * i1)
        - 1];
    }
  }

  emxFree_real_T(&psfValueX);
  emxInit_real_T(&g_x0, 2);

  /* w.r.t. x */
  i0 = d_x0->size[0];
  i1 = g_x0->size[0] * g_x0->size[1];
  g_x0->size[0] = 1;
  g_x0->size[1] = i0;
  emxEnsureCapacity((emxArray__common *)g_x0, i1, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    g_x0->data[g_x0->size[0] * i0] = d_x0->data[i0 + d_x0->size[0] * 3];
  }

  c_repmat(g_x0, (real_T)image->size[0], r12);
  emxFree_real_T(&g_x0);
  if (2U > ((uint32_T)numPSF << 2)) {
    i0 = 0;
    i1 = 1;
  } else {
    i0 = 1;
    i1 = 4;
  }

  loop_ub = r12->size[1] - 1;
  for (nm1 = 0; nm1 <= loop_ub; nm1++) {
    nx = r12->size[0] - 1;
    for (k = 0; k <= nx; k++) {
      J->data[k + J->size[0] * (i0 + i1 * nm1)] = r12->data[k + r12->size[0] *
        nm1] * (psfValueY->data[((int32_T)relIndxY->data[k] + psfValueY->size[0]
                 * nm1) - 1] - psfValueY->data[((int32_T)(relIndxY->data[k] +
                  1.0) + psfValueY->size[0] * nm1) - 1]) * psfIntegX->data
        [((int32_T)relIndxX->data[k] + psfIntegX->size[0] * nm1) - 1] *
        psfIntegZ->data[((int32_T)relIndxZ->data[k] + psfIntegZ->size[0] * nm1)
        - 1];
    }
  }

  emxFree_real_T(&psfValueY);
  emxInit_real_T(&h_x0, 2);

  /* w.r.t. y */
  i0 = d_x0->size[0];
  i1 = h_x0->size[0] * h_x0->size[1];
  h_x0->size[0] = 1;
  h_x0->size[1] = i0;
  emxEnsureCapacity((emxArray__common *)h_x0, i1, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    h_x0->data[h_x0->size[0] * i0] = d_x0->data[i0 + d_x0->size[0] * 3];
  }

  emxFree_real_T(&d_x0);
  c_repmat(h_x0, (real_T)image->size[0], r12);
  emxFree_real_T(&h_x0);
  if (3U > ((uint32_T)numPSF << 2)) {
    i0 = 0;
    i1 = 1;
  } else {
    i0 = 2;
    i1 = 4;
  }

  loop_ub = r12->size[1] - 1;
  for (nm1 = 0; nm1 <= loop_ub; nm1++) {
    nx = r12->size[0] - 1;
    for (k = 0; k <= nx; k++) {
      J->data[k + J->size[0] * (i0 + i1 * nm1)] = r12->data[k + r12->size[0] *
        nm1] * (psfValueZ->data[((int32_T)relIndxZ->data[k] + psfValueZ->size[0]
                 * nm1) - 1] - psfValueZ->data[((int32_T)(relIndxZ->data[k] +
                  1.0) + psfValueZ->size[0] * nm1) - 1]) * psfIntegX->data
        [((int32_T)relIndxX->data[k] + psfIntegX->size[0] * nm1) - 1] *
        psfIntegY->data[((int32_T)relIndxY->data[k] + psfIntegY->size[0] * nm1)
        - 1];
    }
  }

  emxFree_real_T(&r12);
  emxFree_real_T(&psfValueZ);

  /* w.r.t. z */
  if (4U > ((uint32_T)numPSF << 2)) {
    i0 = 0;
    i1 = 1;
  } else {
    i0 = 3;
    i1 = 4;
  }

  nm1 = psfIntegX->size[1] - 1;
  for (k = 0; k <= nm1; k++) {
    loop_ub = relIndxX->size[0] - 1;
    for (nx = 0; nx <= loop_ub; nx++) {
      J->data[nx + J->size[0] * (i0 + i1 * k)] = psfIntegX->data[((int32_T)
        relIndxX->data[nx] + psfIntegX->size[0] * k) - 1] * psfIntegY->data
        [((int32_T)relIndxY->data[nx] + psfIntegY->size[0] * k) - 1] *
        psfIntegZ->data[((int32_T)relIndxZ->data[nx] + psfIntegZ->size[0] * k) -
        1];
    }
  }

  emxFree_real_T(&relIndxZ);
  emxFree_real_T(&relIndxY);
  emxFree_real_T(&relIndxX);
  emxFree_real_T(&psfIntegZ);
  emxFree_real_T(&psfIntegY);
  emxFree_real_T(&psfIntegX);
  emxInit_real_T(&b_J, 2);

  /* w.r.t. amp */
  /* remove pixels with NaN (which means they are out of the cropped image */
  /* area) */
  nm1 = J->size[1];
  i0 = b_J->size[0] * b_J->size[1];
  b_J->size[0] = indxPixel->size[0];
  b_J->size[1] = nm1;
  emxEnsureCapacity((emxArray__common *)b_J, i0, (int32_T)sizeof(real_T));
  loop_ub = nm1 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    nx = indxPixel->size[0] - 1;
    for (i1 = 0; i1 <= nx; i1++) {
      b_J->data[i1 + b_J->size[0] * i0] = J->data[((int32_T)indxPixel->data[i1]
        + J->size[0] * i0) - 1];
    }
  }

  emxFree_real_T(&indxPixel);
  i0 = J->size[0] * J->size[1];
  J->size[0] = b_J->size[0];
  J->size[1] = b_J->size[1];
  emxEnsureCapacity((emxArray__common *)J, i0, (int32_T)sizeof(real_T));
  loop_ub = b_J->size[1] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    nx = b_J->size[0] - 1;
    for (i1 = 0; i1 <= nx; i1++) {
      J->data[i1 + J->size[0] * i0] = b_J->data[i1 + b_J->size[0] * i0];
    }
  }

  emxFree_real_T(&b_J);

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

/* End of code generation (fitNGaussians3D_mexCode.c) */
