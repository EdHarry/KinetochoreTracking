/*
 * fitNGaussians3D_mexCode.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode'
 *
 * C source code generated on: Tue Nov 19 11:16:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "GaussListND_mexCode.h"
#include "fitNGaussians3D_mexCode_emxutil.h"
#include "exp.h"
#include "rdivide.h"
#include "power.h"
#include "repmat.h"
#include "sum.h"

/* Function Definitions */
void fitNGaussians3D_mexCode(emxArray_real_T *x0, const emxArray_real_T *image,
  const emxArray_real_T *b_index, const real_T psfSigma[2], emxArray_real_T *F,
  emxArray_real_T *J)
{
  real_T bgAmp;
  int32_T loop_ub;
  emxArray_real_T *b_x0;
  int32_T i0;
  emxArray_real_T *c_x0;
  real_T numPSF;
  int32_T ixstart;
  emxArray_real_T *d_x0;
  int32_T i;
  emxArray_real_T *psfAmp;
  real_T mtmp;
  int32_T idx;
  int32_T exitg13;
  int32_T exitg12;
  real_T b_mtmp;
  int32_T exitg11;
  int32_T exitg10;
  real_T c_mtmp;
  int32_T exitg9;
  int32_T exitg8;
  real_T d_mtmp;
  int32_T exitg7;
  int32_T exitg6;
  real_T e_mtmp;
  int32_T exitg5;
  int32_T exitg4;
  real_T f_mtmp;
  int32_T exitg3;
  int32_T exitg2;
  emxArray_real_T *psfIntegX;
  emxArray_real_T *temp;
  emxArray_real_T *y;
  emxArray_real_T *b_y;
  int32_T n;
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
  emxArray_real_T *r0;
  emxArray_real_T *r1;
  emxArray_real_T *r2;
  emxArray_boolean_T *x;
  emxArray_int32_T *ii;
  boolean_T exitg1;
  boolean_T guard1 = FALSE;
  emxArray_int32_T *b_ii;
  emxArray_real_T *b_F;
  uint32_T h_y;
  emxArray_real_T *b_psfAmp;
  emxArray_real_T *r3;
  emxArray_real_T *c_psfAmp;
  emxArray_real_T *d_psfAmp;
  emxArray_real_T *b_J;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);

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
    loop_ub = -1;
  } else {
    loop_ub = x0->size[0] - 2;
  }

  emxInit_real_T(&b_x0, 1, TRUE);
  i0 = b_x0->size[0];
  b_x0->size[0] = loop_ub + 1;
  emxEnsureCapacity((emxArray__common *)b_x0, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 <= loop_ub; i0++) {
    b_x0->data[i0] = x0->data[i0];
  }

  i0 = x0->size[0];
  x0->size[0] = b_x0->size[0];
  emxEnsureCapacity((emxArray__common *)x0, i0, (int32_T)sizeof(real_T));
  loop_ub = b_x0->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    x0->data[i0] = b_x0->data[i0];
  }

  emxFree_real_T(&b_x0);
  b_emxInit_real_T(&c_x0, 2, TRUE);

  /* get number of PSFs considered */
  /*  numPSF = length(x0)/3; */
  numPSF = (real_T)x0->size[0] / 4.0;

  /*  %reshape 3nx1 vector x0 into nx3 matrix */
  /*  x0 = reshape(x0,3,numPSF); */
  /*  x0 = x0'; */
  /* reshape 4nx1 vector x0 into nx4 matrix */
  i0 = c_x0->size[0] * c_x0->size[1];
  c_x0->size[0] = 4;
  c_x0->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)c_x0, i0, (int32_T)sizeof(real_T));
  for (ixstart = 0; ixstart + 1 <= x0->size[0]; ixstart++) {
    c_x0->data[ixstart] = x0->data[ixstart];
  }

  b_emxInit_real_T(&d_x0, 2, TRUE);
  i0 = d_x0->size[0] * d_x0->size[1];
  d_x0->size[0] = c_x0->size[1];
  d_x0->size[1] = 4;
  emxEnsureCapacity((emxArray__common *)d_x0, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < 4; i0++) {
    loop_ub = c_x0->size[1];
    for (i = 0; i < loop_ub; i++) {
      d_x0->data[i + d_x0->size[0] * i0] = c_x0->data[i0 + c_x0->size[0] * i];
    }
  }

  emxFree_real_T(&c_x0);
  emxInit_real_T(&psfAmp, 1, TRUE);

  /* extract PSF center positions and amplitudes */
  /*  psfPos = x0(:,1:2); */
  /*  psfAmp = x0(:,3); */
  loop_ub = d_x0->size[0];
  i0 = psfAmp->size[0];
  psfAmp->size[0] = loop_ub;
  emxEnsureCapacity((emxArray__common *)psfAmp, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < loop_ub; i0++) {
    psfAmp->data[i0] = d_x0->data[i0 + d_x0->size[0] * 3];
  }

  /* find minimum and maximum pixel indices */
  ixstart = 1;
  mtmp = b_index->data[0];
  i0 = b_index->size[0];
  if (i0 > 1) {
    if (muDoubleScalarIsNaN(mtmp)) {
      idx = 2;
      do {
        exitg13 = 0;
        i0 = b_index->size[0];
        if (idx <= i0) {
          ixstart = idx;
          if (!muDoubleScalarIsNaN(b_index->data[idx - 1])) {
            mtmp = b_index->data[idx - 1];
            exitg13 = 1;
          } else {
            idx++;
          }
        } else {
          exitg13 = 1;
        }
      } while (exitg13 == 0);
    }

    i0 = b_index->size[0];
    if (ixstart < i0) {
      do {
        exitg12 = 0;
        i0 = b_index->size[0];
        if (ixstart + 1 <= i0) {
          if (b_index->data[ixstart] < mtmp) {
            mtmp = b_index->data[ixstart];
          }

          ixstart++;
        } else {
          exitg12 = 1;
        }
      } while (exitg12 == 0);
    }
  }

  ixstart = 1;
  b_mtmp = b_index->data[0];
  i0 = b_index->size[0];
  if (i0 > 1) {
    if (muDoubleScalarIsNaN(b_mtmp)) {
      idx = 2;
      do {
        exitg11 = 0;
        i0 = b_index->size[0];
        if (idx <= i0) {
          ixstart = idx;
          if (!muDoubleScalarIsNaN(b_index->data[idx - 1])) {
            b_mtmp = b_index->data[idx - 1];
            exitg11 = 1;
          } else {
            idx++;
          }
        } else {
          exitg11 = 1;
        }
      } while (exitg11 == 0);
    }

    i0 = b_index->size[0];
    if (ixstart < i0) {
      do {
        exitg10 = 0;
        i0 = b_index->size[0];
        if (ixstart + 1 <= i0) {
          if (b_index->data[ixstart] > b_mtmp) {
            b_mtmp = b_index->data[ixstart];
          }

          ixstart++;
        } else {
          exitg10 = 1;
        }
      } while (exitg10 == 0);
    }
  }

  ixstart = 1;
  c_mtmp = b_index->data[b_index->size[0]];
  i0 = b_index->size[0];
  if (i0 > 1) {
    if (muDoubleScalarIsNaN(c_mtmp)) {
      idx = 2;
      do {
        exitg9 = 0;
        i0 = b_index->size[0];
        if (idx <= i0) {
          ixstart = idx;
          if (!muDoubleScalarIsNaN(b_index->data[(idx + b_index->size[0]) - 1]))
          {
            c_mtmp = b_index->data[(idx + b_index->size[0]) - 1];
            exitg9 = 1;
          } else {
            idx++;
          }
        } else {
          exitg9 = 1;
        }
      } while (exitg9 == 0);
    }

    i0 = b_index->size[0];
    if (ixstart < i0) {
      do {
        exitg8 = 0;
        i0 = b_index->size[0];
        if (ixstart + 1 <= i0) {
          if (b_index->data[ixstart + b_index->size[0]] < c_mtmp) {
            c_mtmp = b_index->data[ixstart + b_index->size[0]];
          }

          ixstart++;
        } else {
          exitg8 = 1;
        }
      } while (exitg8 == 0);
    }
  }

  ixstart = 1;
  d_mtmp = b_index->data[b_index->size[0]];
  i0 = b_index->size[0];
  if (i0 > 1) {
    if (muDoubleScalarIsNaN(d_mtmp)) {
      idx = 2;
      do {
        exitg7 = 0;
        i0 = b_index->size[0];
        if (idx <= i0) {
          ixstart = idx;
          if (!muDoubleScalarIsNaN(b_index->data[(idx + b_index->size[0]) - 1]))
          {
            d_mtmp = b_index->data[(idx + b_index->size[0]) - 1];
            exitg7 = 1;
          } else {
            idx++;
          }
        } else {
          exitg7 = 1;
        }
      } while (exitg7 == 0);
    }

    i0 = b_index->size[0];
    if (ixstart < i0) {
      do {
        exitg6 = 0;
        i0 = b_index->size[0];
        if (ixstart + 1 <= i0) {
          if (b_index->data[ixstart + b_index->size[0]] > d_mtmp) {
            d_mtmp = b_index->data[ixstart + b_index->size[0]];
          }

          ixstart++;
        } else {
          exitg6 = 1;
        }
      } while (exitg6 == 0);
    }
  }

  ixstart = 1;
  e_mtmp = b_index->data[b_index->size[0] << 1];
  i0 = b_index->size[0];
  if (i0 > 1) {
    if (muDoubleScalarIsNaN(e_mtmp)) {
      idx = 2;
      do {
        exitg5 = 0;
        i0 = b_index->size[0];
        if (idx <= i0) {
          ixstart = idx;
          if (!muDoubleScalarIsNaN(b_index->data[(idx + (b_index->size[0] << 1))
               - 1])) {
            e_mtmp = b_index->data[(idx + (b_index->size[0] << 1)) - 1];
            exitg5 = 1;
          } else {
            idx++;
          }
        } else {
          exitg5 = 1;
        }
      } while (exitg5 == 0);
    }

    i0 = b_index->size[0];
    if (ixstart < i0) {
      do {
        exitg4 = 0;
        i0 = b_index->size[0];
        if (ixstart + 1 <= i0) {
          if (b_index->data[ixstart + (b_index->size[0] << 1)] < e_mtmp) {
            e_mtmp = b_index->data[ixstart + (b_index->size[0] << 1)];
          }

          ixstart++;
        } else {
          exitg4 = 1;
        }
      } while (exitg4 == 0);
    }
  }

  ixstart = 1;
  f_mtmp = b_index->data[b_index->size[0] << 1];
  i0 = b_index->size[0];
  if (i0 > 1) {
    if (muDoubleScalarIsNaN(f_mtmp)) {
      idx = 2;
      do {
        exitg3 = 0;
        i0 = b_index->size[0];
        if (idx <= i0) {
          ixstart = idx;
          if (!muDoubleScalarIsNaN(b_index->data[(idx + (b_index->size[0] << 1))
               - 1])) {
            f_mtmp = b_index->data[(idx + (b_index->size[0] << 1)) - 1];
            exitg3 = 1;
          } else {
            idx++;
          }
        } else {
          exitg3 = 1;
        }
      } while (exitg3 == 0);
    }

    i0 = b_index->size[0];
    if (ixstart < i0) {
      do {
        exitg2 = 0;
        i0 = b_index->size[0];
        if (ixstart + 1 <= i0) {
          if (b_index->data[ixstart + (b_index->size[0] << 1)] > f_mtmp) {
            f_mtmp = b_index->data[ixstart + (b_index->size[0] << 1)];
          }

          ixstart++;
        } else {
          exitg2 = 1;
        }
      } while (exitg2 == 0);
    }
  }

  b_emxInit_real_T(&psfIntegX, 2, TRUE);

  /* determine the contribution of each PSF (assuming amplitude 1) to a */
  /* pixel based on its x-coordinate (needed to calculate F & J) */
  i0 = psfIntegX->size[0] * psfIntegX->size[1];
  psfIntegX->size[0] = (int32_T)((b_mtmp - mtmp) + 1.0);
  psfIntegX->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfIntegX, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((b_mtmp - mtmp) + 1.0) * (int32_T)numPSF;
  for (i0 = 0; i0 < loop_ub; i0++) {
    psfIntegX->data[i0] = 0.0;
  }

  i = 0;
  c_emxInit_real_T(&temp, 3, TRUE);
  b_emxInit_real_T(&y, 2, TRUE);
  emxInit_real_T(&b_y, 1, TRUE);
  while (i <= (int32_T)numPSF - 1) {
    if (muDoubleScalarIsNaN(mtmp) || muDoubleScalarIsNaN(b_mtmp)) {
      n = 0;
      anew = rtNaN;
      apnd = b_mtmp;
    } else if (b_mtmp < mtmp) {
      n = -1;
      anew = mtmp;
      apnd = b_mtmp;
    } else if (muDoubleScalarIsInf(mtmp) || muDoubleScalarIsInf(b_mtmp)) {
      n = 0;
      anew = rtNaN;
      apnd = b_mtmp;
    } else {
      anew = mtmp;
      ndbl = muDoubleScalarFloor((b_mtmp - mtmp) + 0.5);
      apnd = mtmp + ndbl;
      cdiff = apnd - b_mtmp;
      absa = muDoubleScalarAbs(mtmp);
      absb = muDoubleScalarAbs(b_mtmp);
      if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * muDoubleScalarMax
          (absa, absb)) {
        ndbl++;
        apnd = b_mtmp;
      } else if (cdiff > 0.0) {
        apnd = mtmp + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int32_T)ndbl - 1;
      } else {
        n = -1;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n + 1;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (n + 1 > 0) {
      y->data[0] = anew;
      if (n + 1 > 1) {
        y->data[n] = apnd;
        i0 = n + (n < 0);
        if (i0 >= 0) {
          idx = (int32_T)((uint32_T)i0 >> 1);
        } else {
          idx = ~(int32_T)((uint32_T)~i0 >> 1);
        }

        for (ixstart = 1; ixstart < idx; ixstart++) {
          y->data[ixstart] = anew + (real_T)ixstart;
          y->data[n - ixstart] = apnd - (real_T)ixstart;
        }

        if (idx << 1 == n) {
          y->data[idx] = (anew + apnd) / 2.0;
        } else {
          y->data[idx] = anew + (real_T)idx;
          y->data[idx + 1] = apnd - (real_T)idx;
        }
      }
    }

    i0 = b_y->size[0];
    b_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)b_y, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_y->data[i0] = y->data[y->size[0] * i0];
    }

    GaussListND_mexCode(b_y, psfSigma[0], d_x0->data[i], temp);

    /* temp2 = squeeze(temp); */
    /* clear temp */
    loop_ub = temp->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      psfIntegX->data[i0 + psfIntegX->size[0] * i] = temp->data[i0];
    }

    /* clear temp2 */
    i++;
  }

  emxFree_real_T(&b_y);
  b_emxInit_real_T(&psfIntegY, 2, TRUE);

  /* determine the contribution of each PSF (assuming amplitude 1) to a */
  /* pixel based on its y-coordinate (needed to calculate F & J) */
  i0 = psfIntegY->size[0] * psfIntegY->size[1];
  psfIntegY->size[0] = (int32_T)((d_mtmp - c_mtmp) + 1.0);
  psfIntegY->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfIntegY, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((d_mtmp - c_mtmp) + 1.0) * (int32_T)numPSF;
  for (i0 = 0; i0 < loop_ub; i0++) {
    psfIntegY->data[i0] = 0.0;
  }

  i = 0;
  emxInit_real_T(&c_y, 1, TRUE);
  while (i <= (int32_T)numPSF - 1) {
    if (muDoubleScalarIsNaN(c_mtmp) || muDoubleScalarIsNaN(d_mtmp)) {
      n = 0;
      anew = rtNaN;
      apnd = d_mtmp;
    } else if (d_mtmp < c_mtmp) {
      n = -1;
      anew = c_mtmp;
      apnd = d_mtmp;
    } else if (muDoubleScalarIsInf(c_mtmp) || muDoubleScalarIsInf(d_mtmp)) {
      n = 0;
      anew = rtNaN;
      apnd = d_mtmp;
    } else {
      anew = c_mtmp;
      ndbl = muDoubleScalarFloor((d_mtmp - c_mtmp) + 0.5);
      apnd = c_mtmp + ndbl;
      cdiff = apnd - d_mtmp;
      absa = muDoubleScalarAbs(c_mtmp);
      absb = muDoubleScalarAbs(d_mtmp);
      if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * muDoubleScalarMax
          (absa, absb)) {
        ndbl++;
        apnd = d_mtmp;
      } else if (cdiff > 0.0) {
        apnd = c_mtmp + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int32_T)ndbl - 1;
      } else {
        n = -1;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n + 1;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (n + 1 > 0) {
      y->data[0] = anew;
      if (n + 1 > 1) {
        y->data[n] = apnd;
        i0 = n + (n < 0);
        if (i0 >= 0) {
          idx = (int32_T)((uint32_T)i0 >> 1);
        } else {
          idx = ~(int32_T)((uint32_T)~i0 >> 1);
        }

        for (ixstart = 1; ixstart < idx; ixstart++) {
          y->data[ixstart] = anew + (real_T)ixstart;
          y->data[n - ixstart] = apnd - (real_T)ixstart;
        }

        if (idx << 1 == n) {
          y->data[idx] = (anew + apnd) / 2.0;
        } else {
          y->data[idx] = anew + (real_T)idx;
          y->data[idx + 1] = apnd - (real_T)idx;
        }
      }
    }

    i0 = c_y->size[0];
    c_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)c_y, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_y->data[i0] = y->data[y->size[0] * i0];
    }

    GaussListND_mexCode(c_y, psfSigma[0], d_x0->data[i + d_x0->size[0]], temp);

    /* temp2 = squeeze(temp); */
    /* clear temp */
    loop_ub = temp->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      psfIntegY->data[i0 + psfIntegY->size[0] * i] = temp->data[i0];
    }

    /* clear temp2 */
    i++;
  }

  emxFree_real_T(&c_y);
  b_emxInit_real_T(&psfIntegZ, 2, TRUE);

  /* determine the contribution of each PSF (assuming amplitude 1) to a */
  /* pixel based on its z-coordinate (needed to calculate F & J) */
  i0 = psfIntegZ->size[0] * psfIntegZ->size[1];
  psfIntegZ->size[0] = (int32_T)((f_mtmp - e_mtmp) + 1.0);
  psfIntegZ->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfIntegZ, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((f_mtmp - e_mtmp) + 1.0) * (int32_T)numPSF;
  for (i0 = 0; i0 < loop_ub; i0++) {
    psfIntegZ->data[i0] = 0.0;
  }

  i = 0;
  emxInit_real_T(&d_y, 1, TRUE);
  while (i <= (int32_T)numPSF - 1) {
    if (muDoubleScalarIsNaN(e_mtmp) || muDoubleScalarIsNaN(f_mtmp)) {
      n = 0;
      anew = rtNaN;
      apnd = f_mtmp;
    } else if (f_mtmp < e_mtmp) {
      n = -1;
      anew = e_mtmp;
      apnd = f_mtmp;
    } else if (muDoubleScalarIsInf(e_mtmp) || muDoubleScalarIsInf(f_mtmp)) {
      n = 0;
      anew = rtNaN;
      apnd = f_mtmp;
    } else {
      anew = e_mtmp;
      ndbl = muDoubleScalarFloor((f_mtmp - e_mtmp) + 0.5);
      apnd = e_mtmp + ndbl;
      cdiff = apnd - f_mtmp;
      absa = muDoubleScalarAbs(e_mtmp);
      absb = muDoubleScalarAbs(f_mtmp);
      if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * muDoubleScalarMax
          (absa, absb)) {
        ndbl++;
        apnd = f_mtmp;
      } else if (cdiff > 0.0) {
        apnd = e_mtmp + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int32_T)ndbl - 1;
      } else {
        n = -1;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n + 1;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (n + 1 > 0) {
      y->data[0] = anew;
      if (n + 1 > 1) {
        y->data[n] = apnd;
        i0 = n + (n < 0);
        if (i0 >= 0) {
          idx = (int32_T)((uint32_T)i0 >> 1);
        } else {
          idx = ~(int32_T)((uint32_T)~i0 >> 1);
        }

        for (ixstart = 1; ixstart < idx; ixstart++) {
          y->data[ixstart] = anew + (real_T)ixstart;
          y->data[n - ixstart] = apnd - (real_T)ixstart;
        }

        if (idx << 1 == n) {
          y->data[idx] = (anew + apnd) / 2.0;
        } else {
          y->data[idx] = anew + (real_T)idx;
          y->data[idx + 1] = apnd - (real_T)idx;
        }
      }
    }

    i0 = d_y->size[0];
    d_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)d_y, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      d_y->data[i0] = y->data[y->size[0] * i0];
    }

    GaussListND_mexCode(d_y, psfSigma[1], d_x0->data[i + (d_x0->size[0] << 1)],
                        temp);

    /* temp2 = squeeze(temp); */
    /* clear temp */
    loop_ub = temp->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      psfIntegZ->data[i0 + psfIntegZ->size[0] * i] = temp->data[i0];
    }

    /* clear temp2 */
    i++;
  }

  emxFree_real_T(&d_y);
  emxFree_real_T(&temp);
  b_emxInit_real_T(&psfValueX, 2, TRUE);

  /* calculate the value of each PSF (assuming amplitude 1) at the */
  /* x-coordinates of the corners of all pixels (needed to calculate J) */
  i0 = psfValueX->size[0] * psfValueX->size[1];
  psfValueX->size[0] = (int32_T)((b_mtmp - mtmp) + 2.0);
  psfValueX->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfValueX, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((b_mtmp - mtmp) + 2.0) * (int32_T)numPSF;
  for (i0 = 0; i0 < loop_ub; i0++) {
    psfValueX->data[i0] = 0.0;
  }

  i = 0;
  emxInit_real_T(&indxPixel, 1, TRUE);
  emxInit_real_T(&b_indxPixel, 1, TRUE);
  emxInit_real_T(&c_indxPixel, 1, TRUE);
  emxInit_real_T(&e_y, 1, TRUE);
  while (i <= (int32_T)numPSF - 1) {
    if (muDoubleScalarIsNaN(mtmp - 0.5) || muDoubleScalarIsNaN(b_mtmp + 0.5)) {
      n = 0;
      anew = rtNaN;
      apnd = b_mtmp + 0.5;
    } else if (b_mtmp + 0.5 < mtmp - 0.5) {
      n = -1;
      anew = mtmp - 0.5;
      apnd = b_mtmp + 0.5;
    } else if (muDoubleScalarIsInf(mtmp - 0.5) || muDoubleScalarIsInf(b_mtmp +
                0.5)) {
      n = 0;
      anew = rtNaN;
      apnd = b_mtmp + 0.5;
    } else {
      anew = mtmp - 0.5;
      ndbl = muDoubleScalarFloor(((b_mtmp + 0.5) - (mtmp - 0.5)) + 0.5);
      apnd = (mtmp - 0.5) + ndbl;
      cdiff = apnd - (b_mtmp + 0.5);
      absa = muDoubleScalarAbs(mtmp - 0.5);
      absb = muDoubleScalarAbs(b_mtmp + 0.5);
      if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * muDoubleScalarMax
          (absa, absb)) {
        ndbl++;
        apnd = b_mtmp + 0.5;
      } else if (cdiff > 0.0) {
        apnd = (mtmp - 0.5) + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int32_T)ndbl - 1;
      } else {
        n = -1;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n + 1;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (n + 1 > 0) {
      y->data[0] = anew;
      if (n + 1 > 1) {
        y->data[n] = apnd;
        i0 = n + (n < 0);
        if (i0 >= 0) {
          idx = (int32_T)((uint32_T)i0 >> 1);
        } else {
          idx = ~(int32_T)((uint32_T)~i0 >> 1);
        }

        for (ixstart = 1; ixstart < idx; ixstart++) {
          y->data[ixstart] = anew + (real_T)ixstart;
          y->data[n - ixstart] = apnd - (real_T)ixstart;
        }

        if (idx << 1 == n) {
          y->data[idx] = (anew + apnd) / 2.0;
        } else {
          y->data[idx] = anew + (real_T)idx;
          y->data[idx + 1] = apnd - (real_T)idx;
        }
      }
    }

    anew = d_x0->data[i];
    i0 = e_y->size[0];
    e_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)e_y, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      e_y->data[i0] = y->data[y->size[0] * i0] - anew;
    }

    power(e_y, indxPixel);
    i0 = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)indxPixel, i0, (int32_T)sizeof(real_T));
    loop_ub = indxPixel->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      indxPixel->data[i0] = -indxPixel->data[i0];
    }

    i0 = c_indxPixel->size[0];
    c_indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)c_indxPixel, i0, (int32_T)sizeof
                      (real_T));
    loop_ub = indxPixel->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      c_indxPixel->data[i0] = indxPixel->data[i0];
    }

    rdivide(c_indxPixel, 2.0, indxPixel);
    i0 = b_indxPixel->size[0];
    b_indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)b_indxPixel, i0, (int32_T)sizeof
                      (real_T));
    loop_ub = indxPixel->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_indxPixel->data[i0] = indxPixel->data[i0];
    }

    rdivide(b_indxPixel, psfSigma[0] * psfSigma[0], indxPixel);
    b_exp(indxPixel);
    loop_ub = indxPixel->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      psfValueX->data[i0 + psfValueX->size[0] * i] = indxPixel->data[i0];
    }

    i++;
  }

  emxFree_real_T(&e_y);
  emxFree_real_T(&c_indxPixel);
  emxFree_real_T(&b_indxPixel);
  b_emxInit_real_T(&psfValueY, 2, TRUE);

  /* calculate the value of each PSF (assuming amplitude 1) at the */
  /* y-coordinates of the corners of all pixels (needed to calculate J) */
  i0 = psfValueY->size[0] * psfValueY->size[1];
  psfValueY->size[0] = (int32_T)((d_mtmp - c_mtmp) + 2.0);
  psfValueY->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfValueY, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((d_mtmp - c_mtmp) + 2.0) * (int32_T)numPSF;
  for (i0 = 0; i0 < loop_ub; i0++) {
    psfValueY->data[i0] = 0.0;
  }

  i = 0;
  emxInit_real_T(&d_indxPixel, 1, TRUE);
  emxInit_real_T(&e_indxPixel, 1, TRUE);
  emxInit_real_T(&f_y, 1, TRUE);
  while (i <= (int32_T)numPSF - 1) {
    if (muDoubleScalarIsNaN(c_mtmp - 0.5) || muDoubleScalarIsNaN(d_mtmp + 0.5))
    {
      n = 0;
      anew = rtNaN;
      apnd = d_mtmp + 0.5;
    } else if (d_mtmp + 0.5 < c_mtmp - 0.5) {
      n = -1;
      anew = c_mtmp - 0.5;
      apnd = d_mtmp + 0.5;
    } else if (muDoubleScalarIsInf(c_mtmp - 0.5) || muDoubleScalarIsInf(d_mtmp +
                0.5)) {
      n = 0;
      anew = rtNaN;
      apnd = d_mtmp + 0.5;
    } else {
      anew = c_mtmp - 0.5;
      ndbl = muDoubleScalarFloor(((d_mtmp + 0.5) - (c_mtmp - 0.5)) + 0.5);
      apnd = (c_mtmp - 0.5) + ndbl;
      cdiff = apnd - (d_mtmp + 0.5);
      absa = muDoubleScalarAbs(c_mtmp - 0.5);
      absb = muDoubleScalarAbs(d_mtmp + 0.5);
      if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * muDoubleScalarMax
          (absa, absb)) {
        ndbl++;
        apnd = d_mtmp + 0.5;
      } else if (cdiff > 0.0) {
        apnd = (c_mtmp - 0.5) + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int32_T)ndbl - 1;
      } else {
        n = -1;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n + 1;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (n + 1 > 0) {
      y->data[0] = anew;
      if (n + 1 > 1) {
        y->data[n] = apnd;
        i0 = n + (n < 0);
        if (i0 >= 0) {
          idx = (int32_T)((uint32_T)i0 >> 1);
        } else {
          idx = ~(int32_T)((uint32_T)~i0 >> 1);
        }

        for (ixstart = 1; ixstart < idx; ixstart++) {
          y->data[ixstart] = anew + (real_T)ixstart;
          y->data[n - ixstart] = apnd - (real_T)ixstart;
        }

        if (idx << 1 == n) {
          y->data[idx] = (anew + apnd) / 2.0;
        } else {
          y->data[idx] = anew + (real_T)idx;
          y->data[idx + 1] = apnd - (real_T)idx;
        }
      }
    }

    anew = d_x0->data[i + d_x0->size[0]];
    i0 = f_y->size[0];
    f_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)f_y, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      f_y->data[i0] = y->data[y->size[0] * i0] - anew;
    }

    power(f_y, indxPixel);
    i0 = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)indxPixel, i0, (int32_T)sizeof(real_T));
    loop_ub = indxPixel->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      indxPixel->data[i0] = -indxPixel->data[i0];
    }

    i0 = e_indxPixel->size[0];
    e_indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)e_indxPixel, i0, (int32_T)sizeof
                      (real_T));
    loop_ub = indxPixel->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      e_indxPixel->data[i0] = indxPixel->data[i0];
    }

    rdivide(e_indxPixel, 2.0, indxPixel);
    i0 = d_indxPixel->size[0];
    d_indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)d_indxPixel, i0, (int32_T)sizeof
                      (real_T));
    loop_ub = indxPixel->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      d_indxPixel->data[i0] = indxPixel->data[i0];
    }

    rdivide(d_indxPixel, psfSigma[0] * psfSigma[0], indxPixel);
    b_exp(indxPixel);
    loop_ub = indxPixel->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      psfValueY->data[i0 + psfValueY->size[0] * i] = indxPixel->data[i0];
    }

    i++;
  }

  emxFree_real_T(&f_y);
  emxFree_real_T(&e_indxPixel);
  emxFree_real_T(&d_indxPixel);
  b_emxInit_real_T(&psfValueZ, 2, TRUE);

  /* calculate the value of each PSF (assuming amplitude 1) at the */
  /* z-coordinates of the corners of all pixels (needed to calculate J) */
  i0 = psfValueZ->size[0] * psfValueZ->size[1];
  psfValueZ->size[0] = (int32_T)((f_mtmp - e_mtmp) + 2.0);
  psfValueZ->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfValueZ, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((f_mtmp - e_mtmp) + 2.0) * (int32_T)numPSF;
  for (i0 = 0; i0 < loop_ub; i0++) {
    psfValueZ->data[i0] = 0.0;
  }

  i = 0;
  emxInit_real_T(&f_indxPixel, 1, TRUE);
  emxInit_real_T(&g_indxPixel, 1, TRUE);
  emxInit_real_T(&g_y, 1, TRUE);
  while (i <= (int32_T)numPSF - 1) {
    if (muDoubleScalarIsNaN(e_mtmp - 0.5) || muDoubleScalarIsNaN(f_mtmp + 0.5))
    {
      n = 0;
      anew = rtNaN;
      apnd = f_mtmp + 0.5;
    } else if (f_mtmp + 0.5 < e_mtmp - 0.5) {
      n = -1;
      anew = e_mtmp - 0.5;
      apnd = f_mtmp + 0.5;
    } else if (muDoubleScalarIsInf(e_mtmp - 0.5) || muDoubleScalarIsInf(f_mtmp +
                0.5)) {
      n = 0;
      anew = rtNaN;
      apnd = f_mtmp + 0.5;
    } else {
      anew = e_mtmp - 0.5;
      ndbl = muDoubleScalarFloor(((f_mtmp + 0.5) - (e_mtmp - 0.5)) + 0.5);
      apnd = (e_mtmp - 0.5) + ndbl;
      cdiff = apnd - (f_mtmp + 0.5);
      absa = muDoubleScalarAbs(e_mtmp - 0.5);
      absb = muDoubleScalarAbs(f_mtmp + 0.5);
      if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * muDoubleScalarMax
          (absa, absb)) {
        ndbl++;
        apnd = f_mtmp + 0.5;
      } else if (cdiff > 0.0) {
        apnd = (e_mtmp - 0.5) + (ndbl - 1.0);
      } else {
        ndbl++;
      }

      if (ndbl >= 0.0) {
        n = (int32_T)ndbl - 1;
      } else {
        n = -1;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n + 1;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (n + 1 > 0) {
      y->data[0] = anew;
      if (n + 1 > 1) {
        y->data[n] = apnd;
        i0 = n + (n < 0);
        if (i0 >= 0) {
          idx = (int32_T)((uint32_T)i0 >> 1);
        } else {
          idx = ~(int32_T)((uint32_T)~i0 >> 1);
        }

        for (ixstart = 1; ixstart < idx; ixstart++) {
          y->data[ixstart] = anew + (real_T)ixstart;
          y->data[n - ixstart] = apnd - (real_T)ixstart;
        }

        if (idx << 1 == n) {
          y->data[idx] = (anew + apnd) / 2.0;
        } else {
          y->data[idx] = anew + (real_T)idx;
          y->data[idx + 1] = apnd - (real_T)idx;
        }
      }
    }

    anew = d_x0->data[i + (d_x0->size[0] << 1)];
    i0 = g_y->size[0];
    g_y->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)g_y, i0, (int32_T)sizeof(real_T));
    loop_ub = y->size[1];
    for (i0 = 0; i0 < loop_ub; i0++) {
      g_y->data[i0] = y->data[y->size[0] * i0] - anew;
    }

    power(g_y, indxPixel);
    i0 = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)indxPixel, i0, (int32_T)sizeof(real_T));
    loop_ub = indxPixel->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      indxPixel->data[i0] = -indxPixel->data[i0];
    }

    i0 = g_indxPixel->size[0];
    g_indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)g_indxPixel, i0, (int32_T)sizeof
                      (real_T));
    loop_ub = indxPixel->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      g_indxPixel->data[i0] = indxPixel->data[i0];
    }

    rdivide(g_indxPixel, 2.0, indxPixel);
    i0 = f_indxPixel->size[0];
    f_indxPixel->size[0] = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)f_indxPixel, i0, (int32_T)sizeof
                      (real_T));
    loop_ub = indxPixel->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      f_indxPixel->data[i0] = indxPixel->data[i0];
    }

    rdivide(f_indxPixel, psfSigma[1] * psfSigma[1], indxPixel);
    b_exp(indxPixel);
    loop_ub = indxPixel->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      psfValueZ->data[i0 + psfValueZ->size[0] * i] = indxPixel->data[i0];
    }

    i++;
  }

  emxFree_real_T(&g_y);
  emxFree_real_T(&g_indxPixel);
  emxFree_real_T(&f_indxPixel);
  emxFree_real_T(&d_x0);
  emxInit_real_T(&relIndxX, 1, TRUE);

  /* get number of pixels in image */
  /*  %get xy-indices relative to minimum */
  /*  relIndxX = index(:,1) - minIndxX + 1; */
  /*  relIndxY = index(:,2) - minIndxY + 1; */
  /* get xy-indices relative to minimum */
  loop_ub = b_index->size[0];
  i0 = relIndxX->size[0];
  relIndxX->size[0] = loop_ub;
  emxEnsureCapacity((emxArray__common *)relIndxX, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < loop_ub; i0++) {
    relIndxX->data[i0] = (b_index->data[i0] - mtmp) + 1.0;
  }

  emxInit_real_T(&relIndxY, 1, TRUE);
  loop_ub = b_index->size[0];
  i0 = relIndxY->size[0];
  relIndxY->size[0] = loop_ub;
  emxEnsureCapacity((emxArray__common *)relIndxY, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < loop_ub; i0++) {
    relIndxY->data[i0] = (b_index->data[i0 + b_index->size[0]] - c_mtmp) + 1.0;
  }

  emxInit_real_T(&relIndxZ, 1, TRUE);
  loop_ub = b_index->size[0];
  i0 = relIndxZ->size[0];
  relIndxZ->size[0] = loop_ub;
  emxEnsureCapacity((emxArray__common *)relIndxZ, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < loop_ub; i0++) {
    relIndxZ->data[i0] = (b_index->data[i0 + (b_index->size[0] << 1)] - e_mtmp)
      + 1.0;
  }

  b_emxInit_real_T(&r0, 2, TRUE);
  b_emxInit_real_T(&r1, 2, TRUE);

  /*  %calculate the value of F at all pixels */
  /*  F = (sum(repmat(psfAmp,1,numPixel).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)',1))' ... */
  /*      + repmat(bgAmp,numPixel,1) - image; */
  /* calculate the value of F at all pixels */
  b_repmat(psfAmp, image->size[0], r0);
  i0 = r1->size[0] * r1->size[1];
  r1->size[0] = r0->size[0];
  r1->size[1] = r0->size[1];
  emxEnsureCapacity((emxArray__common *)r1, i0, (int32_T)sizeof(real_T));
  loop_ub = r0->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    idx = r0->size[0];
    for (i = 0; i < idx; i++) {
      r1->data[i + r1->size[0] * i0] = r0->data[i + r0->size[0] * i0] *
        psfIntegX->data[((int32_T)relIndxX->data[i0] + psfIntegX->size[0] * i) -
        1] * psfIntegY->data[((int32_T)relIndxY->data[i0] + psfIntegY->size[0] *
        i) - 1] * psfIntegZ->data[((int32_T)relIndxZ->data[i0] + psfIntegZ->
        size[0] * i) - 1];
    }
  }

  sum(r1, y);
  i0 = indxPixel->size[0];
  indxPixel->size[0] = y->size[1];
  emxEnsureCapacity((emxArray__common *)indxPixel, i0, (int32_T)sizeof(real_T));
  loop_ub = y->size[1];
  emxFree_real_T(&r1);
  for (i0 = 0; i0 < loop_ub; i0++) {
    indxPixel->data[i0] = y->data[y->size[0] * i0];
  }

  emxFree_real_T(&y);
  emxInit_real_T(&r2, 1, TRUE);
  c_repmat(bgAmp, image->size[0], r2);
  ixstart = indxPixel->size[0];
  i0 = F->size[0] * F->size[1];
  F->size[0] = ixstart;
  emxEnsureCapacity((emxArray__common *)F, i0, (int32_T)sizeof(real_T));
  i0 = F->size[0] * F->size[1];
  F->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)F, i0, (int32_T)sizeof(real_T));
  loop_ub = indxPixel->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    F->data[i0] = (indxPixel->data[i0] + r2->data[i0]) - image->data[i0];
  }

  emxFree_real_T(&r2);
  emxInit_boolean_T(&x, 1, TRUE);

  /* remove pixels with NaN (which means they are out of the cropped image */
  /* area) */
  i0 = x->size[0];
  x->size[0] = image->size[0];
  emxEnsureCapacity((emxArray__common *)x, i0, (int32_T)sizeof(boolean_T));
  loop_ub = image->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    x->data[i0] = muDoubleScalarIsNaN(image->data[i0]);
  }

  i0 = x->size[0];
  emxEnsureCapacity((emxArray__common *)x, i0, (int32_T)sizeof(boolean_T));
  loop_ub = x->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    x->data[i0] = !x->data[i0];
  }

  emxInit_int32_T(&ii, 1, TRUE);
  idx = 0;
  i0 = ii->size[0];
  ii->size[0] = x->size[0];
  emxEnsureCapacity((emxArray__common *)ii, i0, (int32_T)sizeof(int32_T));
  ixstart = 1;
  exitg1 = FALSE;
  while ((exitg1 == FALSE) && (ixstart <= x->size[0])) {
    guard1 = FALSE;
    if (x->data[ixstart - 1]) {
      idx++;
      ii->data[idx - 1] = ixstart;
      if (idx >= x->size[0]) {
        exitg1 = TRUE;
      } else {
        guard1 = TRUE;
      }
    } else {
      guard1 = TRUE;
    }

    if (guard1 == TRUE) {
      ixstart++;
    }
  }

  if (x->size[0] == 1) {
    if (idx == 0) {
      i0 = ii->size[0];
      ii->size[0] = 0;
      emxEnsureCapacity((emxArray__common *)ii, i0, (int32_T)sizeof(int32_T));
    }
  } else {
    if (1 > idx) {
      loop_ub = 0;
    } else {
      loop_ub = idx;
    }

    emxInit_int32_T(&b_ii, 1, TRUE);
    i0 = b_ii->size[0];
    b_ii->size[0] = loop_ub;
    emxEnsureCapacity((emxArray__common *)b_ii, i0, (int32_T)sizeof(int32_T));
    for (i0 = 0; i0 < loop_ub; i0++) {
      b_ii->data[i0] = ii->data[i0];
    }

    i0 = ii->size[0];
    ii->size[0] = b_ii->size[0];
    emxEnsureCapacity((emxArray__common *)ii, i0, (int32_T)sizeof(int32_T));
    loop_ub = b_ii->size[0];
    for (i0 = 0; i0 < loop_ub; i0++) {
      ii->data[i0] = b_ii->data[i0];
    }

    emxFree_int32_T(&b_ii);
  }

  emxFree_boolean_T(&x);
  i0 = indxPixel->size[0];
  indxPixel->size[0] = ii->size[0];
  emxEnsureCapacity((emxArray__common *)indxPixel, i0, (int32_T)sizeof(real_T));
  loop_ub = ii->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    indxPixel->data[i0] = ii->data[i0];
  }

  emxFree_int32_T(&ii);
  emxInit_real_T(&b_F, 1, TRUE);
  i0 = b_F->size[0];
  b_F->size[0] = indxPixel->size[0];
  emxEnsureCapacity((emxArray__common *)b_F, i0, (int32_T)sizeof(real_T));
  loop_ub = indxPixel->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_F->data[i0] = F->data[(int32_T)indxPixel->data[i0] - 1];
  }

  ixstart = indxPixel->size[0];
  i0 = F->size[0] * F->size[1];
  F->size[0] = ixstart;
  F->size[1] = 1;
  emxEnsureCapacity((emxArray__common *)F, i0, (int32_T)sizeof(real_T));
  i0 = 0;
  while (i0 <= 0) {
    for (i0 = 0; i0 < ixstart; i0++) {
      F->data[i0] = b_F->data[i0];
    }

    i0 = 1;
  }

  emxFree_real_T(&b_F);

  /* calculate the derivative at all pixels */
  /*      J = ones(numPixel,3*numPSF+1); %(last column for background amplitude) */
  /*      J(:,1:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueX(relIndxX,:)-... */
  /*          psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:); %w.r.t. x */
  /*      J(:,2:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueY(relIndxY,:)-... */
  /*          psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:); %w.r.t. y */
  /*      J(:,3:3:3*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); %w.r.t. amp */
  h_y = (uint32_T)numPSF << 2;
  ixstart = image->size[0];
  i0 = J->size[0] * J->size[1];
  J->size[0] = ixstart;
  emxEnsureCapacity((emxArray__common *)J, i0, (int32_T)sizeof(real_T));
  i0 = J->size[0] * J->size[1];
  J->size[1] = (int32_T)((real_T)h_y + 1.0);
  emxEnsureCapacity((emxArray__common *)J, i0, (int32_T)sizeof(real_T));
  loop_ub = image->size[0] * (int32_T)((real_T)h_y + 1.0);
  for (i0 = 0; i0 < loop_ub; i0++) {
    J->data[i0] = 1.0;
  }

  /* (last column for background amplitude) */
  loop_ub = psfValueX->size[1];
  i0 = r0->size[0] * r0->size[1];
  r0->size[0] = relIndxX->size[0];
  r0->size[1] = loop_ub;
  emxEnsureCapacity((emxArray__common *)r0, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < loop_ub; i0++) {
    idx = relIndxX->size[0];
    for (i = 0; i < idx; i++) {
      r0->data[i + r0->size[0] * i0] = psfValueX->data[((int32_T)(relIndxX->
        data[i] + 1.0) + psfValueX->size[0] * i0) - 1];
    }
  }

  b_emxInit_real_T(&b_psfAmp, 2, TRUE);
  i0 = b_psfAmp->size[0] * b_psfAmp->size[1];
  b_psfAmp->size[0] = 1;
  b_psfAmp->size[1] = psfAmp->size[0];
  emxEnsureCapacity((emxArray__common *)b_psfAmp, i0, (int32_T)sizeof(real_T));
  loop_ub = psfAmp->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    b_psfAmp->data[b_psfAmp->size[0] * i0] = psfAmp->data[i0];
  }

  b_emxInit_real_T(&r3, 2, TRUE);
  d_repmat(b_psfAmp, image->size[0], r3);
  emxFree_real_T(&b_psfAmp);
  if (1U > ((uint32_T)numPSF << 2)) {
    i0 = 1;
  } else {
    i0 = 4;
  }

  loop_ub = r3->size[1];
  for (i = 0; i < loop_ub; i++) {
    idx = r3->size[0];
    for (ixstart = 0; ixstart < idx; ixstart++) {
      J->data[ixstart + J->size[0] * (i0 * i)] = r3->data[ixstart + r3->size[0] *
        i] * (psfValueX->data[((int32_T)relIndxX->data[ixstart] +
               psfValueX->size[0] * i) - 1] - r0->data[ixstart + r0->size[0] * i])
        * psfIntegY->data[((int32_T)relIndxY->data[ixstart] + psfIntegY->size[0]
                           * i) - 1] * psfIntegZ->data[((int32_T)relIndxZ->
        data[ixstart] + psfIntegZ->size[0] * i) - 1];
    }
  }

  emxFree_real_T(&psfValueX);

  /* w.r.t. x */
  loop_ub = psfValueY->size[1];
  i0 = r0->size[0] * r0->size[1];
  r0->size[0] = relIndxY->size[0];
  r0->size[1] = loop_ub;
  emxEnsureCapacity((emxArray__common *)r0, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < loop_ub; i0++) {
    idx = relIndxY->size[0];
    for (i = 0; i < idx; i++) {
      r0->data[i + r0->size[0] * i0] = psfValueY->data[((int32_T)(relIndxY->
        data[i] + 1.0) + psfValueY->size[0] * i0) - 1];
    }
  }

  b_emxInit_real_T(&c_psfAmp, 2, TRUE);
  i0 = c_psfAmp->size[0] * c_psfAmp->size[1];
  c_psfAmp->size[0] = 1;
  c_psfAmp->size[1] = psfAmp->size[0];
  emxEnsureCapacity((emxArray__common *)c_psfAmp, i0, (int32_T)sizeof(real_T));
  loop_ub = psfAmp->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    c_psfAmp->data[c_psfAmp->size[0] * i0] = psfAmp->data[i0];
  }

  d_repmat(c_psfAmp, image->size[0], r3);
  emxFree_real_T(&c_psfAmp);
  if (2U > ((uint32_T)numPSF << 2)) {
    i0 = 0;
    i = 1;
  } else {
    i0 = 1;
    i = 4;
  }

  loop_ub = r3->size[1];
  for (ixstart = 0; ixstart < loop_ub; ixstart++) {
    idx = r3->size[0];
    for (n = 0; n < idx; n++) {
      J->data[n + J->size[0] * (i0 + i * ixstart)] = r3->data[n + r3->size[0] *
        ixstart] * (psfValueY->data[((int32_T)relIndxY->data[n] +
        psfValueY->size[0] * ixstart) - 1] - r0->data[n + r0->size[0] * ixstart])
        * psfIntegX->data[((int32_T)relIndxX->data[n] + psfIntegX->size[0] *
                           ixstart) - 1] * psfIntegZ->data[((int32_T)
        relIndxZ->data[n] + psfIntegZ->size[0] * ixstart) - 1];
    }
  }

  emxFree_real_T(&psfValueY);

  /* w.r.t. y */
  loop_ub = psfValueZ->size[1];
  i0 = r0->size[0] * r0->size[1];
  r0->size[0] = relIndxZ->size[0];
  r0->size[1] = loop_ub;
  emxEnsureCapacity((emxArray__common *)r0, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < loop_ub; i0++) {
    idx = relIndxZ->size[0];
    for (i = 0; i < idx; i++) {
      r0->data[i + r0->size[0] * i0] = psfValueZ->data[((int32_T)(relIndxZ->
        data[i] + 1.0) + psfValueZ->size[0] * i0) - 1];
    }
  }

  b_emxInit_real_T(&d_psfAmp, 2, TRUE);
  i0 = d_psfAmp->size[0] * d_psfAmp->size[1];
  d_psfAmp->size[0] = 1;
  d_psfAmp->size[1] = psfAmp->size[0];
  emxEnsureCapacity((emxArray__common *)d_psfAmp, i0, (int32_T)sizeof(real_T));
  loop_ub = psfAmp->size[0];
  for (i0 = 0; i0 < loop_ub; i0++) {
    d_psfAmp->data[d_psfAmp->size[0] * i0] = psfAmp->data[i0];
  }

  emxFree_real_T(&psfAmp);
  d_repmat(d_psfAmp, image->size[0], r3);
  emxFree_real_T(&d_psfAmp);
  if (3U > ((uint32_T)numPSF << 2)) {
    i0 = 0;
    i = 1;
  } else {
    i0 = 2;
    i = 4;
  }

  loop_ub = r3->size[1];
  for (ixstart = 0; ixstart < loop_ub; ixstart++) {
    idx = r3->size[0];
    for (n = 0; n < idx; n++) {
      J->data[n + J->size[0] * (i0 + i * ixstart)] = r3->data[n + r3->size[0] *
        ixstart] * (psfValueZ->data[((int32_T)relIndxZ->data[n] +
        psfValueZ->size[0] * ixstart) - 1] - r0->data[n + r0->size[0] * ixstart])
        * psfIntegX->data[((int32_T)relIndxX->data[n] + psfIntegX->size[0] *
                           ixstart) - 1] * psfIntegY->data[((int32_T)
        relIndxY->data[n] + psfIntegY->size[0] * ixstart) - 1];
    }
  }

  emxFree_real_T(&r3);
  emxFree_real_T(&r0);
  emxFree_real_T(&psfValueZ);

  /* w.r.t. z */
  if (4U > ((uint32_T)numPSF << 2)) {
    i0 = 0;
    i = 1;
  } else {
    i0 = 3;
    i = 4;
  }

  loop_ub = psfIntegX->size[1] - 1;
  for (ixstart = 0; ixstart <= loop_ub; ixstart++) {
    idx = relIndxX->size[0];
    for (n = 0; n < idx; n++) {
      J->data[n + J->size[0] * (i0 + i * ixstart)] = psfIntegX->data[((int32_T)
        relIndxX->data[n] + psfIntegX->size[0] * ixstart) - 1] * psfIntegY->
        data[((int32_T)relIndxY->data[n] + psfIntegY->size[0] * ixstart) - 1] *
        psfIntegZ->data[((int32_T)relIndxZ->data[n] + psfIntegZ->size[0] *
                         ixstart) - 1];
    }
  }

  emxFree_real_T(&relIndxZ);
  emxFree_real_T(&relIndxY);
  emxFree_real_T(&relIndxX);
  emxFree_real_T(&psfIntegZ);
  emxFree_real_T(&psfIntegY);
  emxFree_real_T(&psfIntegX);
  b_emxInit_real_T(&b_J, 2, TRUE);

  /* w.r.t. amp */
  /* remove pixels with NaN (which means they are out of the cropped image */
  /* area) */
  ixstart = J->size[1];
  i0 = b_J->size[0] * b_J->size[1];
  b_J->size[0] = indxPixel->size[0];
  b_J->size[1] = ixstart;
  emxEnsureCapacity((emxArray__common *)b_J, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < ixstart; i0++) {
    loop_ub = indxPixel->size[0];
    for (i = 0; i < loop_ub; i++) {
      b_J->data[i + b_J->size[0] * i0] = J->data[((int32_T)indxPixel->data[i] +
        J->size[0] * i0) - 1];
    }
  }

  emxFree_real_T(&indxPixel);
  i0 = J->size[0] * J->size[1];
  J->size[0] = b_J->size[0];
  J->size[1] = b_J->size[1];
  emxEnsureCapacity((emxArray__common *)J, i0, (int32_T)sizeof(real_T));
  loop_ub = b_J->size[1];
  for (i0 = 0; i0 < loop_ub; i0++) {
    idx = b_J->size[0];
    for (i = 0; i < idx; i++) {
      J->data[i + J->size[0] * i0] = b_J->data[i + b_J->size[0] * i0];
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
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (fitNGaussians3D_mexCode.c) */
