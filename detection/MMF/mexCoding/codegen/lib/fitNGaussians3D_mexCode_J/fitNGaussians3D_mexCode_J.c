/*
 * fitNGaussians3D_mexCode_J.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode_J'
 *
 * C source code generated on: Mon May 28 14:49:14 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_J.h"
#include "squeeze.h"
#include "GaussListND_mexCode.h"
#include "fitNGaussians3D_mexCode_J_emxutil.h"
#include "exp.h"
#include "rdivide.h"
#include "power.h"
#include "repmat.h"
#include "fitNGaussians3D_mexCode_J_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y);

/* Function Definitions */

/*
 *
 */
static void eml_li_find(const emxArray_boolean_T *x, emxArray_int32_T *y)
{
  int32_T n;
  int32_T k;
  int32_T i;
  int32_T j;
  n = x->size[0];
  k = 0;
  for (i = 1; i <= n; i++) {
    if (x->data[i - 1]) {
      k++;
    }
  }

  j = y->size[0];
  y->size[0] = k;
  emxEnsureCapacity((emxArray__common *)y, j, (int32_T)sizeof(int32_T));
  j = 0;
  for (i = 1; i <= n; i++) {
    if (x->data[i - 1]) {
      y->data[j] = i;
      j++;
    }
  }
}

/*
 * function J = fitNGaussians3D_mexCode_J(x0,image,index,psfSigma)
 *  Edit of fitNGaussians2D to work in 3D
 *    EHarry March 2012
 */
void fitNGaussians3D_mexCode_J(emxArray_real_T *x0, const emxArray_real_T *image,
  const emxArray_real_T *b_index, const real_T psfSigma[2], emxArray_real_T *J)
{
  int32_T i0;
  emxArray_real_T *b_x0;
  int32_T n;
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
  boolean_T exitg6;
  emxArray_int32_T *r1;
  real_T b_mtmp;
  boolean_T exitg5;
  emxArray_int32_T *r2;
  real_T c_mtmp;
  boolean_T exitg4;
  emxArray_int32_T *r3;
  real_T d_mtmp;
  boolean_T exitg3;
  emxArray_int32_T *r4;
  real_T e_mtmp;
  boolean_T exitg2;
  emxArray_int32_T *r5;
  real_T f_mtmp;
  boolean_T exitg1;
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
  emxArray_real_T *relIndxX;
  emxArray_real_T *b_relIndxX;
  emxArray_real_T *c_relIndxX;
  emxArray_real_T *e_y;
  emxArray_real_T *psfValueY;
  emxArray_real_T *d_relIndxX;
  emxArray_real_T *e_relIndxX;
  emxArray_real_T *f_y;
  emxArray_real_T *psfValueZ;
  emxArray_real_T *f_relIndxX;
  emxArray_real_T *g_relIndxX;
  emxArray_real_T *g_y;
  emxArray_real_T *relIndxY;
  emxArray_real_T *relIndxZ;
  uint32_T h_y;
  emxArray_real_T *e_x0;
  emxArray_real_T *r6;
  emxArray_real_T *f_x0;
  emxArray_real_T *g_x0;
  emxArray_boolean_T *r7;
  emxArray_boolean_T *r8;
  emxArray_int32_T *r9;
  emxArray_real_T *b_J;

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
  /* 'fitNGaussians3D_mexCode_J:28' J = []; */
  /* % Input */
  /*  %check whether correct number of input arguments was used */
  /*  if nargin ~= 4 */
  /*      disp('--fitNGaussians2D: Incorrect number of input arguments!'); */
  /*      return */
  /*  end */
  /* check whether correct number of input arguments was used */
  /* 'fitNGaussians3D_mexCode_J:38' if nargin ~= 4 */
  /* % Calculating F & J */
  /* extract background intensity from x0 and remove from vector */
  /*  bgAmp = x0(end); */
  /* 'fitNGaussians3D_mexCode_J:47' x0 = x0(1:end-1); */
  if (1 > x0->size[0] - 1) {
    i0 = -1;
  } else {
    i0 = x0->size[0] - 2;
  }

  c_emxInit_real_T(&b_x0, 1);
  n = b_x0->size[0];
  b_x0->size[0] = i0 + 1;
  emxEnsureCapacity((emxArray__common *)b_x0, n, (int32_T)sizeof(real_T));
  for (n = 0; n <= i0; n++) {
    b_x0->data[n] = x0->data[n];
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
  /* 'fitNGaussians3D_mexCode_J:51' numPSF = length(x0)/4; */
  numPSF = (real_T)x0->size[0] / 4.0;

  /*  %reshape 3nx1 vector x0 into nx3 matrix */
  /*  x0 = reshape(x0,3,numPSF); */
  /*  x0 = x0'; */
  /* reshape 4nx1 vector x0 into nx4 matrix */
  /* 'fitNGaussians3D_mexCode_J:58' x0 = reshape(x0,4,numPSF); */
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

  /* 'fitNGaussians3D_mexCode_J:59' x0 = x0'; */
  i0 = d_x0->size[0] * d_x0->size[1];
  d_x0->size[0] = c_x0->size[1];
  d_x0->size[1] = 4;
  emxEnsureCapacity((emxArray__common *)d_x0, i0, (int32_T)sizeof(real_T));
  for (i0 = 0; i0 < 4; i0++) {
    loop_ub = c_x0->size[1] - 1;
    for (n = 0; n <= loop_ub; n++) {
      d_x0->data[n + d_x0->size[0] * i0] = c_x0->data[i0 + c_x0->size[0] * n];
    }
  }

  emxFree_real_T(&c_x0);
  emxInit_int32_T(&r0, 1);

  /* extract PSF center positions and amplitudes */
  /*  psfPos = x0(:,1:2); */
  /*  psfAmp = x0(:,3); */
  /* 'fitNGaussians3D_mexCode_J:65' psfPos = x0(:,1:3); */
  /* 'fitNGaussians3D_mexCode_J:66' psfAmp = x0(:,4); */
  /* find minimum and maximum pixel indices */
  /* 'fitNGaussians3D_mexCode_J:69' minIndxX = min(index(:,1)); */
  ixstart = 1;
  i0 = b_index->size[0];
  n = r0->size[0];
  r0->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r0, n, (int32_T)sizeof(int32_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r0->data[i0] = 1 + i0;
  }

  n = r0->size[0];
  mtmp = b_index->data[0];
  emxFree_int32_T(&r0);
  if (n > 1) {
    if (rtIsNaN(b_index->data[0])) {
      nx = 2;
      exitg6 = FALSE;
      while ((exitg6 == 0U) && (nx <= n)) {
        ixstart = nx;
        if (!rtIsNaN(b_index->data[nx - 1])) {
          mtmp = b_index->data[nx - 1];
          exitg6 = TRUE;
        } else {
          nx++;
        }
      }
    }

    if (ixstart < n) {
      while (ixstart + 1 <= n) {
        if (b_index->data[ixstart] < mtmp) {
          mtmp = b_index->data[ixstart];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r1, 1);

  /* 'fitNGaussians3D_mexCode_J:70' maxIndxX = max(index(:,1)); */
  ixstart = 1;
  i0 = b_index->size[0];
  n = r1->size[0];
  r1->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r1, n, (int32_T)sizeof(int32_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r1->data[i0] = 1 + i0;
  }

  n = r1->size[0];
  b_mtmp = b_index->data[0];
  emxFree_int32_T(&r1);
  if (n > 1) {
    if (rtIsNaN(b_index->data[0])) {
      nx = 2;
      exitg5 = FALSE;
      while ((exitg5 == 0U) && (nx <= n)) {
        ixstart = nx;
        if (!rtIsNaN(b_index->data[nx - 1])) {
          b_mtmp = b_index->data[nx - 1];
          exitg5 = TRUE;
        } else {
          nx++;
        }
      }
    }

    if (ixstart < n) {
      while (ixstart + 1 <= n) {
        if (b_index->data[ixstart] > b_mtmp) {
          b_mtmp = b_index->data[ixstart];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r2, 1);

  /* 'fitNGaussians3D_mexCode_J:71' minIndxY = min(index(:,2)); */
  ixstart = 1;
  i0 = b_index->size[0];
  n = r2->size[0];
  r2->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r2, n, (int32_T)sizeof(int32_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r2->data[i0] = 1 + i0;
  }

  n = r2->size[0];
  c_mtmp = b_index->data[b_index->size[0]];
  emxFree_int32_T(&r2);
  if (n > 1) {
    if (rtIsNaN(b_index->data[b_index->size[0]])) {
      nx = 2;
      exitg4 = FALSE;
      while ((exitg4 == 0U) && (nx <= n)) {
        ixstart = nx;
        if (!rtIsNaN(b_index->data[(nx + b_index->size[0]) - 1])) {
          c_mtmp = b_index->data[(nx + b_index->size[0]) - 1];
          exitg4 = TRUE;
        } else {
          nx++;
        }
      }
    }

    if (ixstart < n) {
      while (ixstart + 1 <= n) {
        if (b_index->data[ixstart + b_index->size[0]] < c_mtmp) {
          c_mtmp = b_index->data[ixstart + b_index->size[0]];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r3, 1);

  /* 'fitNGaussians3D_mexCode_J:72' maxIndxY = max(index(:,2)); */
  ixstart = 1;
  i0 = b_index->size[0];
  n = r3->size[0];
  r3->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r3, n, (int32_T)sizeof(int32_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r3->data[i0] = 1 + i0;
  }

  n = r3->size[0];
  d_mtmp = b_index->data[b_index->size[0]];
  emxFree_int32_T(&r3);
  if (n > 1) {
    if (rtIsNaN(b_index->data[b_index->size[0]])) {
      nx = 2;
      exitg3 = FALSE;
      while ((exitg3 == 0U) && (nx <= n)) {
        ixstart = nx;
        if (!rtIsNaN(b_index->data[(nx + b_index->size[0]) - 1])) {
          d_mtmp = b_index->data[(nx + b_index->size[0]) - 1];
          exitg3 = TRUE;
        } else {
          nx++;
        }
      }
    }

    if (ixstart < n) {
      while (ixstart + 1 <= n) {
        if (b_index->data[ixstart + b_index->size[0]] > d_mtmp) {
          d_mtmp = b_index->data[ixstart + b_index->size[0]];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r4, 1);

  /* 'fitNGaussians3D_mexCode_J:73' minIndxZ = min(index(:,3)); */
  ixstart = 1;
  i0 = b_index->size[0];
  n = r4->size[0];
  r4->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r4, n, (int32_T)sizeof(int32_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r4->data[i0] = 1 + i0;
  }

  n = r4->size[0];
  e_mtmp = b_index->data[b_index->size[0] << 1];
  emxFree_int32_T(&r4);
  if (n > 1) {
    if (rtIsNaN(b_index->data[b_index->size[0] << 1])) {
      nx = 2;
      exitg2 = FALSE;
      while ((exitg2 == 0U) && (nx <= n)) {
        ixstart = nx;
        if (!rtIsNaN(b_index->data[(nx + (b_index->size[0] << 1)) - 1])) {
          e_mtmp = b_index->data[(nx + (b_index->size[0] << 1)) - 1];
          exitg2 = TRUE;
        } else {
          nx++;
        }
      }
    }

    if (ixstart < n) {
      while (ixstart + 1 <= n) {
        if (b_index->data[ixstart + (b_index->size[0] << 1)] < e_mtmp) {
          e_mtmp = b_index->data[ixstart + (b_index->size[0] << 1)];
        }

        ixstart++;
      }
    }
  }

  emxInit_int32_T(&r5, 1);

  /* 'fitNGaussians3D_mexCode_J:74' maxIndxZ = max(index(:,3)); */
  ixstart = 1;
  i0 = b_index->size[0];
  n = r5->size[0];
  r5->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)r5, n, (int32_T)sizeof(int32_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r5->data[i0] = 1 + i0;
  }

  n = r5->size[0];
  f_mtmp = b_index->data[b_index->size[0] << 1];
  emxFree_int32_T(&r5);
  if (n > 1) {
    if (rtIsNaN(b_index->data[b_index->size[0] << 1])) {
      nx = 2;
      exitg1 = FALSE;
      while ((exitg1 == 0U) && (nx <= n)) {
        ixstart = nx;
        if (!rtIsNaN(b_index->data[(nx + (b_index->size[0] << 1)) - 1])) {
          f_mtmp = b_index->data[(nx + (b_index->size[0] << 1)) - 1];
          exitg1 = TRUE;
        } else {
          nx++;
        }
      }
    }

    if (ixstart < n) {
      while (ixstart + 1 <= n) {
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
  /* 'fitNGaussians3D_mexCode_J:78' psfIntegX = zeros(maxIndxX-minIndxX+1,numPSF); */
  i0 = psfIntegX->size[0] * psfIntegX->size[1];
  psfIntegX->size[0] = (int32_T)((b_mtmp - mtmp) + 1.0);
  psfIntegX->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfIntegX, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((b_mtmp - mtmp) + 1.0) * (int32_T)numPSF - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfIntegX->data[i0] = 0.0;
  }

  /* 'fitNGaussians3D_mexCode_J:79' for i=1:numPSF */
  i = 0;
  b_emxInit_real_T(&temp, 3);
  emxInit_real_T(&temp2, 2);
  emxInit_real_T(&y, 2);
  c_emxInit_real_T(&b_y, 1);
  while (i <= (int32_T)numPSF - 1) {
    /* 'fitNGaussians3D_mexCode_J:80' temp = GaussListND_mexCode((minIndxX:maxIndxX)',... */
    /* 'fitNGaussians3D_mexCode_J:81'         psfSigma(1),psfPos(i,1)); */
    if (rtIsNaN(mtmp) || rtIsNaN(b_mtmp)) {
      n = 1;
      anew = rtNaN;
      apnd = b_mtmp;
    } else if (b_mtmp < mtmp) {
      n = 0;
      anew = mtmp;
      apnd = b_mtmp;
    } else if (rtIsInf(mtmp) || rtIsInf(b_mtmp)) {
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

      if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
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
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (n > 0) {
      y->data[0] = anew;
      if (n > 1) {
        y->data[n - 1] = apnd;
        nx = n - 1;
        ixstart = nx / 2;
        for (k = 1; k <= ixstart - 1; k++) {
          y->data[k] = anew + (real_T)k;
          y->data[(n - k) - 1] = apnd - (real_T)k;
        }

        if (ixstart << 1 == nx) {
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

    /* 'fitNGaussians3D_mexCode_J:83' temp2 = squeeze(temp); */
    squeeze(temp, temp2);

    /* clear temp */
    /* 'fitNGaussians3D_mexCode_J:86' psfIntegX(:,i) = temp2(:,1); */
    i0 = temp2->size[0] - 1;
    for (n = 0; n <= i0; n++) {
      psfIntegX->data[n + psfIntegX->size[0] * i] = temp2->data[n];
    }

    /* clear temp2 */
    i++;
  }

  emxFree_real_T(&b_y);
  emxInit_real_T(&psfIntegY, 2);

  /* determine the contribution of each PSF (assuming amplitude 1) to a */
  /* pixel based on its y-coordinate (needed to calculate F & J) */
  /* 'fitNGaussians3D_mexCode_J:92' psfIntegY = zeros(maxIndxY-minIndxY+1,numPSF); */
  i0 = psfIntegY->size[0] * psfIntegY->size[1];
  psfIntegY->size[0] = (int32_T)((d_mtmp - c_mtmp) + 1.0);
  psfIntegY->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfIntegY, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((d_mtmp - c_mtmp) + 1.0) * (int32_T)numPSF - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfIntegY->data[i0] = 0.0;
  }

  /* 'fitNGaussians3D_mexCode_J:93' for i=1:numPSF */
  i = 0;
  c_emxInit_real_T(&c_y, 1);
  while (i <= (int32_T)numPSF - 1) {
    /* 'fitNGaussians3D_mexCode_J:94' temp = GaussListND_mexCode((minIndxY:maxIndxY)',... */
    /* 'fitNGaussians3D_mexCode_J:95'         psfSigma(1),psfPos(i,2)); */
    if (rtIsNaN(c_mtmp) || rtIsNaN(d_mtmp)) {
      n = 1;
      anew = rtNaN;
      apnd = d_mtmp;
    } else if (d_mtmp < c_mtmp) {
      n = 0;
      anew = c_mtmp;
      apnd = d_mtmp;
    } else if (rtIsInf(c_mtmp) || rtIsInf(d_mtmp)) {
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

      if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
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
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (n > 0) {
      y->data[0] = anew;
      if (n > 1) {
        y->data[n - 1] = apnd;
        nx = n - 1;
        ixstart = nx / 2;
        for (k = 1; k <= ixstart - 1; k++) {
          y->data[k] = anew + (real_T)k;
          y->data[(n - k) - 1] = apnd - (real_T)k;
        }

        if (ixstart << 1 == nx) {
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

    /* 'fitNGaussians3D_mexCode_J:97' temp2 = squeeze(temp); */
    squeeze(temp, temp2);

    /* clear temp */
    /* 'fitNGaussians3D_mexCode_J:100' psfIntegY(:,i) = temp2(:,1); */
    i0 = temp2->size[0] - 1;
    for (n = 0; n <= i0; n++) {
      psfIntegY->data[n + psfIntegY->size[0] * i] = temp2->data[n];
    }

    /* clear temp2 */
    i++;
  }

  emxFree_real_T(&c_y);
  emxInit_real_T(&psfIntegZ, 2);

  /* determine the contribution of each PSF (assuming amplitude 1) to a */
  /* pixel based on its z-coordinate (needed to calculate F & J) */
  /* 'fitNGaussians3D_mexCode_J:106' psfIntegZ = zeros(maxIndxZ-minIndxZ+1,numPSF); */
  i0 = psfIntegZ->size[0] * psfIntegZ->size[1];
  psfIntegZ->size[0] = (int32_T)((f_mtmp - e_mtmp) + 1.0);
  psfIntegZ->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfIntegZ, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((f_mtmp - e_mtmp) + 1.0) * (int32_T)numPSF - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfIntegZ->data[i0] = 0.0;
  }

  /* 'fitNGaussians3D_mexCode_J:107' for i=1:numPSF */
  i = 0;
  c_emxInit_real_T(&d_y, 1);
  while (i <= (int32_T)numPSF - 1) {
    /* 'fitNGaussians3D_mexCode_J:108' temp = GaussListND_mexCode((minIndxZ:maxIndxZ)',... */
    /* 'fitNGaussians3D_mexCode_J:109'         psfSigma(2),psfPos(i,3)); */
    if (rtIsNaN(e_mtmp) || rtIsNaN(f_mtmp)) {
      n = 1;
      anew = rtNaN;
      apnd = f_mtmp;
    } else if (f_mtmp < e_mtmp) {
      n = 0;
      anew = e_mtmp;
      apnd = f_mtmp;
    } else if (rtIsInf(e_mtmp) || rtIsInf(f_mtmp)) {
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

      if (fabs(cdiff) < 4.4408920985006262E-16 * absb) {
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
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (n > 0) {
      y->data[0] = anew;
      if (n > 1) {
        y->data[n - 1] = apnd;
        nx = n - 1;
        ixstart = nx / 2;
        for (k = 1; k <= ixstart - 1; k++) {
          y->data[k] = anew + (real_T)k;
          y->data[(n - k) - 1] = apnd - (real_T)k;
        }

        if (ixstart << 1 == nx) {
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

    /* 'fitNGaussians3D_mexCode_J:112' temp2 = squeeze(temp); */
    squeeze(temp, temp2);

    /* clear temp */
    /* 'fitNGaussians3D_mexCode_J:115' psfIntegZ(:,i) = temp2(:,1); */
    i0 = temp2->size[0] - 1;
    for (n = 0; n <= i0; n++) {
      psfIntegZ->data[n + psfIntegZ->size[0] * i] = temp2->data[n];
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
  /* 'fitNGaussians3D_mexCode_J:121' psfValueX = zeros(maxIndxX-minIndxX+2,numPSF); */
  i0 = psfValueX->size[0] * psfValueX->size[1];
  psfValueX->size[0] = (int32_T)((b_mtmp - mtmp) + 2.0);
  psfValueX->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfValueX, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((b_mtmp - mtmp) + 2.0) * (int32_T)numPSF - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfValueX->data[i0] = 0.0;
  }

  /* 'fitNGaussians3D_mexCode_J:122' for i=1:numPSF */
  i = 0;
  c_emxInit_real_T(&relIndxX, 1);
  c_emxInit_real_T(&b_relIndxX, 1);
  c_emxInit_real_T(&c_relIndxX, 1);
  c_emxInit_real_T(&e_y, 1);
  while (i <= (int32_T)numPSF - 1) {
    /* 'fitNGaussians3D_mexCode_J:123' psfValueX(:,i) = exp(-((minIndxX-0.5:maxIndxX+0.5)'... */
    /* 'fitNGaussians3D_mexCode_J:124'         -psfPos(i,1)).^2/2/psfSigma(1)^2); */
    if (rtIsNaN(mtmp - 0.5) || rtIsNaN(b_mtmp + 0.5)) {
      n = 1;
      anew = rtNaN;
      apnd = b_mtmp + 0.5;
    } else if (b_mtmp + 0.5 < mtmp - 0.5) {
      n = 0;
      anew = mtmp - 0.5;
      apnd = b_mtmp + 0.5;
    } else if (rtIsInf(mtmp - 0.5) || rtIsInf(b_mtmp + 0.5)) {
      n = 1;
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
        n = (int32_T)ndbl;
      } else {
        n = 0;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (n > 0) {
      y->data[0] = anew;
      if (n > 1) {
        y->data[n - 1] = apnd;
        nx = n - 1;
        ixstart = nx / 2;
        for (k = 1; k <= ixstart - 1; k++) {
          y->data[k] = anew + (real_T)k;
          y->data[(n - k) - 1] = apnd - (real_T)k;
        }

        if (ixstart << 1 == nx) {
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

    power(e_y, relIndxX);
    i0 = relIndxX->size[0];
    relIndxX->size[0] = relIndxX->size[0];
    emxEnsureCapacity((emxArray__common *)relIndxX, i0, (int32_T)sizeof(real_T));
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      relIndxX->data[i0] = -relIndxX->data[i0];
    }

    i0 = c_relIndxX->size[0];
    c_relIndxX->size[0] = relIndxX->size[0];
    emxEnsureCapacity((emxArray__common *)c_relIndxX, i0, (int32_T)sizeof(real_T));
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      c_relIndxX->data[i0] = relIndxX->data[i0];
    }

    rdivide(c_relIndxX, 2.0, relIndxX);
    i0 = b_relIndxX->size[0];
    b_relIndxX->size[0] = relIndxX->size[0];
    emxEnsureCapacity((emxArray__common *)b_relIndxX, i0, (int32_T)sizeof(real_T));
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      b_relIndxX->data[i0] = relIndxX->data[i0];
    }

    rdivide(b_relIndxX, rt_powd_snf(psfSigma[0], 2.0), relIndxX);
    b_exp(relIndxX);
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      psfValueX->data[i0 + psfValueX->size[0] * i] = relIndxX->data[i0];
    }

    i++;
  }

  emxFree_real_T(&e_y);
  emxFree_real_T(&c_relIndxX);
  emxFree_real_T(&b_relIndxX);
  emxInit_real_T(&psfValueY, 2);

  /* calculate the value of each PSF (assuming amplitude 1) at the */
  /* y-coordinates of the corners of all pixels (needed to calculate J) */
  /* 'fitNGaussians3D_mexCode_J:129' psfValueY = zeros(maxIndxY-minIndxY+2,numPSF); */
  i0 = psfValueY->size[0] * psfValueY->size[1];
  psfValueY->size[0] = (int32_T)((d_mtmp - c_mtmp) + 2.0);
  psfValueY->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfValueY, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((d_mtmp - c_mtmp) + 2.0) * (int32_T)numPSF - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfValueY->data[i0] = 0.0;
  }

  /* 'fitNGaussians3D_mexCode_J:130' for i=1:numPSF */
  i = 0;
  c_emxInit_real_T(&d_relIndxX, 1);
  c_emxInit_real_T(&e_relIndxX, 1);
  c_emxInit_real_T(&f_y, 1);
  while (i <= (int32_T)numPSF - 1) {
    /* 'fitNGaussians3D_mexCode_J:131' psfValueY(:,i) = exp(-((minIndxY-0.5:maxIndxY+0.5)'... */
    /* 'fitNGaussians3D_mexCode_J:132'         -psfPos(i,2)).^2/2/psfSigma(1)^2); */
    if (rtIsNaN(c_mtmp - 0.5) || rtIsNaN(d_mtmp + 0.5)) {
      n = 1;
      anew = rtNaN;
      apnd = d_mtmp + 0.5;
    } else if (d_mtmp + 0.5 < c_mtmp - 0.5) {
      n = 0;
      anew = c_mtmp - 0.5;
      apnd = d_mtmp + 0.5;
    } else if (rtIsInf(c_mtmp - 0.5) || rtIsInf(d_mtmp + 0.5)) {
      n = 1;
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
        n = (int32_T)ndbl;
      } else {
        n = 0;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (n > 0) {
      y->data[0] = anew;
      if (n > 1) {
        y->data[n - 1] = apnd;
        nx = n - 1;
        ixstart = nx / 2;
        for (k = 1; k <= ixstart - 1; k++) {
          y->data[k] = anew + (real_T)k;
          y->data[(n - k) - 1] = apnd - (real_T)k;
        }

        if (ixstart << 1 == nx) {
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

    power(f_y, relIndxX);
    i0 = relIndxX->size[0];
    relIndxX->size[0] = relIndxX->size[0];
    emxEnsureCapacity((emxArray__common *)relIndxX, i0, (int32_T)sizeof(real_T));
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      relIndxX->data[i0] = -relIndxX->data[i0];
    }

    i0 = e_relIndxX->size[0];
    e_relIndxX->size[0] = relIndxX->size[0];
    emxEnsureCapacity((emxArray__common *)e_relIndxX, i0, (int32_T)sizeof(real_T));
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      e_relIndxX->data[i0] = relIndxX->data[i0];
    }

    rdivide(e_relIndxX, 2.0, relIndxX);
    i0 = d_relIndxX->size[0];
    d_relIndxX->size[0] = relIndxX->size[0];
    emxEnsureCapacity((emxArray__common *)d_relIndxX, i0, (int32_T)sizeof(real_T));
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      d_relIndxX->data[i0] = relIndxX->data[i0];
    }

    rdivide(d_relIndxX, rt_powd_snf(psfSigma[0], 2.0), relIndxX);
    b_exp(relIndxX);
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      psfValueY->data[i0 + psfValueY->size[0] * i] = relIndxX->data[i0];
    }

    i++;
  }

  emxFree_real_T(&f_y);
  emxFree_real_T(&e_relIndxX);
  emxFree_real_T(&d_relIndxX);
  emxInit_real_T(&psfValueZ, 2);

  /* calculate the value of each PSF (assuming amplitude 1) at the */
  /* z-coordinates of the corners of all pixels (needed to calculate J) */
  /* 'fitNGaussians3D_mexCode_J:137' psfValueZ = zeros(maxIndxZ-minIndxZ+2,numPSF); */
  i0 = psfValueZ->size[0] * psfValueZ->size[1];
  psfValueZ->size[0] = (int32_T)((f_mtmp - e_mtmp) + 2.0);
  psfValueZ->size[1] = (int32_T)numPSF;
  emxEnsureCapacity((emxArray__common *)psfValueZ, i0, (int32_T)sizeof(real_T));
  loop_ub = (int32_T)((f_mtmp - e_mtmp) + 2.0) * (int32_T)numPSF - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    psfValueZ->data[i0] = 0.0;
  }

  /* 'fitNGaussians3D_mexCode_J:138' for i=1:numPSF */
  i = 0;
  c_emxInit_real_T(&f_relIndxX, 1);
  c_emxInit_real_T(&g_relIndxX, 1);
  c_emxInit_real_T(&g_y, 1);
  while (i <= (int32_T)numPSF - 1) {
    /* 'fitNGaussians3D_mexCode_J:139' psfValueZ(:,i) = exp(-((minIndxZ-0.5:maxIndxZ+0.5)'... */
    /* 'fitNGaussians3D_mexCode_J:140'         -psfPos(i,3)).^2/2/psfSigma(2)^2); */
    if (rtIsNaN(e_mtmp - 0.5) || rtIsNaN(f_mtmp + 0.5)) {
      n = 1;
      anew = rtNaN;
      apnd = f_mtmp + 0.5;
    } else if (f_mtmp + 0.5 < e_mtmp - 0.5) {
      n = 0;
      anew = e_mtmp - 0.5;
      apnd = f_mtmp + 0.5;
    } else if (rtIsInf(e_mtmp - 0.5) || rtIsInf(f_mtmp + 0.5)) {
      n = 1;
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
        n = (int32_T)ndbl;
      } else {
        n = 0;
      }
    }

    i0 = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = n;
    emxEnsureCapacity((emxArray__common *)y, i0, (int32_T)sizeof(real_T));
    if (n > 0) {
      y->data[0] = anew;
      if (n > 1) {
        y->data[n - 1] = apnd;
        nx = n - 1;
        ixstart = nx / 2;
        for (k = 1; k <= ixstart - 1; k++) {
          y->data[k] = anew + (real_T)k;
          y->data[(n - k) - 1] = apnd - (real_T)k;
        }

        if (ixstart << 1 == nx) {
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

    power(g_y, relIndxX);
    i0 = relIndxX->size[0];
    relIndxX->size[0] = relIndxX->size[0];
    emxEnsureCapacity((emxArray__common *)relIndxX, i0, (int32_T)sizeof(real_T));
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      relIndxX->data[i0] = -relIndxX->data[i0];
    }

    i0 = g_relIndxX->size[0];
    g_relIndxX->size[0] = relIndxX->size[0];
    emxEnsureCapacity((emxArray__common *)g_relIndxX, i0, (int32_T)sizeof(real_T));
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      g_relIndxX->data[i0] = relIndxX->data[i0];
    }

    rdivide(g_relIndxX, 2.0, relIndxX);
    i0 = f_relIndxX->size[0];
    f_relIndxX->size[0] = relIndxX->size[0];
    emxEnsureCapacity((emxArray__common *)f_relIndxX, i0, (int32_T)sizeof(real_T));
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      f_relIndxX->data[i0] = relIndxX->data[i0];
    }

    rdivide(f_relIndxX, rt_powd_snf(psfSigma[1], 2.0), relIndxX);
    b_exp(relIndxX);
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
      psfValueZ->data[i0 + psfValueZ->size[0] * i] = relIndxX->data[i0];
    }

    i++;
  }

  emxFree_real_T(&g_y);
  emxFree_real_T(&g_relIndxX);
  emxFree_real_T(&f_relIndxX);
  emxFree_real_T(&y);

  /* get number of pixels in image */
  /* 'fitNGaussians3D_mexCode_J:144' numPixel = length(image); */
  /*  %get xy-indices relative to minimum */
  /*  relIndxX = index(:,1) - minIndxX + 1; */
  /*  relIndxY = index(:,2) - minIndxY + 1; */
  /* get xy-indices relative to minimum */
  /* 'fitNGaussians3D_mexCode_J:150' relIndxX = index(:,1) - minIndxX + 1; */
  i0 = b_index->size[0];
  n = relIndxX->size[0];
  relIndxX->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)relIndxX, n, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    relIndxX->data[i0] = (b_index->data[i0] - mtmp) + 1.0;
  }

  c_emxInit_real_T(&relIndxY, 1);

  /* 'fitNGaussians3D_mexCode_J:151' relIndxY = index(:,2) - minIndxY + 1; */
  i0 = b_index->size[0];
  n = relIndxY->size[0];
  relIndxY->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)relIndxY, n, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    relIndxY->data[i0] = (b_index->data[i0 + b_index->size[0]] - c_mtmp) + 1.0;
  }

  c_emxInit_real_T(&relIndxZ, 1);

  /* 'fitNGaussians3D_mexCode_J:152' relIndxZ = index(:,3) - minIndxZ + 1; */
  i0 = b_index->size[0];
  n = relIndxZ->size[0];
  relIndxZ->size[0] = i0;
  emxEnsureCapacity((emxArray__common *)relIndxZ, n, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    relIndxZ->data[i0] = (b_index->data[i0 + (b_index->size[0] << 1)] - e_mtmp)
      + 1.0;
  }

  /*  %calculate the value of F at all pixels */
  /*  F = (sum(repmat(psfAmp,1,numPixel).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)',1))' ... */
  /*      + repmat(bgAmp,numPixel,1) - image; */
  /* calculate the value of F at all pixels */
  /*  F = (sum(repmat(psfAmp,1,numPixel).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)'.*psfIntegZ(relIndxZ,:)',1))' ... */
  /*      + repmat(bgAmp,numPixel,1) - image; */
  /* remove pixels with NaN (which means they are out of the cropped image */
  /* area) */
  /* 'fitNGaussians3D_mexCode_J:164' indxPixel = ~isnan(image); */
  /*  F = F(indxPixel); */
  /*  if nargout > 1 */
  /* calculate the derivative at all pixels */
  /*      J = ones(numPixel,3*numPSF+1); %(last column for background amplitude) */
  /*      J(:,1:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueX(relIndxX,:)-... */
  /*          psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:); %w.r.t. x */
  /*      J(:,2:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueY(relIndxY,:)-... */
  /*          psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:); %w.r.t. y */
  /*      J(:,3:3:3*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); %w.r.t. amp */
  /* 'fitNGaussians3D_mexCode_J:177' J = ones(numPixel,4*numPSF+1); */
  h_y = (uint32_T)numPSF << 2;
  nx = image->size[0];
  i0 = J->size[0] * J->size[1];
  J->size[0] = nx;
  emxEnsureCapacity((emxArray__common *)J, i0, (int32_T)sizeof(real_T));
  i0 = J->size[0] * J->size[1];
  J->size[1] = (int32_T)((real_T)h_y + 1.0);
  emxEnsureCapacity((emxArray__common *)J, i0, (int32_T)sizeof(real_T));
  loop_ub = image->size[0] * (int32_T)((real_T)h_y + 1.0) - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    J->data[i0] = 1.0;
  }

  emxInit_real_T(&e_x0, 2);

  /* (last column for background amplitude) */
  /* 'fitNGaussians3D_mexCode_J:178' J(:,1:4:4*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueX(relIndxX,:)-... */
  /* 'fitNGaussians3D_mexCode_J:179'     psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:).*psfIntegZ(relIndxZ,:); */
  i0 = d_x0->size[0];
  n = e_x0->size[0] * e_x0->size[1];
  e_x0->size[0] = 1;
  e_x0->size[1] = i0;
  emxEnsureCapacity((emxArray__common *)e_x0, n, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    e_x0->data[e_x0->size[0] * i0] = d_x0->data[i0 + d_x0->size[0] * 3];
  }

  emxInit_real_T(&r6, 2);
  repmat(e_x0, (real_T)image->size[0], r6);
  emxFree_real_T(&e_x0);
  if (1U > ((uint32_T)numPSF << 2)) {
    i0 = 1;
  } else {
    i0 = 4;
  }

  loop_ub = r6->size[1] - 1;
  for (n = 0; n <= loop_ub; n++) {
    nx = r6->size[0] - 1;
    for (k = 0; k <= nx; k++) {
      J->data[k + J->size[0] * (i0 * n)] = r6->data[k + r6->size[0] * n] *
        (psfValueX->data[((int32_T)relIndxX->data[k] + psfValueX->size[0] * n) -
         1] - psfValueX->data[((int32_T)(relIndxX->data[k] + 1.0) +
          psfValueX->size[0] * n) - 1]) * psfIntegY->data[((int32_T)
        relIndxY->data[k] + psfIntegY->size[0] * n) - 1] * psfIntegZ->data
        [((int32_T)relIndxZ->data[k] + psfIntegZ->size[0] * n) - 1];
    }
  }

  emxFree_real_T(&psfValueX);
  emxInit_real_T(&f_x0, 2);

  /* w.r.t. x */
  /* 'fitNGaussians3D_mexCode_J:180' J(:,2:4:4*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueY(relIndxY,:)-... */
  /* 'fitNGaussians3D_mexCode_J:181'     psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:).*psfIntegZ(relIndxZ,:); */
  i0 = d_x0->size[0];
  n = f_x0->size[0] * f_x0->size[1];
  f_x0->size[0] = 1;
  f_x0->size[1] = i0;
  emxEnsureCapacity((emxArray__common *)f_x0, n, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    f_x0->data[f_x0->size[0] * i0] = d_x0->data[i0 + d_x0->size[0] * 3];
  }

  repmat(f_x0, (real_T)image->size[0], r6);
  emxFree_real_T(&f_x0);
  if (2U > ((uint32_T)numPSF << 2)) {
    i0 = 0;
    n = 1;
  } else {
    i0 = 1;
    n = 4;
  }

  loop_ub = r6->size[1] - 1;
  for (k = 0; k <= loop_ub; k++) {
    nx = r6->size[0] - 1;
    for (ixstart = 0; ixstart <= nx; ixstart++) {
      J->data[ixstart + J->size[0] * (i0 + n * k)] = r6->data[ixstart + r6->
        size[0] * k] * (psfValueY->data[((int32_T)relIndxY->data[ixstart] +
        psfValueY->size[0] * k) - 1] - psfValueY->data[((int32_T)(relIndxY->
        data[ixstart] + 1.0) + psfValueY->size[0] * k) - 1]) * psfIntegX->data
        [((int32_T)relIndxX->data[ixstart] + psfIntegX->size[0] * k) - 1] *
        psfIntegZ->data[((int32_T)relIndxZ->data[ixstart] + psfIntegZ->size[0] *
                         k) - 1];
    }
  }

  emxFree_real_T(&psfValueY);
  emxInit_real_T(&g_x0, 2);

  /* w.r.t. y */
  /* 'fitNGaussians3D_mexCode_J:182' J(:,3:4:4*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueZ(relIndxZ,:)-... */
  /* 'fitNGaussians3D_mexCode_J:183'     psfValueZ(relIndxZ+1,:)).*psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); */
  i0 = d_x0->size[0];
  n = g_x0->size[0] * g_x0->size[1];
  g_x0->size[0] = 1;
  g_x0->size[1] = i0;
  emxEnsureCapacity((emxArray__common *)g_x0, n, (int32_T)sizeof(real_T));
  loop_ub = i0 - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    g_x0->data[g_x0->size[0] * i0] = d_x0->data[i0 + d_x0->size[0] * 3];
  }

  emxFree_real_T(&d_x0);
  repmat(g_x0, (real_T)image->size[0], r6);
  emxFree_real_T(&g_x0);
  if (3U > ((uint32_T)numPSF << 2)) {
    i0 = 0;
    n = 1;
  } else {
    i0 = 2;
    n = 4;
  }

  loop_ub = r6->size[1] - 1;
  for (k = 0; k <= loop_ub; k++) {
    nx = r6->size[0] - 1;
    for (ixstart = 0; ixstart <= nx; ixstart++) {
      J->data[ixstart + J->size[0] * (i0 + n * k)] = r6->data[ixstart + r6->
        size[0] * k] * (psfValueZ->data[((int32_T)relIndxZ->data[ixstart] +
        psfValueZ->size[0] * k) - 1] - psfValueZ->data[((int32_T)(relIndxZ->
        data[ixstart] + 1.0) + psfValueZ->size[0] * k) - 1]) * psfIntegX->data
        [((int32_T)relIndxX->data[ixstart] + psfIntegX->size[0] * k) - 1] *
        psfIntegY->data[((int32_T)relIndxY->data[ixstart] + psfIntegY->size[0] *
                         k) - 1];
    }
  }

  emxFree_real_T(&r6);
  emxFree_real_T(&psfValueZ);

  /* w.r.t. z */
  /* 'fitNGaussians3D_mexCode_J:184' J(:,4:4:4*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:).*psfIntegZ(relIndxZ,:); */
  if (4U > ((uint32_T)numPSF << 2)) {
    i0 = 0;
    n = 1;
  } else {
    i0 = 3;
    n = 4;
  }

  k = psfIntegX->size[1] - 1;
  for (ixstart = 0; ixstart <= k; ixstart++) {
    loop_ub = relIndxX->size[0] - 1;
    for (nx = 0; nx <= loop_ub; nx++) {
      J->data[nx + J->size[0] * (i0 + n * ixstart)] = psfIntegX->data[((int32_T)
        relIndxX->data[nx] + psfIntegX->size[0] * ixstart) - 1] *
        psfIntegY->data[((int32_T)relIndxY->data[nx] + psfIntegY->size[0] *
                         ixstart) - 1] * psfIntegZ->data[((int32_T)
        relIndxZ->data[nx] + psfIntegZ->size[0] * ixstart) - 1];
    }
  }

  emxFree_real_T(&relIndxZ);
  emxFree_real_T(&relIndxY);
  emxFree_real_T(&relIndxX);
  emxFree_real_T(&psfIntegZ);
  emxFree_real_T(&psfIntegY);
  emxFree_real_T(&psfIntegX);
  emxInit_boolean_T(&r7, 1);

  /* w.r.t. amp */
  /* remove pixels with NaN (which means they are out of the cropped image */
  /* area) */
  /* 'fitNGaussians3D_mexCode_J:188' J = J(indxPixel,:); */
  i0 = r7->size[0];
  r7->size[0] = image->size[0];
  emxEnsureCapacity((emxArray__common *)r7, i0, (int32_T)sizeof(boolean_T));
  loop_ub = image->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r7->data[i0] = rtIsNaN(image->data[i0]);
  }

  emxInit_boolean_T(&r8, 1);
  i0 = r8->size[0];
  r8->size[0] = r7->size[0];
  emxEnsureCapacity((emxArray__common *)r8, i0, (int32_T)sizeof(boolean_T));
  loop_ub = r7->size[0] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    r8->data[i0] = !r7->data[i0];
  }

  emxFree_boolean_T(&r7);
  emxInit_int32_T(&r9, 1);
  emxInit_real_T(&b_J, 2);
  eml_li_find(r8, r9);
  nx = J->size[1];
  i0 = b_J->size[0] * b_J->size[1];
  b_J->size[0] = r9->size[0];
  b_J->size[1] = nx;
  emxEnsureCapacity((emxArray__common *)b_J, i0, (int32_T)sizeof(real_T));
  emxFree_boolean_T(&r8);
  loop_ub = nx - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    nx = r9->size[0] - 1;
    for (n = 0; n <= nx; n++) {
      b_J->data[n + b_J->size[0] * i0] = J->data[(r9->data[n] + J->size[0] * i0)
        - 1];
    }
  }

  emxFree_int32_T(&r9);
  i0 = J->size[0] * J->size[1];
  J->size[0] = b_J->size[0];
  J->size[1] = b_J->size[1];
  emxEnsureCapacity((emxArray__common *)J, i0, (int32_T)sizeof(real_T));
  loop_ub = b_J->size[1] - 1;
  for (i0 = 0; i0 <= loop_ub; i0++) {
    nx = b_J->size[0] - 1;
    for (n = 0; n <= nx; n++) {
      J->data[n + J->size[0] * i0] = b_J->data[n + b_J->size[0] * i0];
    }
  }

  emxFree_real_T(&b_J);

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

/* End of code generation (fitNGaussians3D_mexCode_J.c) */
