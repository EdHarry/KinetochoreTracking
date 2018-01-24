/*
 * fitNGaussians3D_mexCode_F.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode_F'
 *
 * C source code generated on: Mon May 21 17:33:35 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_F.h"
#include "squeeze.h"
#include "GaussListND_mexCode.h"
#include "fitNGaussians3D_mexCode_F_emxutil.h"
#include "max.h"
#include "min.h"
#include "length.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
static emlrtDCInfo emlrtDCI = { 79, 19, "fitNGaussians3D_mexCode_F", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_F.m", 1 };
static emlrtDCInfo b_emlrtDCI = { 79, 39, "fitNGaussians3D_mexCode_F", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_F.m", 1 };
static emlrtDCInfo c_emlrtDCI = { 93, 19, "fitNGaussians3D_mexCode_F", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_F.m", 1 };
static emlrtDCInfo d_emlrtDCI = { 107, 19, "fitNGaussians3D_mexCode_F", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_F.m", 1 };
static emlrtDCInfo e_emlrtDCI = { 79, 19, "fitNGaussians3D_mexCode_F", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_F.m", 1 };
static emlrtDCInfo f_emlrtDCI = { 79, 39, "fitNGaussians3D_mexCode_F", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_F.m", 1 };
static emlrtDCInfo g_emlrtDCI = { 93, 19, "fitNGaussians3D_mexCode_F", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_F.m", 1 };
static emlrtDCInfo h_emlrtDCI = { 107, 19, "fitNGaussians3D_mexCode_F", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_F.m", 1 };
static emlrtDCInfo i_emlrtDCI = { 160, 47, "fitNGaussians3D_mexCode_F", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_F.m", 1 };
static emlrtDCInfo j_emlrtDCI = { 160, 71, "fitNGaussians3D_mexCode_F", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_F.m", 1 };
static emlrtDCInfo k_emlrtDCI = { 160, 95, "fitNGaussians3D_mexCode_F", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_F.m", 1 };

/* Function Declarations */

/* Function Definitions */

/*
 * function F = fitNGaussians3D_mexCode_F(x0,image,index,psfSigma)
 */
void fitNGaussians3D_mexCode_F(emxArray_real_T *x0, const emxArray_real_T *image, const emxArray_real_T *b_index, const real_T psfSigma[2], emxArray_real_T *F)
{
    real_T bgAmp;
    int32_T i;
    emxArray_real_T *b_x0;
    int32_T nx;
    int32_T ia;
    real_T numPSF;
    int32_T sz[2];
    emxArray_real_T *c_x0;
    real_T minIndxX;
    int32_T k;
    emxArray_real_T *d_x0;
    emxArray_real_T *c_index;
    emxArray_real_T *d_index;
    emxArray_real_T *e_index;
    real_T maxIndxX;
    emxArray_real_T *f_index;
    real_T minIndxY;
    emxArray_real_T *g_index;
    real_T maxIndxY;
    emxArray_real_T *h_index;
    real_T minIndxZ;
    emxArray_real_T *psfIntegX;
    real_T maxIndxZ;
    uint32_T b_i;
    emxArray_real_T *temp;
    emxArray_real_T *temp2;
    emxArray_real_T *y;
    emxArray_real_T *b_y;
    int32_T iy;
    real_T anew;
    real_T apnd;
    real_T ndbl;
    real_T cdiff;
    real_T absa;
    real_T absb;
    int32_T ib;
    emxArray_real_T *psfIntegY;
    emxArray_real_T *c_y;
    emxArray_real_T *psfIntegZ;
    emxArray_real_T *d_y;
    emxArray_int32_T *r0;
    int32_T iv0[2];
    int32_T iv1[2];
    emxArray_real_T *x;
    emxArray_int32_T *r1;
    emxArray_real_T *r2;
    uint32_T b_sz[2];
    emxArray_real_T *r3;
    emxArray_boolean_T *indxPixel;
    emxArray_int32_T *r4;
    emxArray_real_T *b_F;
    emxArray_real_T *c_F;
    int32_T d_F[2];
    emxArray_real_T e_F;
    emlrtHeapReferenceStackEnterFcn();
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
    i = x0->size[0] - 1;
    if (1 > i) {
        i = 0;
    }
    c_emxInit_real_T(&b_x0, 1, TRUE);
    nx = b_x0->size[0];
    b_x0->size[0] = i;
    emxEnsureCapacity((emxArray__common *)b_x0, nx, (int32_T)sizeof(real_T));
    ia = i - 1;
    for (i = 0; i <= ia; i++) {
        b_x0->data[i] = x0->data[i];
    }
    i = x0->size[0];
    x0->size[0] = b_x0->size[0];
    emxEnsureCapacity((emxArray__common *)x0, i, (int32_T)sizeof(real_T));
    ia = b_x0->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        x0->data[i] = b_x0->data[i];
    }
    emxFree_real_T(&b_x0);
    /* get number of PSFs considered */
    /*  numPSF = length(x0)/3; */
    /* 'fitNGaussians3D_mexCode_F:52' numPSF = length(x0)/4; */
    numPSF = length(x0) / 4.0;
    /*  %reshape 3nx1 vector x0 into nx3 matrix */
    /*  x0 = reshape(x0,3,numPSF); */
    /*  x0 = x0'; */
    /* reshape 4nx1 vector x0 into nx4 matrix */
    /* 'fitNGaussians3D_mexCode_F:59' x0 = reshape(x0,4,numPSF); */
    nx = x0->size[0];
    for (i = 0; i < 2; i++) {
        sz[i] = 0;
    }
    emxInit_real_T(&c_x0, 2, TRUE);
    sz[0] = 4;
    minIndxX = numPSF;
    minIndxX = minIndxX < 0.0 ? muDoubleScalarCeil(minIndxX - 0.5) : muDoubleScalarFloor(minIndxX + 0.5);
    if (minIndxX < 2.147483648E+9) {
        if (minIndxX >= -2.147483648E+9) {
            i = (int32_T)minIndxX;
        } else {
            i = MIN_int32_T;
        }
    } else if (minIndxX >= 2.147483648E+9) {
        i = MAX_int32_T;
    } else {
        i = 0;
    }
    sz[1] = i;
    i = c_x0->size[0] * c_x0->size[1];
    c_x0->size[0] = 4;
    c_x0->size[1] = sz[1];
    emxEnsureCapacity((emxArray__common *)c_x0, i, (int32_T)sizeof(real_T));
    for (k = 0; k + 1 <= nx; k++) {
        c_x0->data[k] = x0->data[k];
    }
    emxInit_real_T(&d_x0, 2, TRUE);
    /* 'fitNGaussians3D_mexCode_F:60' x0 = x0'; */
    i = d_x0->size[0] * d_x0->size[1];
    d_x0->size[0] = c_x0->size[1];
    d_x0->size[1] = 4;
    emxEnsureCapacity((emxArray__common *)d_x0, i, (int32_T)sizeof(real_T));
    for (i = 0; i < 4; i++) {
        ia = c_x0->size[1] - 1;
        for (nx = 0; nx <= ia; nx++) {
            d_x0->data[nx + d_x0->size[0] * i] = c_x0->data[i + c_x0->size[0] * nx];
        }
    }
    emxFree_real_T(&c_x0);
    c_emxInit_real_T(&c_index, 1, TRUE);
    /* extract PSF center positions and amplitudes */
    /*  psfPos = x0(:,1:2); */
    /*  psfAmp = x0(:,3); */
    /* 'fitNGaussians3D_mexCode_F:66' psfPos = x0(:,1:3); */
    /* 'fitNGaussians3D_mexCode_F:67' psfAmp = x0(:,4); */
    /* find minimum and maximum pixel indices */
    /* 'fitNGaussians3D_mexCode_F:70' minIndxX = min(index(:,1)); */
    i = c_index->size[0];
    c_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)c_index, i, (int32_T)sizeof(real_T));
    ia = b_index->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        c_index->data[i] = b_index->data[i];
    }
    c_emxInit_real_T(&d_index, 1, TRUE);
    minIndxX = b_min(c_index);
    /* 'fitNGaussians3D_mexCode_F:71' maxIndxX = max(index(:,1)); */
    i = d_index->size[0];
    d_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)d_index, i, (int32_T)sizeof(real_T));
    emxFree_real_T(&c_index);
    ia = b_index->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        d_index->data[i] = b_index->data[i];
    }
    c_emxInit_real_T(&e_index, 1, TRUE);
    maxIndxX = b_max(d_index);
    /* 'fitNGaussians3D_mexCode_F:72' minIndxY = min(index(:,2)); */
    i = e_index->size[0];
    e_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)e_index, i, (int32_T)sizeof(real_T));
    emxFree_real_T(&d_index);
    ia = b_index->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        e_index->data[i] = b_index->data[i + b_index->size[0]];
    }
    c_emxInit_real_T(&f_index, 1, TRUE);
    minIndxY = b_min(e_index);
    /* 'fitNGaussians3D_mexCode_F:73' maxIndxY = max(index(:,2)); */
    i = f_index->size[0];
    f_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)f_index, i, (int32_T)sizeof(real_T));
    emxFree_real_T(&e_index);
    ia = b_index->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        f_index->data[i] = b_index->data[i + b_index->size[0]];
    }
    c_emxInit_real_T(&g_index, 1, TRUE);
    maxIndxY = b_max(f_index);
    /* 'fitNGaussians3D_mexCode_F:74' minIndxZ = min(index(:,3)); */
    i = g_index->size[0];
    g_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)g_index, i, (int32_T)sizeof(real_T));
    emxFree_real_T(&f_index);
    ia = b_index->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        g_index->data[i] = b_index->data[i + (b_index->size[0] << 1)];
    }
    c_emxInit_real_T(&h_index, 1, TRUE);
    minIndxZ = b_min(g_index);
    /* 'fitNGaussians3D_mexCode_F:75' maxIndxZ = max(index(:,3)); */
    i = h_index->size[0];
    h_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)h_index, i, (int32_T)sizeof(real_T));
    emxFree_real_T(&g_index);
    ia = b_index->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        h_index->data[i] = b_index->data[i + (b_index->size[0] << 1)];
    }
    emxInit_real_T(&psfIntegX, 2, TRUE);
    maxIndxZ = b_max(h_index);
    /* determine the contribution of each PSF (assuming amplitude 1) to a */
    /* pixel based on its x-coordinate (needed to calculate F & J) */
    /* 'fitNGaussians3D_mexCode_F:79' psfIntegX = zeros(maxIndxX-minIndxX+1,numPSF); */
    i = psfIntegX->size[0] * psfIntegX->size[1];
    psfIntegX->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxX - minIndxX) + 1.0, &emlrtDCI, &emlrtContextGlobal);
    psfIntegX->size[1] = (int32_T)emlrtIntegerCheckR2011a(numPSF, &b_emlrtDCI, &emlrtContextGlobal);
    emxEnsureCapacity((emxArray__common *)psfIntegX, i, (int32_T)sizeof(real_T));
    emxFree_real_T(&h_index);
    ia = (int32_T)emlrtIntegerCheckR2011a((maxIndxX - minIndxX) + 1.0, &e_emlrtDCI, &emlrtContextGlobal) * (int32_T)emlrtIntegerCheckR2011a(numPSF, &f_emlrtDCI, &emlrtContextGlobal) - 1;
    for (i = 0; i <= ia; i++) {
        psfIntegX->data[i] = 0.0;
    }
    /* 'fitNGaussians3D_mexCode_F:80' for i=1:numPSF */
    b_i = 1U;
    b_emxInit_real_T(&temp, 3, TRUE);
    emxInit_real_T(&temp2, 2, TRUE);
    emxInit_real_T(&y, 2, TRUE);
    c_emxInit_real_T(&b_y, 1, TRUE);
    while ((real_T)b_i <= numPSF) {
        /* 'fitNGaussians3D_mexCode_F:81' temp = GaussListND_mexCode((minIndxX:maxIndxX)',... */
        /* 'fitNGaussians3D_mexCode_F:82'         psfSigma(1),psfPos(i,1)); */
        if (muDoubleScalarIsNaN(minIndxX) || muDoubleScalarIsNaN(maxIndxX)) {
            iy = 1;
            anew = rtNaN;
            apnd = maxIndxX;
        } else if (maxIndxX < minIndxX) {
            iy = 0;
            anew = minIndxX;
            apnd = maxIndxX;
        } else if (muDoubleScalarIsInf(minIndxX) || muDoubleScalarIsInf(maxIndxX)) {
            iy = 1;
            anew = rtNaN;
            apnd = maxIndxX;
        } else {
            anew = minIndxX;
            ndbl = muDoubleScalarFloor((maxIndxX - minIndxX) + 0.5);
            apnd = minIndxX + ndbl;
            cdiff = apnd - maxIndxX;
            absa = muDoubleScalarAbs(minIndxX);
            absb = muDoubleScalarAbs(maxIndxX);
            if (absa > absb) {
                absb = absa;
            }
            if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * absb) {
                ndbl++;
                apnd = maxIndxX;
            } else if (cdiff > 0.0) {
                apnd = minIndxX + (ndbl - 1.0);
            } else {
                ndbl++;
            }
            if (ndbl >= 0.0) {
                iy = (int32_T)ndbl;
            } else {
                iy = 0;
            }
        }
        i = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = iy;
        emxEnsureCapacity((emxArray__common *)y, i, (int32_T)sizeof(real_T));
        if (iy > 0) {
            y->data[0] = anew;
            if (iy > 1) {
                y->data[iy - 1] = apnd;
                ib = iy - 1;
                i = ib;
                nx = (int32_T)((uint32_T)i >> 1);
                ia = nx - 1;
                for (k = 1; k <= ia; k++) {
                    y->data[k] = anew + (real_T)k;
                    y->data[(iy - k) - 1] = apnd - (real_T)k;
                }
                if (nx << 1 == ib) {
                    y->data[nx] = (anew + apnd) / 2.0;
                } else {
                    y->data[nx] = anew + (real_T)nx;
                    y->data[nx + 1] = apnd - (real_T)nx;
                }
            }
        }
        i = b_y->size[0];
        b_y->size[0] = y->size[1];
        emxEnsureCapacity((emxArray__common *)b_y, i, (int32_T)sizeof(real_T));
        ia = y->size[1] - 1;
        for (i = 0; i <= ia; i++) {
            b_y->data[i] = y->data[i];
        }
        GaussListND_mexCode(b_y, psfSigma[0], d_x0->data[(int32_T)b_i - 1], temp);
        /* 'fitNGaussians3D_mexCode_F:84' temp2 = squeeze(temp); */
        squeeze(temp, temp2);
        /* clear temp */
        /* 'fitNGaussians3D_mexCode_F:87' psfIntegX(:,i) = temp2(:,1); */
        i = (int32_T)b_i - 1;
        ia = temp2->size[0] - 1;
        for (nx = 0; nx <= ia; nx++) {
            psfIntegX->data[nx + psfIntegX->size[0] * i] = temp2->data[nx];
        }
        /* clear temp2 */
        b_i++;
    }
    emxFree_real_T(&b_y);
    emxInit_real_T(&psfIntegY, 2, TRUE);
    /* determine the contribution of each PSF (assuming amplitude 1) to a */
    /* pixel based on its y-coordinate (needed to calculate F & J) */
    /* 'fitNGaussians3D_mexCode_F:93' psfIntegY = zeros(maxIndxY-minIndxY+1,numPSF); */
    i = psfIntegY->size[0] * psfIntegY->size[1];
    psfIntegY->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxY - minIndxY) + 1.0, &c_emlrtDCI, &emlrtContextGlobal);
    psfIntegY->size[1] = (int32_T)numPSF;
    emxEnsureCapacity((emxArray__common *)psfIntegY, i, (int32_T)sizeof(real_T));
    ia = (int32_T)emlrtIntegerCheckR2011a((maxIndxY - minIndxY) + 1.0, &g_emlrtDCI, &emlrtContextGlobal) * (int32_T)numPSF - 1;
    for (i = 0; i <= ia; i++) {
        psfIntegY->data[i] = 0.0;
    }
    /* 'fitNGaussians3D_mexCode_F:94' for i=1:numPSF */
    b_i = 1U;
    c_emxInit_real_T(&c_y, 1, TRUE);
    while (b_i <= (uint32_T)numPSF) {
        /* 'fitNGaussians3D_mexCode_F:95' temp = GaussListND_mexCode((minIndxY:maxIndxY)',... */
        /* 'fitNGaussians3D_mexCode_F:96'         psfSigma(1),psfPos(i,2)); */
        if (muDoubleScalarIsNaN(minIndxY) || muDoubleScalarIsNaN(maxIndxY)) {
            iy = 1;
            anew = rtNaN;
            apnd = maxIndxY;
        } else if (maxIndxY < minIndxY) {
            iy = 0;
            anew = minIndxY;
            apnd = maxIndxY;
        } else if (muDoubleScalarIsInf(minIndxY) || muDoubleScalarIsInf(maxIndxY)) {
            iy = 1;
            anew = rtNaN;
            apnd = maxIndxY;
        } else {
            anew = minIndxY;
            ndbl = muDoubleScalarFloor((maxIndxY - minIndxY) + 0.5);
            apnd = minIndxY + ndbl;
            cdiff = apnd - maxIndxY;
            absa = muDoubleScalarAbs(minIndxY);
            absb = muDoubleScalarAbs(maxIndxY);
            if (absa > absb) {
                absb = absa;
            }
            if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * absb) {
                ndbl++;
                apnd = maxIndxY;
            } else if (cdiff > 0.0) {
                apnd = minIndxY + (ndbl - 1.0);
            } else {
                ndbl++;
            }
            if (ndbl >= 0.0) {
                iy = (int32_T)ndbl;
            } else {
                iy = 0;
            }
        }
        i = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = iy;
        emxEnsureCapacity((emxArray__common *)y, i, (int32_T)sizeof(real_T));
        if (iy > 0) {
            y->data[0] = anew;
            if (iy > 1) {
                y->data[iy - 1] = apnd;
                ib = iy - 1;
                i = ib;
                nx = (int32_T)((uint32_T)i >> 1);
                ia = nx - 1;
                for (k = 1; k <= ia; k++) {
                    y->data[k] = anew + (real_T)k;
                    y->data[(iy - k) - 1] = apnd - (real_T)k;
                }
                if (nx << 1 == ib) {
                    y->data[nx] = (anew + apnd) / 2.0;
                } else {
                    y->data[nx] = anew + (real_T)nx;
                    y->data[nx + 1] = apnd - (real_T)nx;
                }
            }
        }
        i = c_y->size[0];
        c_y->size[0] = y->size[1];
        emxEnsureCapacity((emxArray__common *)c_y, i, (int32_T)sizeof(real_T));
        ia = y->size[1] - 1;
        for (i = 0; i <= ia; i++) {
            c_y->data[i] = y->data[i];
        }
        GaussListND_mexCode(c_y, psfSigma[0], d_x0->data[((int32_T)b_i + d_x0->size[0]) - 1], temp);
        /* 'fitNGaussians3D_mexCode_F:98' temp2 = squeeze(temp); */
        squeeze(temp, temp2);
        /* clear temp */
        /* 'fitNGaussians3D_mexCode_F:101' psfIntegY(:,i) = temp2(:,1); */
        i = (int32_T)b_i - 1;
        ia = temp2->size[0] - 1;
        for (nx = 0; nx <= ia; nx++) {
            psfIntegY->data[nx + psfIntegY->size[0] * i] = temp2->data[nx];
        }
        /* clear temp2 */
        b_i++;
    }
    emxFree_real_T(&c_y);
    emxInit_real_T(&psfIntegZ, 2, TRUE);
    /* determine the contribution of each PSF (assuming amplitude 1) to a */
    /* pixel based on its z-coordinate (needed to calculate F & J) */
    /* 'fitNGaussians3D_mexCode_F:107' psfIntegZ = zeros(maxIndxZ-minIndxZ+1,numPSF); */
    i = psfIntegZ->size[0] * psfIntegZ->size[1];
    psfIntegZ->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxZ - minIndxZ) + 1.0, &d_emlrtDCI, &emlrtContextGlobal);
    psfIntegZ->size[1] = (int32_T)numPSF;
    emxEnsureCapacity((emxArray__common *)psfIntegZ, i, (int32_T)sizeof(real_T));
    ia = (int32_T)emlrtIntegerCheckR2011a((maxIndxZ - minIndxZ) + 1.0, &h_emlrtDCI, &emlrtContextGlobal) * (int32_T)numPSF - 1;
    for (i = 0; i <= ia; i++) {
        psfIntegZ->data[i] = 0.0;
    }
    /* 'fitNGaussians3D_mexCode_F:108' for i=1:numPSF */
    b_i = 1U;
    c_emxInit_real_T(&d_y, 1, TRUE);
    while (b_i <= (uint32_T)numPSF) {
        /* 'fitNGaussians3D_mexCode_F:109' temp = GaussListND_mexCode((minIndxZ:maxIndxZ)',... */
        /* 'fitNGaussians3D_mexCode_F:110'         psfSigma(2),psfPos(i,3)); */
        if (muDoubleScalarIsNaN(minIndxZ) || muDoubleScalarIsNaN(maxIndxZ)) {
            iy = 1;
            anew = rtNaN;
            apnd = maxIndxZ;
        } else if (maxIndxZ < minIndxZ) {
            iy = 0;
            anew = minIndxZ;
            apnd = maxIndxZ;
        } else if (muDoubleScalarIsInf(minIndxZ) || muDoubleScalarIsInf(maxIndxZ)) {
            iy = 1;
            anew = rtNaN;
            apnd = maxIndxZ;
        } else {
            anew = minIndxZ;
            ndbl = muDoubleScalarFloor((maxIndxZ - minIndxZ) + 0.5);
            apnd = minIndxZ + ndbl;
            cdiff = apnd - maxIndxZ;
            absa = muDoubleScalarAbs(minIndxZ);
            absb = muDoubleScalarAbs(maxIndxZ);
            if (absa > absb) {
                absb = absa;
            }
            if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * absb) {
                ndbl++;
                apnd = maxIndxZ;
            } else if (cdiff > 0.0) {
                apnd = minIndxZ + (ndbl - 1.0);
            } else {
                ndbl++;
            }
            if (ndbl >= 0.0) {
                iy = (int32_T)ndbl;
            } else {
                iy = 0;
            }
        }
        i = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = iy;
        emxEnsureCapacity((emxArray__common *)y, i, (int32_T)sizeof(real_T));
        if (iy > 0) {
            y->data[0] = anew;
            if (iy > 1) {
                y->data[iy - 1] = apnd;
                ib = iy - 1;
                i = ib;
                nx = (int32_T)((uint32_T)i >> 1);
                ia = nx - 1;
                for (k = 1; k <= ia; k++) {
                    y->data[k] = anew + (real_T)k;
                    y->data[(iy - k) - 1] = apnd - (real_T)k;
                }
                if (nx << 1 == ib) {
                    y->data[nx] = (anew + apnd) / 2.0;
                } else {
                    y->data[nx] = anew + (real_T)nx;
                    y->data[nx + 1] = apnd - (real_T)nx;
                }
            }
        }
        i = d_y->size[0];
        d_y->size[0] = y->size[1];
        emxEnsureCapacity((emxArray__common *)d_y, i, (int32_T)sizeof(real_T));
        ia = y->size[1] - 1;
        for (i = 0; i <= ia; i++) {
            d_y->data[i] = y->data[i];
        }
        GaussListND_mexCode(d_y, psfSigma[1], d_x0->data[((int32_T)b_i + (d_x0->size[0] << 1)) - 1], temp);
        /* 'fitNGaussians3D_mexCode_F:113' temp2 = squeeze(temp); */
        squeeze(temp, temp2);
        /* clear temp */
        /* 'fitNGaussians3D_mexCode_F:116' psfIntegZ(:,i) = temp2(:,1); */
        i = (int32_T)b_i - 1;
        ia = temp2->size[0] - 1;
        for (nx = 0; nx <= ia; nx++) {
            psfIntegZ->data[nx + psfIntegZ->size[0] * i] = temp2->data[nx];
        }
        /* clear temp2 */
        b_i++;
    }
    emxFree_real_T(&d_y);
    emxFree_real_T(&temp);
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
    if (image->size[0] == 0) {
        nx = 0;
    } else if (image->size[0] > 1) {
        nx = image->size[0];
    } else {
        nx = 1;
    }
    emxInit_int32_T(&r0, 1, TRUE);
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
    sz[1] = nx;
    i = r0->size[0];
    r0->size[0] = d_x0->size[0];
    emxEnsureCapacity((emxArray__common *)r0, i, (int32_T)sizeof(int32_T));
    ia = d_x0->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        r0->data[i] = 1 + i;
    }
    iv0[0] = r0->size[0];
    iv0[1] = 1;
    emxFree_int32_T(&r0);
    for (i = 0; i < 2; i++) {
        iv1[i] = iv0[i] * sz[i];
    }
    emxInit_real_T(&x, 2, TRUE);
    i = x->size[0] * x->size[1];
    x->size[0] = iv1[0];
    x->size[1] = iv1[1];
    emxEnsureCapacity((emxArray__common *)x, i, (int32_T)sizeof(real_T));
    if ((x->size[0] == 0) || (x->size[1] == 0)) {
    } else {
        emxInit_int32_T(&r1, 1, TRUE);
        i = r1->size[0];
        r1->size[0] = d_x0->size[0];
        emxEnsureCapacity((emxArray__common *)r1, i, (int32_T)sizeof(int32_T));
        ia = d_x0->size[0] - 1;
        for (i = 0; i <= ia; i++) {
            r1->data[i] = 1 + i;
        }
        nx = r1->size[0];
        ib = 0;
        iy = 1;
        emxFree_int32_T(&r1);
        while (iy <= sz[1]) {
            ia = 0;
            for (k = 1; k <= nx; k++) {
                x->data[ib] = d_x0->data[ia + d_x0->size[0] * 3];
                ia++;
                ib++;
            }
            iy++;
        }
    }
    emxFree_real_T(&d_x0);
    emxInit_real_T(&r2, 2, TRUE);
    i = r2->size[0] * r2->size[1];
    r2->size[0] = psfIntegX->size[1];
    r2->size[1] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)r2, i, (int32_T)sizeof(real_T));
    ia = b_index->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        ib = psfIntegX->size[1] - 1;
        for (nx = 0; nx <= ib; nx++) {
            r2->data[nx + r2->size[0] * i] = psfIntegX->data[((int32_T)emlrtIntegerCheckR2011a((b_index->data[i] - minIndxX) + 1.0, &i_emlrtDCI, &emlrtContextGlobal) + psfIntegX->size[0] * nx) - 1];
        }
    }
    i = psfIntegX->size[0] * psfIntegX->size[1];
    psfIntegX->size[0] = psfIntegY->size[1];
    psfIntegX->size[1] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)psfIntegX, i, (int32_T)sizeof(real_T));
    ia = b_index->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        ib = psfIntegY->size[1] - 1;
        for (nx = 0; nx <= ib; nx++) {
            psfIntegX->data[nx + psfIntegX->size[0] * i] = psfIntegY->data[((int32_T)emlrtIntegerCheckR2011a((b_index->data[i + b_index->size[0]] - minIndxY) + 1.0, &j_emlrtDCI, &emlrtContextGlobal) + psfIntegY->size[0] * nx) - 1];
        }
    }
    emxFree_real_T(&psfIntegY);
    i = temp2->size[0] * temp2->size[1];
    temp2->size[0] = psfIntegZ->size[1];
    temp2->size[1] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)temp2, i, (int32_T)sizeof(real_T));
    ia = b_index->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        ib = psfIntegZ->size[1] - 1;
        for (nx = 0; nx <= ib; nx++) {
            temp2->data[nx + temp2->size[0] * i] = psfIntegZ->data[((int32_T)emlrtIntegerCheckR2011a((b_index->data[i + (b_index->size[0] << 1)] - minIndxZ) + 1.0, &k_emlrtDCI, &emlrtContextGlobal) + psfIntegZ->size[0] * nx) - 1];
        }
    }
    emxFree_real_T(&psfIntegZ);
    i = x->size[0] * x->size[1];
    emxEnsureCapacity((emxArray__common *)x, i, (int32_T)sizeof(real_T));
    nx = x->size[0];
    ib = x->size[1];
    ia = nx * ib - 1;
    for (i = 0; i <= ia; i++) {
        x->data[i] = x->data[i] * r2->data[i] * psfIntegX->data[i] * temp2->data[i];
    }
    emxFree_real_T(&r2);
    emxFree_real_T(&temp2);
    emxFree_real_T(&psfIntegX);
    for (i = 0; i < 2; i++) {
        b_sz[i] = (uint32_T)x->size[i];
    }
    b_sz[0] = 1U;
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = (int32_T)b_sz[1];
    emxEnsureCapacity((emxArray__common *)y, i, (int32_T)sizeof(real_T));
    if ((x->size[0] == 0) || (x->size[1] == 0)) {
        i = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity((emxArray__common *)y, i, (int32_T)sizeof(real_T));
        ia = y->size[1] - 1;
        for (i = 0; i <= ia; i++) {
            y->data[y->size[0] * i] = 0.0;
        }
    } else {
        nx = x->size[0];
        ib = x->size[1];
        ia = -1;
        iy = -1;
        for (i = 1; i <= ib; i++) {
            ia++;
            minIndxX = x->data[ia];
            for (k = 2; k <= nx; k++) {
                ia++;
                minIndxX += x->data[ia];
            }
            iy++;
            y->data[iy] = minIndxX;
        }
    }
    emxFree_real_T(&x);
    c_emxInit_real_T(&r3, 1, TRUE);
    i = r3->size[0];
    r3->size[0] = y->size[1];
    emxEnsureCapacity((emxArray__common *)r3, i, (int32_T)sizeof(real_T));
    ia = y->size[1] - 1;
    for (i = 0; i <= ia; i++) {
        r3->data[i] = y->data[i];
    }
    emxFree_real_T(&y);
    nx = r3->size[0];
    i = F->size[0] * F->size[1];
    F->size[0] = nx;
    emxEnsureCapacity((emxArray__common *)F, i, (int32_T)sizeof(real_T));
    i = F->size[0] * F->size[1];
    F->size[1] = 1;
    emxEnsureCapacity((emxArray__common *)F, i, (int32_T)sizeof(real_T));
    ia = r3->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        F->data[i] = (r3->data[i] + bgAmp) - image->data[i];
    }
    emxFree_real_T(&r3);
    emxInit_boolean_T(&indxPixel, 1, TRUE);
    /* remove pixels with NaN (which means they are out of the cropped image */
    /* area) */
    /* 'fitNGaussians3D_mexCode_F:165' indxPixel = ~isnan(image); */
    i = indxPixel->size[0];
    indxPixel->size[0] = image->size[0];
    emxEnsureCapacity((emxArray__common *)indxPixel, i, (int32_T)sizeof(boolean_T));
    ia = image->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        indxPixel->data[i] = muDoubleScalarIsNaN(image->data[i]);
    }
    i = indxPixel->size[0];
    emxEnsureCapacity((emxArray__common *)indxPixel, i, (int32_T)sizeof(boolean_T));
    ia = indxPixel->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        indxPixel->data[i] = !indxPixel->data[i];
    }
    /* 'fitNGaussians3D_mexCode_F:167' F = F(indxPixel); */
    iy = indxPixel->size[0];
    k = 0;
    for (i = 1; i <= iy; i++) {
        if (indxPixel->data[i - 1]) {
            k++;
        }
    }
    emxInit_int32_T(&r4, 1, TRUE);
    i = r4->size[0];
    r4->size[0] = k;
    emxEnsureCapacity((emxArray__common *)r4, i, (int32_T)sizeof(int32_T));
    ib = 0;
    for (i = 1; i <= iy; i++) {
        if (indxPixel->data[i - 1]) {
            r4->data[ib] = i;
            ib++;
        }
    }
    emxFree_boolean_T(&indxPixel);
    c_emxInit_real_T(&b_F, 1, TRUE);
    i = b_F->size[0];
    b_F->size[0] = r4->size[0];
    emxEnsureCapacity((emxArray__common *)b_F, i, (int32_T)sizeof(real_T));
    ia = r4->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        b_F->data[i] = F->data[r4->data[i] - 1];
    }
    c_emxInit_real_T(&c_F, 1, TRUE);
    i = c_F->size[0];
    c_F->size[0] = r4->size[0];
    emxEnsureCapacity((emxArray__common *)c_F, i, (int32_T)sizeof(real_T));
    ia = r4->size[0] - 1;
    for (i = 0; i <= ia; i++) {
        c_F->data[i] = F->data[r4->data[i] - 1];
    }
    emxFree_int32_T(&r4);
    d_F[0] = b_F->size[0];
    d_F[1] = 1;
    i = F->size[0] * F->size[1];
    F->size[0] = d_F[0];
    F->size[1] = d_F[1];
    emxEnsureCapacity((emxArray__common *)F, i, (int32_T)sizeof(real_T));
    emxFree_real_T(&b_F);
    ia = d_F[1] - 1;
    for (i = 0; i <= ia; i++) {
        ib = d_F[0] - 1;
        for (nx = 0; nx <= ib; nx++) {
            e_F = *c_F;
            e_F.size = (int32_T *)&d_F;
            e_F.numDimensions = 1;
            F->data[nx + F->size[0] * i] = e_F.data[nx + e_F.size[0] * i];
        }
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
    emlrtHeapReferenceStackLeaveFcn();
}
/* End of code generation (fitNGaussians3D_mexCode_F.c) */
