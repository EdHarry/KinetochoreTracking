/*
 * fitNGaussians3D_mexCode_J.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode_J'
 *
 * C source code generated on: Mon May 21 19:42:21 2012
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
#include "length.h"
#include "max.h"
#include "min.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
static emlrtDCInfo emlrtDCI = { 80, 19, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo b_emlrtDCI = { 80, 39, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo c_emlrtDCI = { 94, 19, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo d_emlrtDCI = { 108, 19, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo e_emlrtDCI = { 123, 19, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo f_emlrtDCI = { 131, 19, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo g_emlrtDCI = { 139, 19, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo h_emlrtDCI = { 179, 10, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo i_emlrtDCI = { 80, 19, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo j_emlrtDCI = { 80, 39, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo k_emlrtDCI = { 94, 19, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo l_emlrtDCI = { 108, 19, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo m_emlrtDCI = { 123, 19, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo n_emlrtDCI = { 131, 19, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo o_emlrtDCI = { 139, 19, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo p_emlrtDCI = { 179, 10, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo q_emlrtDCI = { 180, 60, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo r_emlrtDCI = { 181, 41, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };
static emlrtDCInfo s_emlrtDCI = { 181, 64, "fitNGaussians3D_mexCode_J", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/fitNGaussians3D_mexCode_J.m", 1 };

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
 */
void fitNGaussians3D_mexCode_J(emxArray_real_T *x0, const emxArray_real_T *image, const emxArray_real_T *b_index, const real_T psfSigma[2], emxArray_real_T *J)
{
    int32_T i0;
    emxArray_real_T *b_x0;
    int32_T i1;
    int32_T loop_ub;
    real_T numPSF;
    int32_T nx;
    int32_T sz[2];
    emxArray_real_T *c_x0;
    real_T y;
    int32_T k;
    emxArray_real_T *d_x0;
    emxArray_real_T *c_index;
    emxArray_real_T *d_index;
    real_T minIndxX;
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
    uint32_T i;
    emxArray_real_T *temp;
    emxArray_real_T *psfValueX;
    emxArray_real_T *b_y;
    emxArray_real_T *c_y;
    int32_T n;
    real_T apnd;
    real_T ndbl;
    real_T cdiff;
    real_T absa;
    real_T absb;
    int32_T nm1;
    emxArray_real_T *psfIntegY;
    emxArray_real_T *d_y;
    emxArray_real_T *psfIntegZ;
    emxArray_real_T *e_y;
    emxArray_real_T *relIndxX;
    emxArray_real_T *b_relIndxX;
    emxArray_real_T *c_relIndxX;
    emxArray_real_T *f_y;
    real_T numPixel;
    emxArray_real_T *psfValueY;
    emxArray_real_T *d_relIndxX;
    emxArray_real_T *e_relIndxX;
    emxArray_real_T *g_y;
    emxArray_real_T *psfValueZ;
    emxArray_real_T *f_relIndxX;
    emxArray_real_T *g_relIndxX;
    emxArray_real_T *h_y;
    emxArray_real_T *relIndxY;
    emxArray_real_T *relIndxZ;
    emxArray_real_T *e_x0;
    emxArray_real_T *r0;
    emxArray_real_T *f_x0;
    emxArray_real_T *g_x0;
    emxArray_boolean_T *r1;
    emxArray_boolean_T *r2;
    emxArray_int32_T *r3;
    emxArray_real_T *b_J;
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
    /* coder.extrinsic('sparse'); */
    /* 'fitNGaussians3D_mexCode_J:30' J = []; */
    /* % Input */
    /*  %check whether correct number of input arguments was used */
    /*  if nargin ~= 4 */
    /*      disp('--fitNGaussians2D: Incorrect number of input arguments!'); */
    /*      return */
    /*  end */
    /* check whether correct number of input arguments was used */
    /* 'fitNGaussians3D_mexCode_J:40' if nargin ~= 4 */
    /* % Calculating F & J */
    /* extract background intensity from x0 and remove from vector */
    /*  bgAmp = x0(end); */
    /* 'fitNGaussians3D_mexCode_J:49' x0 = x0(1:end-1); */
    i0 = x0->size[0] - 1;
    if (1 > i0) {
        i0 = 0;
    }
    c_emxInit_real_T(&b_x0, 1, TRUE);
    i1 = b_x0->size[0];
    b_x0->size[0] = i0;
    emxEnsureCapacity((emxArray__common *)b_x0, i1, (int32_T)sizeof(real_T));
    loop_ub = i0 - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        b_x0->data[i0] = x0->data[i0];
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
    /* 'fitNGaussians3D_mexCode_J:53' numPSF = length(x0)/4; */
    numPSF = length(x0) / 4.0;
    /*  %reshape 3nx1 vector x0 into nx3 matrix */
    /*  x0 = reshape(x0,3,numPSF); */
    /*  x0 = x0'; */
    /* reshape 4nx1 vector x0 into nx4 matrix */
    /* 'fitNGaussians3D_mexCode_J:60' x0 = reshape(x0,4,numPSF); */
    nx = x0->size[0];
    for (i0 = 0; i0 < 2; i0++) {
        sz[i0] = 0;
    }
    emxInit_real_T(&c_x0, 2, TRUE);
    sz[0] = 4;
    y = numPSF;
    y = y < 0.0 ? muDoubleScalarCeil(y - 0.5) : muDoubleScalarFloor(y + 0.5);
    if (y < 2.147483648E+9) {
        if (y >= -2.147483648E+9) {
            i0 = (int32_T)y;
        } else {
            i0 = MIN_int32_T;
        }
    } else if (y >= 2.147483648E+9) {
        i0 = MAX_int32_T;
    } else {
        i0 = 0;
    }
    sz[1] = i0;
    i0 = c_x0->size[0] * c_x0->size[1];
    c_x0->size[0] = 4;
    c_x0->size[1] = sz[1];
    emxEnsureCapacity((emxArray__common *)c_x0, i0, (int32_T)sizeof(real_T));
    for (k = 0; k + 1 <= nx; k++) {
        c_x0->data[k] = x0->data[k];
    }
    emxInit_real_T(&d_x0, 2, TRUE);
    /* 'fitNGaussians3D_mexCode_J:61' x0 = x0'; */
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
    c_emxInit_real_T(&c_index, 1, TRUE);
    /* extract PSF center positions and amplitudes */
    /*  psfPos = x0(:,1:2); */
    /*  psfAmp = x0(:,3); */
    /* 'fitNGaussians3D_mexCode_J:67' psfPos = x0(:,1:3); */
    /* 'fitNGaussians3D_mexCode_J:68' psfAmp = x0(:,4); */
    /* find minimum and maximum pixel indices */
    /* 'fitNGaussians3D_mexCode_J:71' minIndxX = min(index(:,1)); */
    i0 = c_index->size[0];
    c_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)c_index, i0, (int32_T)sizeof(real_T));
    loop_ub = b_index->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        c_index->data[i0] = b_index->data[i0];
    }
    c_emxInit_real_T(&d_index, 1, TRUE);
    minIndxX = b_min(c_index);
    /* 'fitNGaussians3D_mexCode_J:72' maxIndxX = max(index(:,1)); */
    i0 = d_index->size[0];
    d_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)d_index, i0, (int32_T)sizeof(real_T));
    emxFree_real_T(&c_index);
    loop_ub = b_index->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        d_index->data[i0] = b_index->data[i0];
    }
    c_emxInit_real_T(&e_index, 1, TRUE);
    maxIndxX = b_max(d_index);
    /* 'fitNGaussians3D_mexCode_J:73' minIndxY = min(index(:,2)); */
    i0 = e_index->size[0];
    e_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)e_index, i0, (int32_T)sizeof(real_T));
    emxFree_real_T(&d_index);
    loop_ub = b_index->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        e_index->data[i0] = b_index->data[i0 + b_index->size[0]];
    }
    c_emxInit_real_T(&f_index, 1, TRUE);
    minIndxY = b_min(e_index);
    /* 'fitNGaussians3D_mexCode_J:74' maxIndxY = max(index(:,2)); */
    i0 = f_index->size[0];
    f_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)f_index, i0, (int32_T)sizeof(real_T));
    emxFree_real_T(&e_index);
    loop_ub = b_index->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        f_index->data[i0] = b_index->data[i0 + b_index->size[0]];
    }
    c_emxInit_real_T(&g_index, 1, TRUE);
    maxIndxY = b_max(f_index);
    /* 'fitNGaussians3D_mexCode_J:75' minIndxZ = min(index(:,3)); */
    i0 = g_index->size[0];
    g_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)g_index, i0, (int32_T)sizeof(real_T));
    emxFree_real_T(&f_index);
    loop_ub = b_index->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        g_index->data[i0] = b_index->data[i0 + (b_index->size[0] << 1)];
    }
    c_emxInit_real_T(&h_index, 1, TRUE);
    minIndxZ = b_min(g_index);
    /* 'fitNGaussians3D_mexCode_J:76' maxIndxZ = max(index(:,3)); */
    i0 = h_index->size[0];
    h_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)h_index, i0, (int32_T)sizeof(real_T));
    emxFree_real_T(&g_index);
    loop_ub = b_index->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        h_index->data[i0] = b_index->data[i0 + (b_index->size[0] << 1)];
    }
    emxInit_real_T(&psfIntegX, 2, TRUE);
    maxIndxZ = b_max(h_index);
    /* determine the contribution of each PSF (assuming amplitude 1) to a */
    /* pixel based on its x-coordinate (needed to calculate F & J) */
    /* 'fitNGaussians3D_mexCode_J:80' psfIntegX = zeros(maxIndxX-minIndxX+1,numPSF); */
    i0 = psfIntegX->size[0] * psfIntegX->size[1];
    psfIntegX->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxX - minIndxX) + 1.0, &emlrtDCI, &emlrtContextGlobal);
    psfIntegX->size[1] = (int32_T)emlrtIntegerCheckR2011a(numPSF, &b_emlrtDCI, &emlrtContextGlobal);
    emxEnsureCapacity((emxArray__common *)psfIntegX, i0, (int32_T)sizeof(real_T));
    emxFree_real_T(&h_index);
    loop_ub = (int32_T)emlrtIntegerCheckR2011a((maxIndxX - minIndxX) + 1.0, &i_emlrtDCI, &emlrtContextGlobal) * (int32_T)emlrtIntegerCheckR2011a(numPSF, &j_emlrtDCI, &emlrtContextGlobal) - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        psfIntegX->data[i0] = 0.0;
    }
    /* 'fitNGaussians3D_mexCode_J:81' for i=1:numPSF */
    i = 1U;
    b_emxInit_real_T(&temp, 3, TRUE);
    emxInit_real_T(&psfValueX, 2, TRUE);
    emxInit_real_T(&b_y, 2, TRUE);
    c_emxInit_real_T(&c_y, 1, TRUE);
    while ((real_T)i <= numPSF) {
        /* 'fitNGaussians3D_mexCode_J:82' temp = GaussListND_mexCode((minIndxX:maxIndxX)',... */
        /* 'fitNGaussians3D_mexCode_J:83'         psfSigma(1),psfPos(i,1)); */
        if (muDoubleScalarIsNaN(minIndxX) || muDoubleScalarIsNaN(maxIndxX)) {
            n = 1;
            y = rtNaN;
            apnd = maxIndxX;
        } else if (maxIndxX < minIndxX) {
            n = 0;
            y = minIndxX;
            apnd = maxIndxX;
        } else if (muDoubleScalarIsInf(minIndxX) || muDoubleScalarIsInf(maxIndxX)) {
            n = 1;
            y = rtNaN;
            apnd = maxIndxX;
        } else {
            y = minIndxX;
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
                n = (int32_T)ndbl;
            } else {
                n = 0;
            }
        }
        i0 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = 1;
        b_y->size[1] = n;
        emxEnsureCapacity((emxArray__common *)b_y, i0, (int32_T)sizeof(real_T));
        if (n > 0) {
            b_y->data[0] = y;
            if (n > 1) {
                b_y->data[n - 1] = apnd;
                nm1 = n - 1;
                i0 = nm1;
                nx = (int32_T)((uint32_T)i0 >> 1);
                loop_ub = nx - 1;
                for (k = 1; k <= loop_ub; k++) {
                    b_y->data[k] = y + (real_T)k;
                    b_y->data[(n - k) - 1] = apnd - (real_T)k;
                }
                if (nx << 1 == nm1) {
                    b_y->data[nx] = (y + apnd) / 2.0;
                } else {
                    b_y->data[nx] = y + (real_T)nx;
                    b_y->data[nx + 1] = apnd - (real_T)nx;
                }
            }
        }
        i0 = c_y->size[0];
        c_y->size[0] = b_y->size[1];
        emxEnsureCapacity((emxArray__common *)c_y, i0, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i0 = 0; i0 <= loop_ub; i0++) {
            c_y->data[i0] = b_y->data[i0];
        }
        GaussListND_mexCode(c_y, psfSigma[0], d_x0->data[(int32_T)i - 1], temp);
        /* 'fitNGaussians3D_mexCode_J:85' temp2 = squeeze(temp); */
        squeeze(temp, psfValueX);
        /* clear temp */
        /* 'fitNGaussians3D_mexCode_J:88' psfIntegX(:,i) = temp2(:,1); */
        i0 = (int32_T)i - 1;
        loop_ub = psfValueX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            psfIntegX->data[i1 + psfIntegX->size[0] * i0] = psfValueX->data[i1];
        }
        /* clear temp2 */
        i++;
    }
    emxFree_real_T(&c_y);
    emxInit_real_T(&psfIntegY, 2, TRUE);
    /* determine the contribution of each PSF (assuming amplitude 1) to a */
    /* pixel based on its y-coordinate (needed to calculate F & J) */
    /* 'fitNGaussians3D_mexCode_J:94' psfIntegY = zeros(maxIndxY-minIndxY+1,numPSF); */
    i0 = psfIntegY->size[0] * psfIntegY->size[1];
    psfIntegY->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxY - minIndxY) + 1.0, &c_emlrtDCI, &emlrtContextGlobal);
    psfIntegY->size[1] = (int32_T)numPSF;
    emxEnsureCapacity((emxArray__common *)psfIntegY, i0, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)emlrtIntegerCheckR2011a((maxIndxY - minIndxY) + 1.0, &k_emlrtDCI, &emlrtContextGlobal) * (int32_T)numPSF - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        psfIntegY->data[i0] = 0.0;
    }
    /* 'fitNGaussians3D_mexCode_J:95' for i=1:numPSF */
    i = 1U;
    c_emxInit_real_T(&d_y, 1, TRUE);
    while (i <= (uint32_T)numPSF) {
        /* 'fitNGaussians3D_mexCode_J:96' temp = GaussListND_mexCode((minIndxY:maxIndxY)',... */
        /* 'fitNGaussians3D_mexCode_J:97'         psfSigma(1),psfPos(i,2)); */
        if (muDoubleScalarIsNaN(minIndxY) || muDoubleScalarIsNaN(maxIndxY)) {
            n = 1;
            y = rtNaN;
            apnd = maxIndxY;
        } else if (maxIndxY < minIndxY) {
            n = 0;
            y = minIndxY;
            apnd = maxIndxY;
        } else if (muDoubleScalarIsInf(minIndxY) || muDoubleScalarIsInf(maxIndxY)) {
            n = 1;
            y = rtNaN;
            apnd = maxIndxY;
        } else {
            y = minIndxY;
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
                n = (int32_T)ndbl;
            } else {
                n = 0;
            }
        }
        i0 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = 1;
        b_y->size[1] = n;
        emxEnsureCapacity((emxArray__common *)b_y, i0, (int32_T)sizeof(real_T));
        if (n > 0) {
            b_y->data[0] = y;
            if (n > 1) {
                b_y->data[n - 1] = apnd;
                nm1 = n - 1;
                i0 = nm1;
                nx = (int32_T)((uint32_T)i0 >> 1);
                loop_ub = nx - 1;
                for (k = 1; k <= loop_ub; k++) {
                    b_y->data[k] = y + (real_T)k;
                    b_y->data[(n - k) - 1] = apnd - (real_T)k;
                }
                if (nx << 1 == nm1) {
                    b_y->data[nx] = (y + apnd) / 2.0;
                } else {
                    b_y->data[nx] = y + (real_T)nx;
                    b_y->data[nx + 1] = apnd - (real_T)nx;
                }
            }
        }
        i0 = d_y->size[0];
        d_y->size[0] = b_y->size[1];
        emxEnsureCapacity((emxArray__common *)d_y, i0, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i0 = 0; i0 <= loop_ub; i0++) {
            d_y->data[i0] = b_y->data[i0];
        }
        GaussListND_mexCode(d_y, psfSigma[0], d_x0->data[((int32_T)i + d_x0->size[0]) - 1], temp);
        /* 'fitNGaussians3D_mexCode_J:99' temp2 = squeeze(temp); */
        squeeze(temp, psfValueX);
        /* clear temp */
        /* 'fitNGaussians3D_mexCode_J:102' psfIntegY(:,i) = temp2(:,1); */
        i0 = (int32_T)i - 1;
        loop_ub = psfValueX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            psfIntegY->data[i1 + psfIntegY->size[0] * i0] = psfValueX->data[i1];
        }
        /* clear temp2 */
        i++;
    }
    emxFree_real_T(&d_y);
    emxInit_real_T(&psfIntegZ, 2, TRUE);
    /* determine the contribution of each PSF (assuming amplitude 1) to a */
    /* pixel based on its z-coordinate (needed to calculate F & J) */
    /* 'fitNGaussians3D_mexCode_J:108' psfIntegZ = zeros(maxIndxZ-minIndxZ+1,numPSF); */
    i0 = psfIntegZ->size[0] * psfIntegZ->size[1];
    psfIntegZ->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxZ - minIndxZ) + 1.0, &d_emlrtDCI, &emlrtContextGlobal);
    psfIntegZ->size[1] = (int32_T)numPSF;
    emxEnsureCapacity((emxArray__common *)psfIntegZ, i0, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)emlrtIntegerCheckR2011a((maxIndxZ - minIndxZ) + 1.0, &l_emlrtDCI, &emlrtContextGlobal) * (int32_T)numPSF - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        psfIntegZ->data[i0] = 0.0;
    }
    /* 'fitNGaussians3D_mexCode_J:109' for i=1:numPSF */
    i = 1U;
    c_emxInit_real_T(&e_y, 1, TRUE);
    while (i <= (uint32_T)numPSF) {
        /* 'fitNGaussians3D_mexCode_J:110' temp = GaussListND_mexCode((minIndxZ:maxIndxZ)',... */
        /* 'fitNGaussians3D_mexCode_J:111'         psfSigma(2),psfPos(i,3)); */
        if (muDoubleScalarIsNaN(minIndxZ) || muDoubleScalarIsNaN(maxIndxZ)) {
            n = 1;
            y = rtNaN;
            apnd = maxIndxZ;
        } else if (maxIndxZ < minIndxZ) {
            n = 0;
            y = minIndxZ;
            apnd = maxIndxZ;
        } else if (muDoubleScalarIsInf(minIndxZ) || muDoubleScalarIsInf(maxIndxZ)) {
            n = 1;
            y = rtNaN;
            apnd = maxIndxZ;
        } else {
            y = minIndxZ;
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
                n = (int32_T)ndbl;
            } else {
                n = 0;
            }
        }
        i0 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = 1;
        b_y->size[1] = n;
        emxEnsureCapacity((emxArray__common *)b_y, i0, (int32_T)sizeof(real_T));
        if (n > 0) {
            b_y->data[0] = y;
            if (n > 1) {
                b_y->data[n - 1] = apnd;
                nm1 = n - 1;
                i0 = nm1;
                nx = (int32_T)((uint32_T)i0 >> 1);
                loop_ub = nx - 1;
                for (k = 1; k <= loop_ub; k++) {
                    b_y->data[k] = y + (real_T)k;
                    b_y->data[(n - k) - 1] = apnd - (real_T)k;
                }
                if (nx << 1 == nm1) {
                    b_y->data[nx] = (y + apnd) / 2.0;
                } else {
                    b_y->data[nx] = y + (real_T)nx;
                    b_y->data[nx + 1] = apnd - (real_T)nx;
                }
            }
        }
        i0 = e_y->size[0];
        e_y->size[0] = b_y->size[1];
        emxEnsureCapacity((emxArray__common *)e_y, i0, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i0 = 0; i0 <= loop_ub; i0++) {
            e_y->data[i0] = b_y->data[i0];
        }
        GaussListND_mexCode(e_y, psfSigma[1], d_x0->data[((int32_T)i + (d_x0->size[0] << 1)) - 1], temp);
        /* 'fitNGaussians3D_mexCode_J:114' temp2 = squeeze(temp); */
        squeeze(temp, psfValueX);
        /* clear temp */
        /* 'fitNGaussians3D_mexCode_J:117' psfIntegZ(:,i) = temp2(:,1); */
        i0 = (int32_T)i - 1;
        loop_ub = psfValueX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            psfIntegZ->data[i1 + psfIntegZ->size[0] * i0] = psfValueX->data[i1];
        }
        /* clear temp2 */
        i++;
    }
    emxFree_real_T(&e_y);
    emxFree_real_T(&temp);
    /* calculate the value of each PSF (assuming amplitude 1) at the */
    /* x-coordinates of the corners of all pixels (needed to calculate J) */
    /* 'fitNGaussians3D_mexCode_J:123' psfValueX = zeros(maxIndxX-minIndxX+2,numPSF); */
    i0 = psfValueX->size[0] * psfValueX->size[1];
    psfValueX->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxX - minIndxX) + 2.0, &e_emlrtDCI, &emlrtContextGlobal);
    psfValueX->size[1] = (int32_T)numPSF;
    emxEnsureCapacity((emxArray__common *)psfValueX, i0, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)emlrtIntegerCheckR2011a((maxIndxX - minIndxX) + 2.0, &m_emlrtDCI, &emlrtContextGlobal) * (int32_T)numPSF - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        psfValueX->data[i0] = 0.0;
    }
    /* 'fitNGaussians3D_mexCode_J:124' for i=1:numPSF */
    i = 1U;
    c_emxInit_real_T(&relIndxX, 1, TRUE);
    c_emxInit_real_T(&b_relIndxX, 1, TRUE);
    c_emxInit_real_T(&c_relIndxX, 1, TRUE);
    c_emxInit_real_T(&f_y, 1, TRUE);
    while (i <= (uint32_T)numPSF) {
        /* 'fitNGaussians3D_mexCode_J:125' psfValueX(:,i) = exp(-((minIndxX-0.5:maxIndxX+0.5)'... */
        /* 'fitNGaussians3D_mexCode_J:126'         -psfPos(i,1)).^2/2/psfSigma(1)^2); */
        i0 = (int32_T)i - 1;
        y = minIndxX - 0.5;
        numPixel = maxIndxX + 0.5;
        if (muDoubleScalarIsNaN(y) || muDoubleScalarIsNaN(numPixel)) {
            n = 1;
            y = rtNaN;
            apnd = numPixel;
        } else if (numPixel < y) {
            n = 0;
            apnd = numPixel;
        } else if (muDoubleScalarIsInf(y) || muDoubleScalarIsInf(numPixel)) {
            n = 1;
            y = rtNaN;
            apnd = numPixel;
        } else {
            ndbl = muDoubleScalarFloor((numPixel - y) + 0.5);
            apnd = y + ndbl;
            cdiff = apnd - numPixel;
            absa = muDoubleScalarAbs(y);
            absb = muDoubleScalarAbs(numPixel);
            if (absa > absb) {
                absb = absa;
            }
            if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * absb) {
                ndbl++;
                apnd = numPixel;
            } else if (cdiff > 0.0) {
                apnd = y + (ndbl - 1.0);
            } else {
                ndbl++;
            }
            if (ndbl >= 0.0) {
                n = (int32_T)ndbl;
            } else {
                n = 0;
            }
        }
        i1 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = 1;
        b_y->size[1] = n;
        emxEnsureCapacity((emxArray__common *)b_y, i1, (int32_T)sizeof(real_T));
        if (n > 0) {
            b_y->data[0] = y;
            if (n > 1) {
                b_y->data[n - 1] = apnd;
                nm1 = n - 1;
                i1 = nm1;
                nx = (int32_T)((uint32_T)i1 >> 1);
                loop_ub = nx - 1;
                for (k = 1; k <= loop_ub; k++) {
                    b_y->data[k] = y + (real_T)k;
                    b_y->data[(n - k) - 1] = apnd - (real_T)k;
                }
                if (nx << 1 == nm1) {
                    b_y->data[nx] = (y + apnd) / 2.0;
                } else {
                    b_y->data[nx] = y + (real_T)nx;
                    b_y->data[nx + 1] = apnd - (real_T)nx;
                }
            }
        }
        y = d_x0->data[(int32_T)i - 1];
        i1 = f_y->size[0];
        f_y->size[0] = b_y->size[1];
        emxEnsureCapacity((emxArray__common *)f_y, i1, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            f_y->data[i1] = b_y->data[i1] - y;
        }
        power(f_y, relIndxX);
        i1 = relIndxX->size[0];
        emxEnsureCapacity((emxArray__common *)relIndxX, i1, (int32_T)sizeof(real_T));
        loop_ub = relIndxX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            relIndxX->data[i1] = -relIndxX->data[i1];
        }
        i1 = c_relIndxX->size[0];
        c_relIndxX->size[0] = relIndxX->size[0];
        emxEnsureCapacity((emxArray__common *)c_relIndxX, i1, (int32_T)sizeof(real_T));
        loop_ub = relIndxX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            c_relIndxX->data[i1] = relIndxX->data[i1];
        }
        rdivide(c_relIndxX, 2.0, relIndxX);
        i1 = b_relIndxX->size[0];
        b_relIndxX->size[0] = relIndxX->size[0];
        emxEnsureCapacity((emxArray__common *)b_relIndxX, i1, (int32_T)sizeof(real_T));
        loop_ub = relIndxX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            b_relIndxX->data[i1] = relIndxX->data[i1];
        }
        rdivide(b_relIndxX, muDoubleScalarPower(psfSigma[0], 2.0), relIndxX);
        b_exp(relIndxX);
        loop_ub = relIndxX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            psfValueX->data[i1 + psfValueX->size[0] * i0] = relIndxX->data[i1];
        }
        i++;
    }
    emxFree_real_T(&f_y);
    emxFree_real_T(&c_relIndxX);
    emxFree_real_T(&b_relIndxX);
    emxInit_real_T(&psfValueY, 2, TRUE);
    /* calculate the value of each PSF (assuming amplitude 1) at the */
    /* y-coordinates of the corners of all pixels (needed to calculate J) */
    /* 'fitNGaussians3D_mexCode_J:131' psfValueY = zeros(maxIndxY-minIndxY+2,numPSF); */
    i0 = psfValueY->size[0] * psfValueY->size[1];
    psfValueY->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxY - minIndxY) + 2.0, &f_emlrtDCI, &emlrtContextGlobal);
    psfValueY->size[1] = (int32_T)numPSF;
    emxEnsureCapacity((emxArray__common *)psfValueY, i0, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)emlrtIntegerCheckR2011a((maxIndxY - minIndxY) + 2.0, &n_emlrtDCI, &emlrtContextGlobal) * (int32_T)numPSF - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        psfValueY->data[i0] = 0.0;
    }
    /* 'fitNGaussians3D_mexCode_J:132' for i=1:numPSF */
    i = 1U;
    c_emxInit_real_T(&d_relIndxX, 1, TRUE);
    c_emxInit_real_T(&e_relIndxX, 1, TRUE);
    c_emxInit_real_T(&g_y, 1, TRUE);
    while (i <= (uint32_T)numPSF) {
        /* 'fitNGaussians3D_mexCode_J:133' psfValueY(:,i) = exp(-((minIndxY-0.5:maxIndxY+0.5)'... */
        /* 'fitNGaussians3D_mexCode_J:134'         -psfPos(i,2)).^2/2/psfSigma(1)^2); */
        i0 = (int32_T)i - 1;
        y = minIndxY - 0.5;
        numPixel = maxIndxY + 0.5;
        if (muDoubleScalarIsNaN(y) || muDoubleScalarIsNaN(numPixel)) {
            n = 1;
            y = rtNaN;
            apnd = numPixel;
        } else if (numPixel < y) {
            n = 0;
            apnd = numPixel;
        } else if (muDoubleScalarIsInf(y) || muDoubleScalarIsInf(numPixel)) {
            n = 1;
            y = rtNaN;
            apnd = numPixel;
        } else {
            ndbl = muDoubleScalarFloor((numPixel - y) + 0.5);
            apnd = y + ndbl;
            cdiff = apnd - numPixel;
            absa = muDoubleScalarAbs(y);
            absb = muDoubleScalarAbs(numPixel);
            if (absa > absb) {
                absb = absa;
            }
            if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * absb) {
                ndbl++;
                apnd = numPixel;
            } else if (cdiff > 0.0) {
                apnd = y + (ndbl - 1.0);
            } else {
                ndbl++;
            }
            if (ndbl >= 0.0) {
                n = (int32_T)ndbl;
            } else {
                n = 0;
            }
        }
        i1 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = 1;
        b_y->size[1] = n;
        emxEnsureCapacity((emxArray__common *)b_y, i1, (int32_T)sizeof(real_T));
        if (n > 0) {
            b_y->data[0] = y;
            if (n > 1) {
                b_y->data[n - 1] = apnd;
                nm1 = n - 1;
                i1 = nm1;
                nx = (int32_T)((uint32_T)i1 >> 1);
                loop_ub = nx - 1;
                for (k = 1; k <= loop_ub; k++) {
                    b_y->data[k] = y + (real_T)k;
                    b_y->data[(n - k) - 1] = apnd - (real_T)k;
                }
                if (nx << 1 == nm1) {
                    b_y->data[nx] = (y + apnd) / 2.0;
                } else {
                    b_y->data[nx] = y + (real_T)nx;
                    b_y->data[nx + 1] = apnd - (real_T)nx;
                }
            }
        }
        y = d_x0->data[((int32_T)i + d_x0->size[0]) - 1];
        i1 = g_y->size[0];
        g_y->size[0] = b_y->size[1];
        emxEnsureCapacity((emxArray__common *)g_y, i1, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            g_y->data[i1] = b_y->data[i1] - y;
        }
        power(g_y, relIndxX);
        i1 = relIndxX->size[0];
        emxEnsureCapacity((emxArray__common *)relIndxX, i1, (int32_T)sizeof(real_T));
        loop_ub = relIndxX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            relIndxX->data[i1] = -relIndxX->data[i1];
        }
        i1 = e_relIndxX->size[0];
        e_relIndxX->size[0] = relIndxX->size[0];
        emxEnsureCapacity((emxArray__common *)e_relIndxX, i1, (int32_T)sizeof(real_T));
        loop_ub = relIndxX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            e_relIndxX->data[i1] = relIndxX->data[i1];
        }
        rdivide(e_relIndxX, 2.0, relIndxX);
        i1 = d_relIndxX->size[0];
        d_relIndxX->size[0] = relIndxX->size[0];
        emxEnsureCapacity((emxArray__common *)d_relIndxX, i1, (int32_T)sizeof(real_T));
        loop_ub = relIndxX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            d_relIndxX->data[i1] = relIndxX->data[i1];
        }
        rdivide(d_relIndxX, muDoubleScalarPower(psfSigma[0], 2.0), relIndxX);
        b_exp(relIndxX);
        loop_ub = relIndxX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            psfValueY->data[i1 + psfValueY->size[0] * i0] = relIndxX->data[i1];
        }
        i++;
    }
    emxFree_real_T(&g_y);
    emxFree_real_T(&e_relIndxX);
    emxFree_real_T(&d_relIndxX);
    emxInit_real_T(&psfValueZ, 2, TRUE);
    /* calculate the value of each PSF (assuming amplitude 1) at the */
    /* z-coordinates of the corners of all pixels (needed to calculate J) */
    /* 'fitNGaussians3D_mexCode_J:139' psfValueZ = zeros(maxIndxZ-minIndxZ+2,numPSF); */
    i0 = psfValueZ->size[0] * psfValueZ->size[1];
    psfValueZ->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxZ - minIndxZ) + 2.0, &g_emlrtDCI, &emlrtContextGlobal);
    psfValueZ->size[1] = (int32_T)numPSF;
    emxEnsureCapacity((emxArray__common *)psfValueZ, i0, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)emlrtIntegerCheckR2011a((maxIndxZ - minIndxZ) + 2.0, &o_emlrtDCI, &emlrtContextGlobal) * (int32_T)numPSF - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        psfValueZ->data[i0] = 0.0;
    }
    /* 'fitNGaussians3D_mexCode_J:140' for i=1:numPSF */
    i = 1U;
    c_emxInit_real_T(&f_relIndxX, 1, TRUE);
    c_emxInit_real_T(&g_relIndxX, 1, TRUE);
    c_emxInit_real_T(&h_y, 1, TRUE);
    while (i <= (uint32_T)numPSF) {
        /* 'fitNGaussians3D_mexCode_J:141' psfValueZ(:,i) = exp(-((minIndxZ-0.5:maxIndxZ+0.5)'... */
        /* 'fitNGaussians3D_mexCode_J:142'         -psfPos(i,3)).^2/2/psfSigma(2)^2); */
        i0 = (int32_T)i - 1;
        y = minIndxZ - 0.5;
        numPixel = maxIndxZ + 0.5;
        if (muDoubleScalarIsNaN(y) || muDoubleScalarIsNaN(numPixel)) {
            n = 1;
            y = rtNaN;
            apnd = numPixel;
        } else if (numPixel < y) {
            n = 0;
            apnd = numPixel;
        } else if (muDoubleScalarIsInf(y) || muDoubleScalarIsInf(numPixel)) {
            n = 1;
            y = rtNaN;
            apnd = numPixel;
        } else {
            ndbl = muDoubleScalarFloor((numPixel - y) + 0.5);
            apnd = y + ndbl;
            cdiff = apnd - numPixel;
            absa = muDoubleScalarAbs(y);
            absb = muDoubleScalarAbs(numPixel);
            if (absa > absb) {
                absb = absa;
            }
            if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * absb) {
                ndbl++;
                apnd = numPixel;
            } else if (cdiff > 0.0) {
                apnd = y + (ndbl - 1.0);
            } else {
                ndbl++;
            }
            if (ndbl >= 0.0) {
                n = (int32_T)ndbl;
            } else {
                n = 0;
            }
        }
        i1 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = 1;
        b_y->size[1] = n;
        emxEnsureCapacity((emxArray__common *)b_y, i1, (int32_T)sizeof(real_T));
        if (n > 0) {
            b_y->data[0] = y;
            if (n > 1) {
                b_y->data[n - 1] = apnd;
                nm1 = n - 1;
                i1 = nm1;
                nx = (int32_T)((uint32_T)i1 >> 1);
                loop_ub = nx - 1;
                for (k = 1; k <= loop_ub; k++) {
                    b_y->data[k] = y + (real_T)k;
                    b_y->data[(n - k) - 1] = apnd - (real_T)k;
                }
                if (nx << 1 == nm1) {
                    b_y->data[nx] = (y + apnd) / 2.0;
                } else {
                    b_y->data[nx] = y + (real_T)nx;
                    b_y->data[nx + 1] = apnd - (real_T)nx;
                }
            }
        }
        y = d_x0->data[((int32_T)i + (d_x0->size[0] << 1)) - 1];
        i1 = h_y->size[0];
        h_y->size[0] = b_y->size[1];
        emxEnsureCapacity((emxArray__common *)h_y, i1, (int32_T)sizeof(real_T));
        loop_ub = b_y->size[1] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            h_y->data[i1] = b_y->data[i1] - y;
        }
        power(h_y, relIndxX);
        i1 = relIndxX->size[0];
        emxEnsureCapacity((emxArray__common *)relIndxX, i1, (int32_T)sizeof(real_T));
        loop_ub = relIndxX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            relIndxX->data[i1] = -relIndxX->data[i1];
        }
        i1 = g_relIndxX->size[0];
        g_relIndxX->size[0] = relIndxX->size[0];
        emxEnsureCapacity((emxArray__common *)g_relIndxX, i1, (int32_T)sizeof(real_T));
        loop_ub = relIndxX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            g_relIndxX->data[i1] = relIndxX->data[i1];
        }
        rdivide(g_relIndxX, 2.0, relIndxX);
        i1 = f_relIndxX->size[0];
        f_relIndxX->size[0] = relIndxX->size[0];
        emxEnsureCapacity((emxArray__common *)f_relIndxX, i1, (int32_T)sizeof(real_T));
        loop_ub = relIndxX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            f_relIndxX->data[i1] = relIndxX->data[i1];
        }
        rdivide(f_relIndxX, muDoubleScalarPower(psfSigma[1], 2.0), relIndxX);
        b_exp(relIndxX);
        loop_ub = relIndxX->size[0] - 1;
        for (i1 = 0; i1 <= loop_ub; i1++) {
            psfValueZ->data[i1 + psfValueZ->size[0] * i0] = relIndxX->data[i1];
        }
        i++;
    }
    emxFree_real_T(&h_y);
    emxFree_real_T(&g_relIndxX);
    emxFree_real_T(&f_relIndxX);
    emxFree_real_T(&b_y);
    /* get number of pixels in image */
    /* 'fitNGaussians3D_mexCode_J:146' numPixel = length(image); */
    numPixel = length(image);
    /*  %get xy-indices relative to minimum */
    /*  relIndxX = index(:,1) - minIndxX + 1; */
    /*  relIndxY = index(:,2) - minIndxY + 1; */
    /* get xy-indices relative to minimum */
    /* 'fitNGaussians3D_mexCode_J:152' relIndxX = index(:,1) - minIndxX + 1; */
    i0 = relIndxX->size[0];
    relIndxX->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)relIndxX, i0, (int32_T)sizeof(real_T));
    loop_ub = b_index->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        relIndxX->data[i0] = (b_index->data[i0] - minIndxX) + 1.0;
    }
    c_emxInit_real_T(&relIndxY, 1, TRUE);
    /* 'fitNGaussians3D_mexCode_J:153' relIndxY = index(:,2) - minIndxY + 1; */
    i0 = relIndxY->size[0];
    relIndxY->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)relIndxY, i0, (int32_T)sizeof(real_T));
    loop_ub = b_index->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        relIndxY->data[i0] = (b_index->data[i0 + b_index->size[0]] - minIndxY) + 1.0;
    }
    c_emxInit_real_T(&relIndxZ, 1, TRUE);
    /* 'fitNGaussians3D_mexCode_J:154' relIndxZ = index(:,3) - minIndxZ + 1; */
    i0 = relIndxZ->size[0];
    relIndxZ->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)relIndxZ, i0, (int32_T)sizeof(real_T));
    loop_ub = b_index->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        relIndxZ->data[i0] = (b_index->data[i0 + (b_index->size[0] << 1)] - minIndxZ) + 1.0;
    }
    /*  %calculate the value of F at all pixels */
    /*  F = (sum(repmat(psfAmp,1,numPixel).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)',1))' ... */
    /*      + repmat(bgAmp,numPixel,1) - image; */
    /* calculate the value of F at all pixels */
    /*  F = (sum(repmat(psfAmp,1,numPixel).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)'.*psfIntegZ(relIndxZ,:)',1))' ... */
    /*      + repmat(bgAmp,numPixel,1) - image; */
    /* remove pixels with NaN (which means they are out of the cropped image */
    /* area) */
    /* 'fitNGaussians3D_mexCode_J:166' indxPixel = ~isnan(image); */
    /*  F = F(indxPixel); */
    /*  if nargout > 1 */
    /* calculate the derivative at all pixels */
    /*      J = ones(numPixel,3*numPSF+1); %(last column for background amplitude) */
    /*      J(:,1:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueX(relIndxX,:)-... */
    /*          psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:); %w.r.t. x */
    /*      J(:,2:3:3*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueY(relIndxY,:)-... */
    /*          psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:); %w.r.t. y */
    /*      J(:,3:3:3*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); %w.r.t. amp */
    /* 'fitNGaussians3D_mexCode_J:179' J = ones(numPixel,4*numPSF+1); */
    y = 4.0 * numPSF;
    i0 = J->size[0] * J->size[1];
    J->size[0] = (int32_T)emlrtIntegerCheckR2011a(numPixel, &h_emlrtDCI, &emlrtContextGlobal);
    J->size[1] = (int32_T)(y + 1.0);
    emxEnsureCapacity((emxArray__common *)J, i0, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)emlrtIntegerCheckR2011a(numPixel, &p_emlrtDCI, &emlrtContextGlobal) * (int32_T)(y + 1.0) - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        J->data[i0] = 1.0;
    }
    /* (last column for background amplitude) */
    /* 'fitNGaussians3D_mexCode_J:180' J(:,1:4:4*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueX(relIndxX,:)-... */
    /* 'fitNGaussians3D_mexCode_J:181'     psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:).*psfIntegZ(relIndxZ,:); */
    loop_ub = relIndxX->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        emlrtIntegerCheckR2011a(relIndxX->data[i0], &q_emlrtDCI, &emlrtContextGlobal);
    }
    emxInit_real_T(&e_x0, 2, TRUE);
    i0 = e_x0->size[0] * e_x0->size[1];
    e_x0->size[0] = 1;
    e_x0->size[1] = d_x0->size[0];
    emxEnsureCapacity((emxArray__common *)e_x0, i0, (int32_T)sizeof(real_T));
    loop_ub = d_x0->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        e_x0->data[e_x0->size[0] * i0] = d_x0->data[i0 + d_x0->size[0] * 3];
    }
    emxInit_real_T(&r0, 2, TRUE);
    repmat(e_x0, numPixel, r0);
    emxFree_real_T(&e_x0);
    loop_ub = relIndxY->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        emlrtIntegerCheckR2011a(relIndxY->data[i0], &r_emlrtDCI, &emlrtContextGlobal);
    }
    loop_ub = relIndxZ->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        emlrtIntegerCheckR2011a(relIndxZ->data[i0], &s_emlrtDCI, &emlrtContextGlobal);
    }
    if (1.0 > 4.0 * numPSF) {
        i0 = 1;
    } else {
        i0 = 4;
    }
    loop_ub = r0->size[1] - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
        k = r0->size[0] - 1;
        for (nx = 0; nx <= k; nx++) {
            J->data[nx + J->size[0] * (i0 * i1)] = r0->data[nx + r0->size[0] * i1] * (psfValueX->data[((int32_T)relIndxX->data[nx] + psfValueX->size[0] * i1) - 1] - psfValueX->data[((int32_T)(relIndxX->data[nx] + 1.0) + psfValueX->size[0] * i1) - 1]) * psfIntegY->data[((int32_T)relIndxY->data[nx] + psfIntegY->size[0] * i1) - 1] * psfIntegZ->data[((int32_T)relIndxZ->data[nx] + psfIntegZ->size[0] * i1) - 1];
        }
    }
    emxFree_real_T(&psfValueX);
    emxInit_real_T(&f_x0, 2, TRUE);
    /* w.r.t. x */
    /* 'fitNGaussians3D_mexCode_J:182' J(:,2:4:4*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueY(relIndxY,:)-... */
    /* 'fitNGaussians3D_mexCode_J:183'     psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:).*psfIntegZ(relIndxZ,:); */
    i0 = f_x0->size[0] * f_x0->size[1];
    f_x0->size[0] = 1;
    f_x0->size[1] = d_x0->size[0];
    emxEnsureCapacity((emxArray__common *)f_x0, i0, (int32_T)sizeof(real_T));
    loop_ub = d_x0->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        f_x0->data[f_x0->size[0] * i0] = d_x0->data[i0 + d_x0->size[0] * 3];
    }
    repmat(f_x0, numPixel, r0);
    emxFree_real_T(&f_x0);
    if (2.0 > 4.0 * numPSF) {
        i0 = 0;
        i1 = 1;
    } else {
        i0 = 1;
        i1 = 4;
    }
    loop_ub = r0->size[1] - 1;
    for (nx = 0; nx <= loop_ub; nx++) {
        k = r0->size[0] - 1;
        for (n = 0; n <= k; n++) {
            J->data[n + J->size[0] * (i0 + i1 * nx)] = r0->data[n + r0->size[0] * nx] * (psfValueY->data[((int32_T)relIndxY->data[n] + psfValueY->size[0] * nx) - 1] - psfValueY->data[((int32_T)(relIndxY->data[n] + 1.0) + psfValueY->size[0] * nx) - 1]) * psfIntegX->data[((int32_T)relIndxX->data[n] + psfIntegX->size[0] * nx) - 1] * psfIntegZ->data[((int32_T)relIndxZ->data[n] + psfIntegZ->size[0] * nx) - 1];
        }
    }
    emxFree_real_T(&psfValueY);
    emxInit_real_T(&g_x0, 2, TRUE);
    /* w.r.t. y */
    /* 'fitNGaussians3D_mexCode_J:184' J(:,3:4:4*numPSF) = repmat(psfAmp',numPixel,1).*(psfValueZ(relIndxZ,:)-... */
    /* 'fitNGaussians3D_mexCode_J:185'     psfValueZ(relIndxZ+1,:)).*psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); */
    i0 = g_x0->size[0] * g_x0->size[1];
    g_x0->size[0] = 1;
    g_x0->size[1] = d_x0->size[0];
    emxEnsureCapacity((emxArray__common *)g_x0, i0, (int32_T)sizeof(real_T));
    loop_ub = d_x0->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        g_x0->data[g_x0->size[0] * i0] = d_x0->data[i0 + d_x0->size[0] * 3];
    }
    emxFree_real_T(&d_x0);
    repmat(g_x0, numPixel, r0);
    emxFree_real_T(&g_x0);
    if (3.0 > 4.0 * numPSF) {
        i0 = 0;
        i1 = 1;
    } else {
        i0 = 2;
        i1 = 4;
    }
    loop_ub = r0->size[1] - 1;
    for (nx = 0; nx <= loop_ub; nx++) {
        k = r0->size[0] - 1;
        for (n = 0; n <= k; n++) {
            J->data[n + J->size[0] * (i0 + i1 * nx)] = r0->data[n + r0->size[0] * nx] * (psfValueZ->data[((int32_T)relIndxZ->data[n] + psfValueZ->size[0] * nx) - 1] - psfValueZ->data[((int32_T)(relIndxZ->data[n] + 1.0) + psfValueZ->size[0] * nx) - 1]) * psfIntegX->data[((int32_T)relIndxX->data[n] + psfIntegX->size[0] * nx) - 1] * psfIntegY->data[((int32_T)relIndxY->data[n] + psfIntegY->size[0] * nx) - 1];
        }
    }
    emxFree_real_T(&r0);
    emxFree_real_T(&psfValueZ);
    /* w.r.t. z */
    /* 'fitNGaussians3D_mexCode_J:186' J(:,4:4:4*numPSF) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:).*psfIntegZ(relIndxZ,:); */
    if (4.0 > 4.0 * numPSF) {
        i0 = 0;
        i1 = 1;
    } else {
        i0 = 3;
        i1 = 4;
    }
    loop_ub = psfIntegX->size[1] - 1;
    for (nx = 0; nx <= loop_ub; nx++) {
        k = relIndxX->size[0] - 1;
        for (n = 0; n <= k; n++) {
            J->data[n + J->size[0] * (i0 + i1 * nx)] = psfIntegX->data[((int32_T)relIndxX->data[n] + psfIntegX->size[0] * nx) - 1] * psfIntegY->data[((int32_T)relIndxY->data[n] + psfIntegY->size[0] * nx) - 1] * psfIntegZ->data[((int32_T)relIndxZ->data[n] + psfIntegZ->size[0] * nx) - 1];
        }
    }
    emxFree_real_T(&relIndxZ);
    emxFree_real_T(&relIndxY);
    emxFree_real_T(&relIndxX);
    emxFree_real_T(&psfIntegZ);
    emxFree_real_T(&psfIntegY);
    emxFree_real_T(&psfIntegX);
    emxInit_boolean_T(&r1, 1, TRUE);
    /* w.r.t. amp */
    /* remove pixels with NaN (which means they are out of the cropped image */
    /* area) */
    /* 'fitNGaussians3D_mexCode_J:190' J = J(indxPixel,:); */
    i0 = r1->size[0];
    r1->size[0] = image->size[0];
    emxEnsureCapacity((emxArray__common *)r1, i0, (int32_T)sizeof(boolean_T));
    loop_ub = image->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        r1->data[i0] = muDoubleScalarIsNaN(image->data[i0]);
    }
    emxInit_boolean_T(&r2, 1, TRUE);
    i0 = r2->size[0];
    r2->size[0] = r1->size[0];
    emxEnsureCapacity((emxArray__common *)r2, i0, (int32_T)sizeof(boolean_T));
    loop_ub = r1->size[0] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        r2->data[i0] = !r1->data[i0];
    }
    emxFree_boolean_T(&r1);
    emxInit_int32_T(&r3, 1, TRUE);
    emxInit_real_T(&b_J, 2, TRUE);
    eml_li_find(r2, r3);
    nx = J->size[1];
    i0 = b_J->size[0] * b_J->size[1];
    b_J->size[0] = r3->size[0];
    b_J->size[1] = nx;
    emxEnsureCapacity((emxArray__common *)b_J, i0, (int32_T)sizeof(real_T));
    emxFree_boolean_T(&r2);
    loop_ub = nx - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        k = r3->size[0] - 1;
        for (i1 = 0; i1 <= k; i1++) {
            b_J->data[i1 + b_J->size[0] * i0] = J->data[(r3->data[i1] + J->size[0] * i0) - 1];
        }
    }
    emxFree_int32_T(&r3);
    i0 = J->size[0] * J->size[1];
    J->size[0] = b_J->size[0];
    J->size[1] = b_J->size[1];
    emxEnsureCapacity((emxArray__common *)J, i0, (int32_T)sizeof(real_T));
    loop_ub = b_J->size[1] - 1;
    for (i0 = 0; i0 <= loop_ub; i0++) {
        k = b_J->size[0] - 1;
        for (i1 = 0; i1 <= k; i1++) {
            J->data[i1 + J->size[0] * i0] = b_J->data[i1 + b_J->size[0] * i0];
        }
    }
    emxFree_real_T(&b_J);
    /* J = sparse(J); */
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
/* End of code generation (fitNGaussians3D_mexCode_J.c) */
