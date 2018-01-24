/*
 * objFcn.c
 *
 * Code generation for function 'objFcn'
 *
 * C source code generated on: Fri May 25 21:48:50 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "objFcn.h"
#include "squeeze.h"
#include "GaussListND_mexCode.h"
#include "objFcn_emxutil.h"
#include "exp.h"
#include "rdivide.h"
#include "power.h"
#include "repmat.h"
#include "sum.h"
#include "max.h"
#include "min.h"
#include "mrdivide.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */
static emlrtDCInfo emlrtDCI = { 48, 19, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo b_emlrtDCI = { 61, 19, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo c_emlrtDCI = { 74, 19, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo d_emlrtDCI = { 90, 23, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo e_emlrtDCI = { 98, 23, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo f_emlrtDCI = { 106, 23, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo g_emlrtDCI = { 48, 19, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo h_emlrtDCI = { 61, 19, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo i_emlrtDCI = { 74, 19, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo j_emlrtDCI = { 90, 23, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo k_emlrtDCI = { 98, 23, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo l_emlrtDCI = { 106, 23, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo m_emlrtDCI = { 114, 55, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo n_emlrtDCI = { 115, 45, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo o_emlrtDCI = { 115, 68, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo p_emlrtDCI = { 130, 48, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo q_emlrtDCI = { 130, 72, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo r_emlrtDCI = { 130, 96, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo s_emlrtDCI = { 128, 56, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo t_emlrtDCI = { 128, 80, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };
static emlrtDCInfo u_emlrtDCI = { 128, 104, "objFcn", "/gpfs/home/sysbio/msrfby/MATLAB/MATLAB/newTracking/detection/MMF/mexCoding/nag/objFcn.m", 1 };

/* Function Declarations */
static int32_T _s32_add__(int32_T b, int32_T c);
static int32_T _s32_minus__(int32_T b, int32_T c);
static int32_T _s32_shl_s32_(int32_T b, int32_T c);

/* Function Definitions */

static int32_T _s32_add__(int32_T b, int32_T c)
{
    int32_T a;
    a = b + c;
    if (((a ^ b) & (a ^ c)) < 0) {
        emlrtIntegerOverflowErrorR2008a(0);
    }
    return a;
}

static int32_T _s32_minus__(int32_T b, int32_T c)
{
    int32_T a;
    a = b - c;
    if (((b ^ a) & (b ^ c)) < 0) {
        emlrtIntegerOverflowErrorR2008a(0);
    }
    return a;
}

static int32_T _s32_shl_s32_(int32_T b, int32_T c)
{
    int32_T a;
    a = b << c;
    if ((b != a >> c) || ((uint32_T)c > 31)) {
        emlrtIntegerOverflowErrorR2008a(0);
    }
    return a;
}

/*
 * function [mode, f, fjac, user] = objFcn(mode, m, n, ldfj, needfi, x, fjac, nstate, user)
 */
void objFcn(const int32_T *mode, int32_T m, int32_T n, int32_T ldfj, int32_T needfi, emxArray_real_T *x, emxArray_real_T *fjac, int32_T nstate, const struct_T *user, emxArray_real_T *f)
{
    int32_T i8;
    real_T bgAmp;
    emxArray_real_T *b_x;
    int32_T nx;
    int32_T loop_ub;
    int32_T sz[2];
    emxArray_real_T *c_x;
    int32_T k;
    emxArray_real_T *d_x;
    emxArray_real_T *b_index;
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
    emxArray_real_T *relIndxX;
    real_T maxIndxZ;
    emxArray_real_T *relIndxY;
    emxArray_real_T *relIndxZ;
    emxArray_real_T *psfIntegX;
    int32_T i;
    emxArray_real_T *temp;
    emxArray_real_T *psfValueX;
    emxArray_real_T *y;
    emxArray_real_T *b_y;
    int32_T b_n;
    real_T a;
    real_T apnd;
    real_T ndbl;
    real_T cdiff;
    real_T absa;
    real_T absb;
    int32_T nm1;
    emxArray_real_T *psfIntegY;
    emxArray_real_T *c_y;
    emxArray_real_T *psfIntegZ;
    emxArray_real_T *d_y;
    emxArray_real_T *psfValueY;
    emxArray_real_T *psfValueZ;
    emxArray_real_T *r0;
    emxArray_real_T *r1;
    emxArray_real_T *A;
    emxArray_real_T *r2;
    emxArray_real_T *e_y;
    real_T d;
    emxArray_real_T *r3;
    emxArray_real_T *f_y;
    emxArray_real_T *r4;
    emxArray_real_T *g_y;
    emxArray_real_T *e_x;
    emxArray_real_T *f_x;
    emxArray_real_T *g_x;
    emxArray_real_T *h_x;
    emxArray_real_T *r5;
    emxArray_real_T *i_x;
    emxArray_real_T *r6;
    emlrtHeapReferenceStackEnterFcn();
    /*  defaults */
    /* 'objFcn:4' f = []; */
    i8 = f->size[0] * f->size[1];
    f->size[0] = 0;
    f->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)f, i8, (int32_T)sizeof(real_T));
    /* 'objFcn:5' fjac = []; */
    i8 = fjac->size[0] * fjac->size[1];
    fjac->size[0] = 0;
    fjac->size[1] = 0;
    emxEnsureCapacity((emxArray__common *)fjac, i8, (int32_T)sizeof(real_T));
    /* % core */
    /* extract background intensity from x0 and remove from vector */
    /* 'objFcn:9' bgAmp = x(end); */
    bgAmp = x->data[x->size[0] - 1];
    /* 'objFcn:10' x = x(1:end-1); */
    i8 = x->size[0] - 1;
    if (1 > i8) {
        i8 = 0;
    }
    emxInit_real_T(&b_x, 1, TRUE);
    nx = b_x->size[0];
    b_x->size[0] = i8;
    emxEnsureCapacity((emxArray__common *)b_x, nx, (int32_T)sizeof(real_T));
    loop_ub = i8 - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        b_x->data[i8] = x->data[i8];
    }
    i8 = x->size[0];
    x->size[0] = b_x->size[0];
    emxEnsureCapacity((emxArray__common *)x, i8, (int32_T)sizeof(real_T));
    loop_ub = b_x->size[0] - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        x->data[i8] = b_x->data[i8];
    }
    emxFree_real_T(&b_x);
    /* reshape 4nx1 vector x0 into nx4 matrix */
    /* 'objFcn:13' n = (n-1)/4; */
    n = mrdivide(_s32_minus__(n, 1), 4.0);
    /* 'objFcn:14' x = reshape(x,4,n); */
    nx = x->size[0];
    for (i8 = 0; i8 < 2; i8++) {
        sz[i8] = 0;
    }
    b_emxInit_real_T(&c_x, 2, TRUE);
    sz[0] = 4;
    sz[1] = n;
    i8 = c_x->size[0] * c_x->size[1];
    c_x->size[0] = 4;
    c_x->size[1] = sz[1];
    emxEnsureCapacity((emxArray__common *)c_x, i8, (int32_T)sizeof(real_T));
    for (k = 0; k + 1 <= nx; k++) {
        c_x->data[k] = x->data[k];
    }
    b_emxInit_real_T(&d_x, 2, TRUE);
    /* 'objFcn:15' x = x'; */
    i8 = d_x->size[0] * d_x->size[1];
    d_x->size[0] = c_x->size[1];
    d_x->size[1] = 4;
    emxEnsureCapacity((emxArray__common *)d_x, i8, (int32_T)sizeof(real_T));
    for (i8 = 0; i8 < 4; i8++) {
        loop_ub = c_x->size[1] - 1;
        for (nx = 0; nx <= loop_ub; nx++) {
            d_x->data[nx + d_x->size[0] * i8] = c_x->data[i8 + c_x->size[0] * nx];
        }
    }
    emxFree_real_T(&c_x);
    b_emxInit_real_T(&b_index, 2, TRUE);
    /* extract PSF center positions and amplitudes */
    /* 'objFcn:18' psfPos = x(:,1:3); */
    /* 'objFcn:19' psfAmp = x(:,4); */
    /* % index mod */
    /* 'objFcn:22' index = user.index; */
    i8 = b_index->size[0] * b_index->size[1];
    b_index->size[0] = user->index->size[0];
    b_index->size[1] = 3;
    emxEnsureCapacity((emxArray__common *)b_index, i8, (int32_T)sizeof(real_T));
    loop_ub = user->index->size[0] * user->index->size[1] - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        b_index->data[i8] = user->index->data[i8];
    }
    /* 'objFcn:23' psfSigma = user.psfSigma; */
    /* 'objFcn:25' if needfi > 0 */
    if (needfi > 0) {
        /* 'objFcn:26' index = index(needfi,:); */
        i8 = b_index->size[0] * b_index->size[1];
        b_index->size[0] = 1;
        b_index->size[1] = 3;
        emxEnsureCapacity((emxArray__common *)b_index, i8, (int32_T)sizeof(real_T));
        for (i8 = 0; i8 < 3; i8++) {
            b_index->data[b_index->size[0] * i8] = user->index->data[(needfi + user->index->size[0] * i8) - 1];
        }
    }
    emxInit_real_T(&c_index, 1, TRUE);
    /* % pixel integration calculations */
    /* find minimum and maximum pixel indices */
    /* 'objFcn:31' minIndxX = min(index(:,1)); */
    i8 = c_index->size[0];
    c_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)c_index, i8, (int32_T)sizeof(real_T));
    loop_ub = b_index->size[0] - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        c_index->data[i8] = b_index->data[i8];
    }
    emxInit_real_T(&d_index, 1, TRUE);
    minIndxX = b_min(c_index);
    /* 'objFcn:32' maxIndxX = max(index(:,1)); */
    i8 = d_index->size[0];
    d_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)d_index, i8, (int32_T)sizeof(real_T));
    emxFree_real_T(&c_index);
    loop_ub = b_index->size[0] - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        d_index->data[i8] = b_index->data[i8];
    }
    emxInit_real_T(&e_index, 1, TRUE);
    maxIndxX = b_max(d_index);
    /* 'objFcn:33' minIndxY = min(index(:,2)); */
    i8 = e_index->size[0];
    e_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)e_index, i8, (int32_T)sizeof(real_T));
    emxFree_real_T(&d_index);
    loop_ub = b_index->size[0] - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        e_index->data[i8] = b_index->data[i8 + b_index->size[0]];
    }
    emxInit_real_T(&f_index, 1, TRUE);
    minIndxY = b_min(e_index);
    /* 'objFcn:34' maxIndxY = max(index(:,2)); */
    i8 = f_index->size[0];
    f_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)f_index, i8, (int32_T)sizeof(real_T));
    emxFree_real_T(&e_index);
    loop_ub = b_index->size[0] - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        f_index->data[i8] = b_index->data[i8 + b_index->size[0]];
    }
    emxInit_real_T(&g_index, 1, TRUE);
    maxIndxY = b_max(f_index);
    /* 'objFcn:35' minIndxZ = min(index(:,3)); */
    i8 = g_index->size[0];
    g_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)g_index, i8, (int32_T)sizeof(real_T));
    emxFree_real_T(&f_index);
    loop_ub = b_index->size[0] - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        g_index->data[i8] = b_index->data[i8 + (b_index->size[0] << 1)];
    }
    emxInit_real_T(&h_index, 1, TRUE);
    minIndxZ = b_min(g_index);
    /* 'objFcn:36' maxIndxZ = max(index(:,3)); */
    i8 = h_index->size[0];
    h_index->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)h_index, i8, (int32_T)sizeof(real_T));
    emxFree_real_T(&g_index);
    loop_ub = b_index->size[0] - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        h_index->data[i8] = b_index->data[i8 + (b_index->size[0] << 1)];
    }
    emxInit_real_T(&relIndxX, 1, TRUE);
    maxIndxZ = b_max(h_index);
    /*  get xy-indices relative to minimum */
    /*  relIndxX = index(:,1) - minIndxX + 1; */
    /*  relIndxY = index(:,2) - minIndxY + 1; */
    /*  get xy-indices relative to minimum */
    /* 'objFcn:42' relIndxX = index(:,1) - minIndxX + 1; */
    i8 = relIndxX->size[0];
    relIndxX->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)relIndxX, i8, (int32_T)sizeof(real_T));
    emxFree_real_T(&h_index);
    loop_ub = b_index->size[0] - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        relIndxX->data[i8] = (b_index->data[i8] - minIndxX) + 1.0;
    }
    emxInit_real_T(&relIndxY, 1, TRUE);
    /* 'objFcn:43' relIndxY = index(:,2) - minIndxY + 1; */
    i8 = relIndxY->size[0];
    relIndxY->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)relIndxY, i8, (int32_T)sizeof(real_T));
    loop_ub = b_index->size[0] - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        relIndxY->data[i8] = (b_index->data[i8 + b_index->size[0]] - minIndxY) + 1.0;
    }
    emxInit_real_T(&relIndxZ, 1, TRUE);
    /* 'objFcn:44' relIndxZ = index(:,3) - minIndxZ + 1; */
    i8 = relIndxZ->size[0];
    relIndxZ->size[0] = b_index->size[0];
    emxEnsureCapacity((emxArray__common *)relIndxZ, i8, (int32_T)sizeof(real_T));
    loop_ub = b_index->size[0] - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        relIndxZ->data[i8] = (b_index->data[i8 + (b_index->size[0] << 1)] - minIndxZ) + 1.0;
    }
    emxFree_real_T(&b_index);
    b_emxInit_real_T(&psfIntegX, 2, TRUE);
    /* determine the contribution of each PSF (assuming amplitude 1) to a */
    /* pixel based on its x-coordinate (needed to calculate F & J) */
    /* 'objFcn:48' psfIntegX = zeros(maxIndxX-minIndxX+1,n); */
    i8 = psfIntegX->size[0] * psfIntegX->size[1];
    psfIntegX->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxX - minIndxX) + 1.0, &emlrtDCI, &emlrtContextGlobal);
    psfIntegX->size[1] = n;
    emxEnsureCapacity((emxArray__common *)psfIntegX, i8, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)emlrtIntegerCheckR2011a((maxIndxX - minIndxX) + 1.0, &g_emlrtDCI, &emlrtContextGlobal) * n - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        psfIntegX->data[i8] = 0.0;
    }
    /* 'objFcn:49' for i=1:n */
    i = 0;
    c_emxInit_real_T(&temp, 3, TRUE);
    b_emxInit_real_T(&psfValueX, 2, TRUE);
    b_emxInit_real_T(&y, 2, TRUE);
    emxInit_real_T(&b_y, 1, TRUE);
    while (i + 1 <= n) {
        /* 'objFcn:50' temp = GaussListND_mexCode((minIndxX:maxIndxX)',... */
        /* 'objFcn:51'         psfSigma(1),psfPos(i,1)); */
        if (muDoubleScalarIsNaN(minIndxX) || muDoubleScalarIsNaN(maxIndxX)) {
            b_n = 1;
            a = rtNaN;
            apnd = maxIndxX;
        } else if (maxIndxX < minIndxX) {
            b_n = 0;
            a = minIndxX;
            apnd = maxIndxX;
        } else if (muDoubleScalarIsInf(minIndxX) || muDoubleScalarIsInf(maxIndxX)) {
            b_n = 1;
            a = rtNaN;
            apnd = maxIndxX;
        } else {
            a = minIndxX;
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
                b_n = (int32_T)ndbl;
            } else {
                b_n = 0;
            }
        }
        i8 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = b_n;
        emxEnsureCapacity((emxArray__common *)y, i8, (int32_T)sizeof(real_T));
        if (b_n > 0) {
            y->data[0] = a;
            if (b_n > 1) {
                y->data[b_n - 1] = apnd;
                nm1 = b_n - 1;
                i8 = nm1;
                nx = (int32_T)((uint32_T)i8 >> 1);
                loop_ub = nx - 1;
                for (k = 1; k <= loop_ub; k++) {
                    y->data[k] = a + (real_T)k;
                    y->data[(b_n - k) - 1] = apnd - (real_T)k;
                }
                if (nx << 1 == nm1) {
                    y->data[nx] = (a + apnd) / 2.0;
                } else {
                    y->data[nx] = a + (real_T)nx;
                    y->data[nx + 1] = apnd - (real_T)nx;
                }
            }
        }
        i8 = b_y->size[0];
        b_y->size[0] = y->size[1];
        emxEnsureCapacity((emxArray__common *)b_y, i8, (int32_T)sizeof(real_T));
        loop_ub = y->size[1] - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            b_y->data[i8] = y->data[i8];
        }
        GaussListND_mexCode(b_y, user->psfSigma[0], d_x->data[i], temp);
        /* 'objFcn:53' temp2 = squeeze(temp); */
        squeeze(temp, psfValueX);
        /* clear temp */
        /* 'objFcn:56' psfIntegX(:,i) = temp2(:,1); */
        loop_ub = psfValueX->size[0] - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            psfIntegX->data[i8 + psfIntegX->size[0] * i] = psfValueX->data[i8];
        }
        i++;
    }
    emxFree_real_T(&b_y);
    b_emxInit_real_T(&psfIntegY, 2, TRUE);
    /* determine the contribution of each PSF (assuming amplitude 1) to a */
    /* pixel based on its y-coordinate (needed to calculate F & J) */
    /* 'objFcn:61' psfIntegY = zeros(maxIndxY-minIndxY+1,n); */
    i8 = psfIntegY->size[0] * psfIntegY->size[1];
    psfIntegY->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxY - minIndxY) + 1.0, &b_emlrtDCI, &emlrtContextGlobal);
    psfIntegY->size[1] = n;
    emxEnsureCapacity((emxArray__common *)psfIntegY, i8, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)emlrtIntegerCheckR2011a((maxIndxY - minIndxY) + 1.0, &h_emlrtDCI, &emlrtContextGlobal) * n - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        psfIntegY->data[i8] = 0.0;
    }
    /* 'objFcn:62' for i=1:n */
    i = 0;
    emxInit_real_T(&c_y, 1, TRUE);
    while (i + 1 <= n) {
        /* 'objFcn:63' temp = GaussListND_mexCode((minIndxY:maxIndxY)',... */
        /* 'objFcn:64'         psfSigma(1),psfPos(i,2)); */
        if (muDoubleScalarIsNaN(minIndxY) || muDoubleScalarIsNaN(maxIndxY)) {
            b_n = 1;
            a = rtNaN;
            apnd = maxIndxY;
        } else if (maxIndxY < minIndxY) {
            b_n = 0;
            a = minIndxY;
            apnd = maxIndxY;
        } else if (muDoubleScalarIsInf(minIndxY) || muDoubleScalarIsInf(maxIndxY)) {
            b_n = 1;
            a = rtNaN;
            apnd = maxIndxY;
        } else {
            a = minIndxY;
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
                b_n = (int32_T)ndbl;
            } else {
                b_n = 0;
            }
        }
        i8 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = b_n;
        emxEnsureCapacity((emxArray__common *)y, i8, (int32_T)sizeof(real_T));
        if (b_n > 0) {
            y->data[0] = a;
            if (b_n > 1) {
                y->data[b_n - 1] = apnd;
                nm1 = b_n - 1;
                i8 = nm1;
                nx = (int32_T)((uint32_T)i8 >> 1);
                loop_ub = nx - 1;
                for (k = 1; k <= loop_ub; k++) {
                    y->data[k] = a + (real_T)k;
                    y->data[(b_n - k) - 1] = apnd - (real_T)k;
                }
                if (nx << 1 == nm1) {
                    y->data[nx] = (a + apnd) / 2.0;
                } else {
                    y->data[nx] = a + (real_T)nx;
                    y->data[nx + 1] = apnd - (real_T)nx;
                }
            }
        }
        i8 = c_y->size[0];
        c_y->size[0] = y->size[1];
        emxEnsureCapacity((emxArray__common *)c_y, i8, (int32_T)sizeof(real_T));
        loop_ub = y->size[1] - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            c_y->data[i8] = y->data[i8];
        }
        GaussListND_mexCode(c_y, user->psfSigma[0], d_x->data[i + d_x->size[0]], temp);
        /* 'objFcn:66' temp2 = squeeze(temp); */
        squeeze(temp, psfValueX);
        /* clear temp */
        /* 'objFcn:69' psfIntegY(:,i) = temp2(:,1); */
        loop_ub = psfValueX->size[0] - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            psfIntegY->data[i8 + psfIntegY->size[0] * i] = psfValueX->data[i8];
        }
        i++;
    }
    emxFree_real_T(&c_y);
    b_emxInit_real_T(&psfIntegZ, 2, TRUE);
    /* determine the contribution of each PSF (assuming amplitude 1) to a */
    /* pixel based on its z-coordinate (needed to calculate F & J) */
    /* 'objFcn:74' psfIntegZ = zeros(maxIndxZ-minIndxZ+1,n); */
    i8 = psfIntegZ->size[0] * psfIntegZ->size[1];
    psfIntegZ->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxZ - minIndxZ) + 1.0, &c_emlrtDCI, &emlrtContextGlobal);
    psfIntegZ->size[1] = n;
    emxEnsureCapacity((emxArray__common *)psfIntegZ, i8, (int32_T)sizeof(real_T));
    loop_ub = (int32_T)emlrtIntegerCheckR2011a((maxIndxZ - minIndxZ) + 1.0, &i_emlrtDCI, &emlrtContextGlobal) * n - 1;
    for (i8 = 0; i8 <= loop_ub; i8++) {
        psfIntegZ->data[i8] = 0.0;
    }
    /* 'objFcn:75' for i=1:n */
    i = 0;
    emxInit_real_T(&d_y, 1, TRUE);
    while (i + 1 <= n) {
        /* 'objFcn:76' temp = GaussListND_mexCode((minIndxZ:maxIndxZ)',... */
        /* 'objFcn:77'         psfSigma(2),psfPos(i,3)); */
        if (muDoubleScalarIsNaN(minIndxZ) || muDoubleScalarIsNaN(maxIndxZ)) {
            b_n = 1;
            a = rtNaN;
            apnd = maxIndxZ;
        } else if (maxIndxZ < minIndxZ) {
            b_n = 0;
            a = minIndxZ;
            apnd = maxIndxZ;
        } else if (muDoubleScalarIsInf(minIndxZ) || muDoubleScalarIsInf(maxIndxZ)) {
            b_n = 1;
            a = rtNaN;
            apnd = maxIndxZ;
        } else {
            a = minIndxZ;
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
                b_n = (int32_T)ndbl;
            } else {
                b_n = 0;
            }
        }
        i8 = y->size[0] * y->size[1];
        y->size[0] = 1;
        y->size[1] = b_n;
        emxEnsureCapacity((emxArray__common *)y, i8, (int32_T)sizeof(real_T));
        if (b_n > 0) {
            y->data[0] = a;
            if (b_n > 1) {
                y->data[b_n - 1] = apnd;
                nm1 = b_n - 1;
                i8 = nm1;
                nx = (int32_T)((uint32_T)i8 >> 1);
                loop_ub = nx - 1;
                for (k = 1; k <= loop_ub; k++) {
                    y->data[k] = a + (real_T)k;
                    y->data[(b_n - k) - 1] = apnd - (real_T)k;
                }
                if (nx << 1 == nm1) {
                    y->data[nx] = (a + apnd) / 2.0;
                } else {
                    y->data[nx] = a + (real_T)nx;
                    y->data[nx + 1] = apnd - (real_T)nx;
                }
            }
        }
        i8 = d_y->size[0];
        d_y->size[0] = y->size[1];
        emxEnsureCapacity((emxArray__common *)d_y, i8, (int32_T)sizeof(real_T));
        loop_ub = y->size[1] - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            d_y->data[i8] = y->data[i8];
        }
        GaussListND_mexCode(d_y, user->psfSigma[1], d_x->data[i + (d_x->size[0] << 1)], temp);
        /* 'objFcn:80' temp2 = squeeze(temp); */
        squeeze(temp, psfValueX);
        /* clear temp */
        /* 'objFcn:83' psfIntegZ(:,i) = temp2(:,1); */
        loop_ub = psfValueX->size[0] - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            psfIntegZ->data[i8 + psfIntegZ->size[0] * i] = psfValueX->data[i8];
        }
        i++;
    }
    emxFree_real_T(&d_y);
    emxFree_real_T(&temp);
    /* % PSF values for J */
    /* 'objFcn:87' if mode > 0 */
    b_emxInit_real_T(&psfValueY, 2, TRUE);
    b_emxInit_real_T(&psfValueZ, 2, TRUE);
    b_emxInit_real_T(&r0, 2, TRUE);
    emxInit_real_T(&r1, 1, TRUE);
    emxInit_real_T(&A, 1, TRUE);
    if (*mode > 0) {
        /* calculate the value of each PSF (assuming amplitude 1) at the */
        /* x-coordinates of the corners of all pixels (needed to calculate J) */
        /* 'objFcn:90' psfValueX = zeros(maxIndxX-minIndxX+2,n); */
        i8 = psfValueX->size[0] * psfValueX->size[1];
        psfValueX->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxX - minIndxX) + 2.0, &d_emlrtDCI, &emlrtContextGlobal);
        psfValueX->size[1] = n;
        emxEnsureCapacity((emxArray__common *)psfValueX, i8, (int32_T)sizeof(real_T));
        loop_ub = (int32_T)emlrtIntegerCheckR2011a((maxIndxX - minIndxX) + 2.0, &j_emlrtDCI, &emlrtContextGlobal) * n - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            psfValueX->data[i8] = 0.0;
        }
        /* 'objFcn:91' for i=1:n */
        i = 0;
        emxInit_real_T(&r2, 1, TRUE);
        emxInit_real_T(&e_y, 1, TRUE);
        while (i + 1 <= n) {
            /* 'objFcn:92' psfValueX(:,i) = exp(-((minIndxX-0.5:maxIndxX+0.5)'... */
            /* 'objFcn:93'             -psfPos(i,1)).^2/2/psfSigma(1)^2); */
            a = minIndxX - 0.5;
            d = maxIndxX + 0.5;
            if (muDoubleScalarIsNaN(a) || muDoubleScalarIsNaN(d)) {
                b_n = 1;
                a = rtNaN;
                apnd = d;
            } else if (d < a) {
                b_n = 0;
                apnd = d;
            } else if (muDoubleScalarIsInf(a) || muDoubleScalarIsInf(d)) {
                b_n = 1;
                a = rtNaN;
                apnd = d;
            } else {
                ndbl = muDoubleScalarFloor((d - a) + 0.5);
                apnd = a + ndbl;
                cdiff = apnd - d;
                absa = muDoubleScalarAbs(a);
                absb = muDoubleScalarAbs(d);
                if (absa > absb) {
                    absb = absa;
                }
                if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * absb) {
                    ndbl++;
                    apnd = d;
                } else if (cdiff > 0.0) {
                    apnd = a + (ndbl - 1.0);
                } else {
                    ndbl++;
                }
                if (ndbl >= 0.0) {
                    b_n = (int32_T)ndbl;
                } else {
                    b_n = 0;
                }
            }
            i8 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = b_n;
            emxEnsureCapacity((emxArray__common *)y, i8, (int32_T)sizeof(real_T));
            if (b_n > 0) {
                y->data[0] = a;
                if (b_n > 1) {
                    y->data[b_n - 1] = apnd;
                    nm1 = b_n - 1;
                    i8 = nm1;
                    nx = (int32_T)((uint32_T)i8 >> 1);
                    loop_ub = nx - 1;
                    for (k = 1; k <= loop_ub; k++) {
                        y->data[k] = a + (real_T)k;
                        y->data[(b_n - k) - 1] = apnd - (real_T)k;
                    }
                    if (nx << 1 == nm1) {
                        y->data[nx] = (a + apnd) / 2.0;
                    } else {
                        y->data[nx] = a + (real_T)nx;
                        y->data[nx + 1] = apnd - (real_T)nx;
                    }
                }
            }
            a = d_x->data[i];
            i8 = e_y->size[0];
            e_y->size[0] = y->size[1];
            emxEnsureCapacity((emxArray__common *)e_y, i8, (int32_T)sizeof(real_T));
            loop_ub = y->size[1] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                e_y->data[i8] = y->data[i8] - a;
            }
            power(e_y, A);
            i8 = A->size[0];
            emxEnsureCapacity((emxArray__common *)A, i8, (int32_T)sizeof(real_T));
            loop_ub = A->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                A->data[i8] = -A->data[i8];
            }
            rdivide(A, 2.0, r1);
            i8 = r2->size[0];
            r2->size[0] = r1->size[0];
            emxEnsureCapacity((emxArray__common *)r2, i8, (int32_T)sizeof(real_T));
            loop_ub = r1->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                r2->data[i8] = r1->data[i8];
            }
            rdivide(r2, muDoubleScalarPower(user->psfSigma[0], 2.0), r1);
            b_exp(r1);
            loop_ub = r1->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                psfValueX->data[i8 + psfValueX->size[0] * i] = r1->data[i8];
            }
            i++;
        }
        emxFree_real_T(&e_y);
        emxFree_real_T(&r2);
        /* calculate the value of each PSF (assuming amplitude 1) at the */
        /* y-coordinates of the corners of all pixels (needed to calculate J) */
        /* 'objFcn:98' psfValueY = zeros(maxIndxY-minIndxY+2,n); */
        i8 = psfValueY->size[0] * psfValueY->size[1];
        psfValueY->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxY - minIndxY) + 2.0, &e_emlrtDCI, &emlrtContextGlobal);
        psfValueY->size[1] = n;
        emxEnsureCapacity((emxArray__common *)psfValueY, i8, (int32_T)sizeof(real_T));
        loop_ub = (int32_T)emlrtIntegerCheckR2011a((maxIndxY - minIndxY) + 2.0, &k_emlrtDCI, &emlrtContextGlobal) * n - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            psfValueY->data[i8] = 0.0;
        }
        /* 'objFcn:99' for i=1:n */
        i = 0;
        emxInit_real_T(&r3, 1, TRUE);
        emxInit_real_T(&f_y, 1, TRUE);
        while (i + 1 <= n) {
            /* 'objFcn:100' psfValueY(:,i) = exp(-((minIndxY-0.5:maxIndxY+0.5)'... */
            /* 'objFcn:101'             -psfPos(i,2)).^2/2/psfSigma(1)^2); */
            a = minIndxY - 0.5;
            d = maxIndxY + 0.5;
            if (muDoubleScalarIsNaN(a) || muDoubleScalarIsNaN(d)) {
                b_n = 1;
                a = rtNaN;
                apnd = d;
            } else if (d < a) {
                b_n = 0;
                apnd = d;
            } else if (muDoubleScalarIsInf(a) || muDoubleScalarIsInf(d)) {
                b_n = 1;
                a = rtNaN;
                apnd = d;
            } else {
                ndbl = muDoubleScalarFloor((d - a) + 0.5);
                apnd = a + ndbl;
                cdiff = apnd - d;
                absa = muDoubleScalarAbs(a);
                absb = muDoubleScalarAbs(d);
                if (absa > absb) {
                    absb = absa;
                }
                if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * absb) {
                    ndbl++;
                    apnd = d;
                } else if (cdiff > 0.0) {
                    apnd = a + (ndbl - 1.0);
                } else {
                    ndbl++;
                }
                if (ndbl >= 0.0) {
                    b_n = (int32_T)ndbl;
                } else {
                    b_n = 0;
                }
            }
            i8 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = b_n;
            emxEnsureCapacity((emxArray__common *)y, i8, (int32_T)sizeof(real_T));
            if (b_n > 0) {
                y->data[0] = a;
                if (b_n > 1) {
                    y->data[b_n - 1] = apnd;
                    nm1 = b_n - 1;
                    i8 = nm1;
                    nx = (int32_T)((uint32_T)i8 >> 1);
                    loop_ub = nx - 1;
                    for (k = 1; k <= loop_ub; k++) {
                        y->data[k] = a + (real_T)k;
                        y->data[(b_n - k) - 1] = apnd - (real_T)k;
                    }
                    if (nx << 1 == nm1) {
                        y->data[nx] = (a + apnd) / 2.0;
                    } else {
                        y->data[nx] = a + (real_T)nx;
                        y->data[nx + 1] = apnd - (real_T)nx;
                    }
                }
            }
            a = d_x->data[i + d_x->size[0]];
            i8 = f_y->size[0];
            f_y->size[0] = y->size[1];
            emxEnsureCapacity((emxArray__common *)f_y, i8, (int32_T)sizeof(real_T));
            loop_ub = y->size[1] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                f_y->data[i8] = y->data[i8] - a;
            }
            power(f_y, A);
            i8 = A->size[0];
            emxEnsureCapacity((emxArray__common *)A, i8, (int32_T)sizeof(real_T));
            loop_ub = A->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                A->data[i8] = -A->data[i8];
            }
            rdivide(A, 2.0, r1);
            i8 = r3->size[0];
            r3->size[0] = r1->size[0];
            emxEnsureCapacity((emxArray__common *)r3, i8, (int32_T)sizeof(real_T));
            loop_ub = r1->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                r3->data[i8] = r1->data[i8];
            }
            rdivide(r3, muDoubleScalarPower(user->psfSigma[0], 2.0), r1);
            b_exp(r1);
            loop_ub = r1->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                psfValueY->data[i8 + psfValueY->size[0] * i] = r1->data[i8];
            }
            i++;
        }
        emxFree_real_T(&f_y);
        emxFree_real_T(&r3);
        /* calculate the value of each PSF (assuming amplitude 1) at the */
        /* z-coordinates of the corners of all pixels (needed to calculate J) */
        /* 'objFcn:106' psfValueZ = zeros(maxIndxZ-minIndxZ+2,n); */
        i8 = psfValueZ->size[0] * psfValueZ->size[1];
        psfValueZ->size[0] = (int32_T)emlrtIntegerCheckR2011a((maxIndxZ - minIndxZ) + 2.0, &f_emlrtDCI, &emlrtContextGlobal);
        psfValueZ->size[1] = n;
        emxEnsureCapacity((emxArray__common *)psfValueZ, i8, (int32_T)sizeof(real_T));
        loop_ub = (int32_T)emlrtIntegerCheckR2011a((maxIndxZ - minIndxZ) + 2.0, &l_emlrtDCI, &emlrtContextGlobal) * n - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            psfValueZ->data[i8] = 0.0;
        }
        /* 'objFcn:107' for i=1:n */
        i = 0;
        emxInit_real_T(&r4, 1, TRUE);
        emxInit_real_T(&g_y, 1, TRUE);
        while (i + 1 <= n) {
            /* 'objFcn:108' psfValueZ(:,i) = exp(-((minIndxZ-0.5:maxIndxZ+0.5)'... */
            /* 'objFcn:109'             -psfPos(i,3)).^2/2/psfSigma(2)^2); */
            a = minIndxZ - 0.5;
            d = maxIndxZ + 0.5;
            if (muDoubleScalarIsNaN(a) || muDoubleScalarIsNaN(d)) {
                b_n = 1;
                a = rtNaN;
                apnd = d;
            } else if (d < a) {
                b_n = 0;
                apnd = d;
            } else if (muDoubleScalarIsInf(a) || muDoubleScalarIsInf(d)) {
                b_n = 1;
                a = rtNaN;
                apnd = d;
            } else {
                ndbl = muDoubleScalarFloor((d - a) + 0.5);
                apnd = a + ndbl;
                cdiff = apnd - d;
                absa = muDoubleScalarAbs(a);
                absb = muDoubleScalarAbs(d);
                if (absa > absb) {
                    absb = absa;
                }
                if (muDoubleScalarAbs(cdiff) < 4.4408920985006262E-16 * absb) {
                    ndbl++;
                    apnd = d;
                } else if (cdiff > 0.0) {
                    apnd = a + (ndbl - 1.0);
                } else {
                    ndbl++;
                }
                if (ndbl >= 0.0) {
                    b_n = (int32_T)ndbl;
                } else {
                    b_n = 0;
                }
            }
            i8 = y->size[0] * y->size[1];
            y->size[0] = 1;
            y->size[1] = b_n;
            emxEnsureCapacity((emxArray__common *)y, i8, (int32_T)sizeof(real_T));
            if (b_n > 0) {
                y->data[0] = a;
                if (b_n > 1) {
                    y->data[b_n - 1] = apnd;
                    nm1 = b_n - 1;
                    i8 = nm1;
                    nx = (int32_T)((uint32_T)i8 >> 1);
                    loop_ub = nx - 1;
                    for (k = 1; k <= loop_ub; k++) {
                        y->data[k] = a + (real_T)k;
                        y->data[(b_n - k) - 1] = apnd - (real_T)k;
                    }
                    if (nx << 1 == nm1) {
                        y->data[nx] = (a + apnd) / 2.0;
                    } else {
                        y->data[nx] = a + (real_T)nx;
                        y->data[nx + 1] = apnd - (real_T)nx;
                    }
                }
            }
            a = d_x->data[i + (d_x->size[0] << 1)];
            i8 = g_y->size[0];
            g_y->size[0] = y->size[1];
            emxEnsureCapacity((emxArray__common *)g_y, i8, (int32_T)sizeof(real_T));
            loop_ub = y->size[1] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                g_y->data[i8] = y->data[i8] - a;
            }
            power(g_y, A);
            i8 = A->size[0];
            emxEnsureCapacity((emxArray__common *)A, i8, (int32_T)sizeof(real_T));
            loop_ub = A->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                A->data[i8] = -A->data[i8];
            }
            rdivide(A, 2.0, r1);
            i8 = r4->size[0];
            r4->size[0] = r1->size[0];
            emxEnsureCapacity((emxArray__common *)r4, i8, (int32_T)sizeof(real_T));
            loop_ub = r1->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                r4->data[i8] = r1->data[i8];
            }
            rdivide(r4, muDoubleScalarPower(user->psfSigma[1], 2.0), r1);
            b_exp(r1);
            loop_ub = r1->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                psfValueZ->data[i8 + psfValueZ->size[0] * i] = r1->data[i8];
            }
            i++;
        }
        emxFree_real_T(&g_y);
        emxFree_real_T(&r4);
        /*     %% j calculation */
        /* 'objFcn:113' fjac = ones(m,4*n+1); */
        nx = _s32_add__(_s32_shl_s32_(n, 2), 1);
        i8 = fjac->size[0] * fjac->size[1];
        fjac->size[0] = m;
        fjac->size[1] = nx;
        emxEnsureCapacity((emxArray__common *)fjac, i8, (int32_T)sizeof(real_T));
        loop_ub = nx - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            i = m - 1;
            for (nx = 0; nx <= i; nx++) {
                fjac->data[nx + fjac->size[0] * i8] = 1.0;
            }
        }
        /* (last column for background amplitude) */
        /* 'objFcn:114' fjac(:,1:4:4*n) = repmat(psfAmp',m,1).*(psfValueX(relIndxX,:)-... */
        /* 'objFcn:115'         psfValueX(relIndxX+1,:)).*psfIntegY(relIndxY,:).*psfIntegZ(relIndxZ,:); */
        loop_ub = relIndxX->size[0] - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            emlrtIntegerCheckR2011a(relIndxX->data[i8], &m_emlrtDCI, &emlrtContextGlobal);
        }
        b_emxInit_real_T(&e_x, 2, TRUE);
        i8 = e_x->size[0] * e_x->size[1];
        e_x->size[0] = 1;
        e_x->size[1] = d_x->size[0];
        emxEnsureCapacity((emxArray__common *)e_x, i8, (int32_T)sizeof(real_T));
        loop_ub = d_x->size[0] - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            e_x->data[e_x->size[0] * i8] = d_x->data[i8 + d_x->size[0] * 3];
        }
        b_repmat(e_x, m, r0);
        emxFree_real_T(&e_x);
        loop_ub = relIndxY->size[0] - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            emlrtIntegerCheckR2011a(relIndxY->data[i8], &n_emlrtDCI, &emlrtContextGlobal);
        }
        loop_ub = relIndxZ->size[0] - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            emlrtIntegerCheckR2011a(relIndxZ->data[i8], &o_emlrtDCI, &emlrtContextGlobal);
        }
        if (1 > _s32_shl_s32_(n, 2)) {
            i8 = 1;
        } else {
            i8 = 4;
        }
        loop_ub = r0->size[1] - 1;
        for (nx = 0; nx <= loop_ub; nx++) {
            i = r0->size[0] - 1;
            for (k = 0; k <= i; k++) {
                fjac->data[k + fjac->size[0] * (i8 * nx)] = r0->data[k + r0->size[0] * nx] * (psfValueX->data[((int32_T)relIndxX->data[k] + psfValueX->size[0] * nx) - 1] - psfValueX->data[((int32_T)(relIndxX->data[k] + 1.0) + psfValueX->size[0] * nx) - 1]) * psfIntegY->data[((int32_T)relIndxY->data[k] + psfIntegY->size[0] * nx) - 1] * psfIntegZ->data[((int32_T)relIndxZ->data[k] + psfIntegZ->size[0] * nx) - 1];
            }
        }
        b_emxInit_real_T(&f_x, 2, TRUE);
        /* w.r.t. x */
        /* 'objFcn:116' fjac(:,2:4:4*n) = repmat(psfAmp',m,1).*(psfValueY(relIndxY,:)-... */
        /* 'objFcn:117'         psfValueY(relIndxY+1,:)).*psfIntegX(relIndxX,:).*psfIntegZ(relIndxZ,:); */
        i8 = f_x->size[0] * f_x->size[1];
        f_x->size[0] = 1;
        f_x->size[1] = d_x->size[0];
        emxEnsureCapacity((emxArray__common *)f_x, i8, (int32_T)sizeof(real_T));
        loop_ub = d_x->size[0] - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            f_x->data[f_x->size[0] * i8] = d_x->data[i8 + d_x->size[0] * 3];
        }
        b_repmat(f_x, m, r0);
        emxFree_real_T(&f_x);
        if (2 > _s32_shl_s32_(n, 2)) {
            i8 = 0;
            nx = 1;
        } else {
            i8 = 1;
            nx = 4;
        }
        loop_ub = r0->size[1] - 1;
        for (k = 0; k <= loop_ub; k++) {
            i = r0->size[0] - 1;
            for (b_n = 0; b_n <= i; b_n++) {
                fjac->data[b_n + fjac->size[0] * (i8 + nx * k)] = r0->data[b_n + r0->size[0] * k] * (psfValueY->data[((int32_T)relIndxY->data[b_n] + psfValueY->size[0] * k) - 1] - psfValueY->data[((int32_T)(relIndxY->data[b_n] + 1.0) + psfValueY->size[0] * k) - 1]) * psfIntegX->data[((int32_T)relIndxX->data[b_n] + psfIntegX->size[0] * k) - 1] * psfIntegZ->data[((int32_T)relIndxZ->data[b_n] + psfIntegZ->size[0] * k) - 1];
            }
        }
        b_emxInit_real_T(&g_x, 2, TRUE);
        /* w.r.t. y */
        /* 'objFcn:118' fjac(:,3:4:4*n) = repmat(psfAmp',m,1).*(psfValueZ(relIndxZ,:)-... */
        /* 'objFcn:119'         psfValueZ(relIndxZ+1,:)).*psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:); */
        i8 = g_x->size[0] * g_x->size[1];
        g_x->size[0] = 1;
        g_x->size[1] = d_x->size[0];
        emxEnsureCapacity((emxArray__common *)g_x, i8, (int32_T)sizeof(real_T));
        loop_ub = d_x->size[0] - 1;
        for (i8 = 0; i8 <= loop_ub; i8++) {
            g_x->data[g_x->size[0] * i8] = d_x->data[i8 + d_x->size[0] * 3];
        }
        b_repmat(g_x, m, r0);
        emxFree_real_T(&g_x);
        if (3 > _s32_shl_s32_(n, 2)) {
            i8 = 0;
            nx = 1;
        } else {
            i8 = 2;
            nx = 4;
        }
        loop_ub = r0->size[1] - 1;
        for (k = 0; k <= loop_ub; k++) {
            i = r0->size[0] - 1;
            for (b_n = 0; b_n <= i; b_n++) {
                fjac->data[b_n + fjac->size[0] * (i8 + nx * k)] = r0->data[b_n + r0->size[0] * k] * (psfValueZ->data[((int32_T)relIndxZ->data[b_n] + psfValueZ->size[0] * k) - 1] - psfValueZ->data[((int32_T)(relIndxZ->data[b_n] + 1.0) + psfValueZ->size[0] * k) - 1]) * psfIntegX->data[((int32_T)relIndxX->data[b_n] + psfIntegX->size[0] * k) - 1] * psfIntegY->data[((int32_T)relIndxY->data[b_n] + psfIntegY->size[0] * k) - 1];
            }
        }
        /* w.r.t. z */
        /* 'objFcn:120' fjac(:,4:4:4*n) = psfIntegX(relIndxX,:).*psfIntegY(relIndxY,:).*psfIntegZ(relIndxZ,:); */
        if (4 > _s32_shl_s32_(n, 2)) {
            i8 = 0;
            nx = 1;
        } else {
            i8 = 3;
            nx = 4;
        }
        loop_ub = psfIntegX->size[1] - 1;
        for (k = 0; k <= loop_ub; k++) {
            i = relIndxX->size[0] - 1;
            for (b_n = 0; b_n <= i; b_n++) {
                fjac->data[b_n + fjac->size[0] * (i8 + nx * k)] = psfIntegX->data[((int32_T)relIndxX->data[b_n] + psfIntegX->size[0] * k) - 1] * psfIntegY->data[((int32_T)relIndxY->data[b_n] + psfIntegY->size[0] * k) - 1] * psfIntegZ->data[((int32_T)relIndxZ->data[b_n] + psfIntegZ->size[0] * k) - 1];
            }
        }
        /* w.r.t. amp */
    }
    /* % f calulation */
    /* 'objFcn:125' if mode ~= 1 */
    if (*mode != 1) {
        /* 'objFcn:126' if needfi > 0 */
        if (needfi > 0) {
            /* 'objFcn:127' f = zeros(m,1); */
            i8 = f->size[0] * f->size[1];
            f->size[0] = m;
            f->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)f, i8, (int32_T)sizeof(real_T));
            loop_ub = m - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                f->data[i8] = 0.0;
            }
            emxInit_real_T(&h_x, 1, TRUE);
            /* 'objFcn:128' f(needfi) = (sum(repmat(psfAmp,1,m).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)'.*psfIntegZ(relIndxZ,:)',1))' + repmat(bgAmp,m,1); */
            i8 = h_x->size[0];
            h_x->size[0] = d_x->size[0];
            emxEnsureCapacity((emxArray__common *)h_x, i8, (int32_T)sizeof(real_T));
            loop_ub = d_x->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                h_x->data[i8] = d_x->data[i8 + d_x->size[0] * 3];
            }
            c_repmat(h_x, m, r0);
            i8 = psfValueZ->size[0] * psfValueZ->size[1];
            psfValueZ->size[0] = psfIntegX->size[1];
            psfValueZ->size[1] = relIndxX->size[0];
            emxEnsureCapacity((emxArray__common *)psfValueZ, i8, (int32_T)sizeof(real_T));
            emxFree_real_T(&h_x);
            loop_ub = relIndxX->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                i = psfIntegX->size[1] - 1;
                for (nx = 0; nx <= i; nx++) {
                    psfValueZ->data[nx + psfValueZ->size[0] * i8] = psfIntegX->data[((int32_T)emlrtIntegerCheckR2011a(relIndxX->data[i8], &s_emlrtDCI, &emlrtContextGlobal) + psfIntegX->size[0] * nx) - 1];
                }
            }
            i8 = psfValueY->size[0] * psfValueY->size[1];
            psfValueY->size[0] = psfIntegY->size[1];
            psfValueY->size[1] = relIndxY->size[0];
            emxEnsureCapacity((emxArray__common *)psfValueY, i8, (int32_T)sizeof(real_T));
            loop_ub = relIndxY->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                i = psfIntegY->size[1] - 1;
                for (nx = 0; nx <= i; nx++) {
                    psfValueY->data[nx + psfValueY->size[0] * i8] = psfIntegY->data[((int32_T)emlrtIntegerCheckR2011a(relIndxY->data[i8], &t_emlrtDCI, &emlrtContextGlobal) + psfIntegY->size[0] * nx) - 1];
                }
            }
            i8 = psfValueX->size[0] * psfValueX->size[1];
            psfValueX->size[0] = psfIntegZ->size[1];
            psfValueX->size[1] = relIndxZ->size[0];
            emxEnsureCapacity((emxArray__common *)psfValueX, i8, (int32_T)sizeof(real_T));
            loop_ub = relIndxZ->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                i = psfIntegZ->size[1] - 1;
                for (nx = 0; nx <= i; nx++) {
                    psfValueX->data[nx + psfValueX->size[0] * i8] = psfIntegZ->data[((int32_T)emlrtIntegerCheckR2011a(relIndxZ->data[i8], &u_emlrtDCI, &emlrtContextGlobal) + psfIntegZ->size[0] * nx) - 1];
                }
            }
            b_emxInit_real_T(&r5, 2, TRUE);
            i8 = r5->size[0] * r5->size[1];
            r5->size[0] = r0->size[0];
            r5->size[1] = r0->size[1];
            emxEnsureCapacity((emxArray__common *)r5, i8, (int32_T)sizeof(real_T));
            loop_ub = r0->size[0] * r0->size[1] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                r5->data[i8] = r0->data[i8] * psfValueZ->data[i8] * psfValueY->data[i8] * psfValueX->data[i8];
            }
            sum(r5, y);
            i8 = r1->size[0];
            r1->size[0] = y->size[1];
            emxEnsureCapacity((emxArray__common *)r1, i8, (int32_T)sizeof(real_T));
            emxFree_real_T(&r5);
            loop_ub = y->size[1] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                r1->data[i8] = y->data[i8];
            }
            d_repmat(bgAmp, m, A);
            f->data[needfi - 1] = r1->data[-1] + A->data[-1];
        } else {
            emxInit_real_T(&i_x, 1, TRUE);
            /* 'objFcn:129' else */
            /* 'objFcn:130' f = (sum(repmat(psfAmp,1,m).*psfIntegX(relIndxX,:)'.*psfIntegY(relIndxY,:)'.*psfIntegZ(relIndxZ,:)',1))' + repmat(bgAmp,m,1); */
            i8 = i_x->size[0];
            i_x->size[0] = d_x->size[0];
            emxEnsureCapacity((emxArray__common *)i_x, i8, (int32_T)sizeof(real_T));
            loop_ub = d_x->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                i_x->data[i8] = d_x->data[i8 + d_x->size[0] * 3];
            }
            c_repmat(i_x, m, r0);
            i8 = psfValueZ->size[0] * psfValueZ->size[1];
            psfValueZ->size[0] = psfIntegX->size[1];
            psfValueZ->size[1] = relIndxX->size[0];
            emxEnsureCapacity((emxArray__common *)psfValueZ, i8, (int32_T)sizeof(real_T));
            emxFree_real_T(&i_x);
            loop_ub = relIndxX->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                i = psfIntegX->size[1] - 1;
                for (nx = 0; nx <= i; nx++) {
                    psfValueZ->data[nx + psfValueZ->size[0] * i8] = psfIntegX->data[((int32_T)emlrtIntegerCheckR2011a(relIndxX->data[i8], &p_emlrtDCI, &emlrtContextGlobal) + psfIntegX->size[0] * nx) - 1];
                }
            }
            i8 = psfValueY->size[0] * psfValueY->size[1];
            psfValueY->size[0] = psfIntegY->size[1];
            psfValueY->size[1] = relIndxY->size[0];
            emxEnsureCapacity((emxArray__common *)psfValueY, i8, (int32_T)sizeof(real_T));
            loop_ub = relIndxY->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                i = psfIntegY->size[1] - 1;
                for (nx = 0; nx <= i; nx++) {
                    psfValueY->data[nx + psfValueY->size[0] * i8] = psfIntegY->data[((int32_T)emlrtIntegerCheckR2011a(relIndxY->data[i8], &q_emlrtDCI, &emlrtContextGlobal) + psfIntegY->size[0] * nx) - 1];
                }
            }
            i8 = psfValueX->size[0] * psfValueX->size[1];
            psfValueX->size[0] = psfIntegZ->size[1];
            psfValueX->size[1] = relIndxZ->size[0];
            emxEnsureCapacity((emxArray__common *)psfValueX, i8, (int32_T)sizeof(real_T));
            loop_ub = relIndxZ->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                i = psfIntegZ->size[1] - 1;
                for (nx = 0; nx <= i; nx++) {
                    psfValueX->data[nx + psfValueX->size[0] * i8] = psfIntegZ->data[((int32_T)emlrtIntegerCheckR2011a(relIndxZ->data[i8], &r_emlrtDCI, &emlrtContextGlobal) + psfIntegZ->size[0] * nx) - 1];
                }
            }
            b_emxInit_real_T(&r6, 2, TRUE);
            i8 = r6->size[0] * r6->size[1];
            r6->size[0] = r0->size[0];
            r6->size[1] = r0->size[1];
            emxEnsureCapacity((emxArray__common *)r6, i8, (int32_T)sizeof(real_T));
            loop_ub = r0->size[0] * r0->size[1] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                r6->data[i8] = r0->data[i8] * psfValueZ->data[i8] * psfValueY->data[i8] * psfValueX->data[i8];
            }
            sum(r6, y);
            i8 = r1->size[0];
            r1->size[0] = y->size[1];
            emxEnsureCapacity((emxArray__common *)r1, i8, (int32_T)sizeof(real_T));
            emxFree_real_T(&r6);
            loop_ub = y->size[1] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                r1->data[i8] = y->data[i8];
            }
            d_repmat(bgAmp, m, A);
            nx = r1->size[0];
            i8 = f->size[0] * f->size[1];
            f->size[0] = nx;
            emxEnsureCapacity((emxArray__common *)f, i8, (int32_T)sizeof(real_T));
            i8 = f->size[0] * f->size[1];
            f->size[1] = 1;
            emxEnsureCapacity((emxArray__common *)f, i8, (int32_T)sizeof(real_T));
            loop_ub = r1->size[0] - 1;
            for (i8 = 0; i8 <= loop_ub; i8++) {
                f->data[i8] = r1->data[i8] + A->data[i8];
            }
        }
    }
    emxFree_real_T(&A);
    emxFree_real_T(&y);
    emxFree_real_T(&r1);
    emxFree_real_T(&r0);
    emxFree_real_T(&d_x);
    emxFree_real_T(&psfValueZ);
    emxFree_real_T(&psfValueY);
    emxFree_real_T(&psfValueX);
    emxFree_real_T(&psfIntegZ);
    emxFree_real_T(&psfIntegY);
    emxFree_real_T(&psfIntegX);
    emxFree_real_T(&relIndxZ);
    emxFree_real_T(&relIndxY);
    emxFree_real_T(&relIndxX);
    emlrtHeapReferenceStackLeaveFcn();
}
/* End of code generation (objFcn.c) */
