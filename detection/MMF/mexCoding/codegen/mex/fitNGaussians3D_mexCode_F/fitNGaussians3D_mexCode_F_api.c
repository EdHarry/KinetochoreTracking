/*
 * fitNGaussians3D_mexCode_F_api.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode_F_api'
 *
 * C source code generated on: Mon May 21 17:33:35 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_F.h"
#include "fitNGaussians3D_mexCode_F_api.h"
#include "fitNGaussians3D_mexCode_F_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void c_emlrt_marshallIn(const mxArray *b_index, const char_T *identifier, emxArray_real_T *y);
static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void e_emlrt_marshallIn(const mxArray *psfSigma, const char_T *identifier, real_T y[2]);
static void emlrt_marshallIn(const mxArray *x0, const char_T *identifier, emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(emxArray_real_T *u);
static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, real_T y[2]);
static void g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, real_T ret[2]);

/* Function Definitions */

static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
    g_emlrt_marshallIn(emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

static void c_emlrt_marshallIn(const mxArray *b_index, const char_T *identifier, emxArray_real_T *y)
{
    emlrtMsgIdentifier thisId;
    thisId.fIdentifier = identifier;
    thisId.fParent = NULL;
    d_emlrt_marshallIn(emlrtAlias(b_index), &thisId, y);
    emlrtDestroyArray(&b_index);
}

static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
    h_emlrt_marshallIn(emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

static void e_emlrt_marshallIn(const mxArray *psfSigma, const char_T *identifier, real_T y[2])
{
    emlrtMsgIdentifier thisId;
    thisId.fIdentifier = identifier;
    thisId.fParent = NULL;
    f_emlrt_marshallIn(emlrtAlias(psfSigma), &thisId, y);
    emlrtDestroyArray(&psfSigma);
}

static void emlrt_marshallIn(const mxArray *x0, const char_T *identifier, emxArray_real_T *y)
{
    emlrtMsgIdentifier thisId;
    thisId.fIdentifier = identifier;
    thisId.fParent = NULL;
    b_emlrt_marshallIn(emlrtAlias(x0), &thisId, y);
    emlrtDestroyArray(&x0);
}

static const mxArray *emlrt_marshallOut(emxArray_real_T *u)
{
    const mxArray *y;
    const mxArray *m0;
    real_T (*pData)[];
    int32_T i0;
    int32_T i;
    int32_T b_i;
    y = NULL;
    m0 = mxCreateNumericArray(2, u->size, mxDOUBLE_CLASS, mxREAL);
    pData = (real_T (*)[])mxGetPr(m0);
    i0 = 0;
    for (i = 0; i < u->size[1]; i++) {
        for (b_i = 0; b_i < u->size[0]; b_i++) {
            (*pData)[i0] = u->data[b_i + u->size[0] * i];
            i0++;
        }
    }
    emlrtAssign(&y, m0);
    return y;
}

static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, real_T y[2])
{
    i_emlrt_marshallIn(emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
    int32_T iv3[1];
    boolean_T bv0[1];
    int32_T i1;
    iv3[0] = -1;
    bv0[0] = TRUE;
    emlrtCheckVsBuiltInR2011a(msgId, src, "double", FALSE, 1U, iv3, bv0, ret->size);
    i1 = ret->size[0];
    emxEnsureCapacity((emxArray__common *)ret, i1, (int32_T)sizeof(real_T));
    emlrtImportArrayR2008b(src, ret->data, 8);
    emlrtDestroyArray(&src);
}

static void h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
    int32_T i2;
    int32_T iv4[2];
    static const boolean_T bv1[2] = { TRUE, FALSE };
    boolean_T bv2[2];
    for (i2 = 0; i2 < 2; i2++) {
        iv4[i2] = (i2 << 2) - 1;
        bv2[i2] = bv1[i2];
    }
    emlrtCheckVsBuiltInR2011a(msgId, src, "double", FALSE, 2U, iv4, bv2, ret->size);
    i2 = ret->size[0] * ret->size[1];
    ret->size[1] = 3;
    emxEnsureCapacity((emxArray__common *)ret, i2, (int32_T)sizeof(real_T));
    emlrtImportArrayR2008b(src, ret->data, 8);
    emlrtDestroyArray(&src);
}

static void i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, real_T ret[2])
{
    int32_T i3;
    int32_T iv5[2];
    for (i3 = 0; i3 < 2; i3++) {
        iv5[i3] = 1 + i3;
    }
    emlrtCheckBuiltInR2011a(msgId, src, "double", FALSE, 2U, iv5);
    for (i3 = 0; i3 < 2; i3++) {
        ret[i3] = (*(real_T (*)[2])mxGetData(src))[i3];
    }
    emlrtDestroyArray(&src);
}

void fitNGaussians3D_mexCode_F_api(const mxArray * const prhs[4], const mxArray *plhs[1])
{
    emxArray_real_T *x0;
    emxArray_real_T *image;
    emxArray_real_T *b_index;
    emxArray_real_T *F;
    real_T psfSigma[2];
    emlrtHeapReferenceStackEnterFcn();
    c_emxInit_real_T(&x0, 1, TRUE);
    c_emxInit_real_T(&image, 1, TRUE);
    emxInit_real_T(&b_index, 2, TRUE);
    emxInit_real_T(&F, 2, TRUE);
    /* Marshall function inputs */
    emlrt_marshallIn(emlrtAliasP(prhs[0]), "x0", x0);
    emlrt_marshallIn(emlrtAliasP(prhs[1]), "image", image);
    c_emlrt_marshallIn(emlrtAliasP(prhs[2]), "index", b_index);
    e_emlrt_marshallIn(emlrtAliasP(prhs[3]), "psfSigma", psfSigma);
    /* Invoke the target function */
    fitNGaussians3D_mexCode_F(x0, image, b_index, psfSigma, F);
    /* Marshall function outputs */
    plhs[0] = emlrt_marshallOut(F);
    emxFree_real_T(&F);
    emxFree_real_T(&b_index);
    emxFree_real_T(&image);
    emxFree_real_T(&x0);
    emlrtHeapReferenceStackLeaveFcn();
}
/* End of code generation (fitNGaussians3D_mexCode_F_api.c) */
