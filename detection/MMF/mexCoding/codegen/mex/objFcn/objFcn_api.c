/*
 * objFcn_api.c
 *
 * Code generation for function 'objFcn_api'
 *
 * C source code generated on: Fri May 25 21:48:51 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "objFcn.h"
#include "objFcn_api.h"
#include "objFcn_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static int32_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId);
static const mxArray *b_emlrt_marshallOut(emxArray_real_T *u);
static void c_emlrt_marshallIn(const mxArray *x, const char_T *identifier, emxArray_real_T *y);
static const mxArray *c_emlrt_marshallOut(const emxArray_real_T *u_index, const real_T u_psfSigma[2]);
static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void e_emlrt_marshallIn(const mxArray *fjac, const char_T *identifier, emxArray_real_T *y);
static int32_T emlrt_marshallIn(const mxArray *mode, const char_T *identifier);
static const mxArray *emlrt_marshallOut(int32_T u);
static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void g_emlrt_marshallIn(const mxArray *user, const char_T *identifier, emxArray_real_T *y_index, real_T y_psfSigma[2]);
static void h_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y_index, real_T y_psfSigma[2]);
static void i_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void j_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, real_T y[2]);
static int32_T k_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId);
static void l_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void m_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void n_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void o_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, real_T ret[2]);

/* Function Definitions */

static int32_T b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId)
{
    int32_T y;
    y = k_emlrt_marshallIn(emlrtAlias(u), parentId);
    emlrtDestroyArray(&u);
    return y;
}

static const mxArray *b_emlrt_marshallOut(emxArray_real_T *u)
{
    const mxArray *y;
    const mxArray *m1;
    real_T (*pData)[];
    int32_T i3;
    int32_T i;
    int32_T b_i;
    y = NULL;
    m1 = mxCreateNumericArray(2, u->size, mxDOUBLE_CLASS, mxREAL);
    pData = (real_T (*)[])mxGetPr(m1);
    i3 = 0;
    for (i = 0; i < u->size[1]; i++) {
        for (b_i = 0; b_i < u->size[0]; b_i++) {
            (*pData)[i3] = u->data[b_i + u->size[0] * i];
            i3++;
        }
    }
    emlrtAssign(&y, m1);
    return y;
}

static void c_emlrt_marshallIn(const mxArray *x, const char_T *identifier, emxArray_real_T *y)
{
    emlrtMsgIdentifier thisId;
    thisId.fIdentifier = identifier;
    thisId.fParent = NULL;
    d_emlrt_marshallIn(emlrtAlias(x), &thisId, y);
    emlrtDestroyArray(&x);
}

static const mxArray *c_emlrt_marshallOut(const emxArray_real_T *u_index, const real_T u_psfSigma[2])
{
    const mxArray *y;
    emxArray_real_T *u;
    int32_T i4;
    int32_T loop_ub;
    const mxArray *b_y;
    const mxArray *m2;
    real_T (*pData)[];
    int32_T i;
    const mxArray *c_y;
    static const int32_T iv2[2] = { 1, 2 };
    emlrtHeapReferenceStackEnterFcn();
    b_emxInit_real_T(&u, 2, TRUE);
    y = NULL;
    emlrtAssign(&y, mxCreateStructMatrix(1, 1, 0, NULL));
    i4 = u->size[0] * u->size[1];
    u->size[0] = u_index->size[0];
    u->size[1] = 3;
    emxEnsureCapacity((emxArray__common *)u, i4, (int32_T)sizeof(real_T));
    loop_ub = u_index->size[0] * u_index->size[1] - 1;
    for (i4 = 0; i4 <= loop_ub; i4++) {
        u->data[i4] = u_index->data[i4];
    }
    b_y = NULL;
    m2 = mxCreateNumericArray(2, u->size, mxDOUBLE_CLASS, mxREAL);
    pData = (real_T (*)[])mxGetPr(m2);
    i4 = 0;
    for (loop_ub = 0; loop_ub < 3; loop_ub++) {
        for (i = 0; i < u->size[0]; i++) {
            (*pData)[i4] = u->data[i + u->size[0] * loop_ub];
            i4++;
        }
    }
    emxFree_real_T(&u);
    emlrtAssign(&b_y, m2);
    emlrtAddField(y, b_y, "index", 0);
    c_y = NULL;
    m2 = mxCreateNumericArray(2, (int32_T *)&iv2, mxDOUBLE_CLASS, mxREAL);
    pData = (real_T (*)[])mxGetPr(m2);
    for (loop_ub = 0; loop_ub < 2; loop_ub++) {
        (*pData)[loop_ub] = u_psfSigma[loop_ub];
    }
    emlrtAssign(&c_y, m2);
    emlrtAddField(y, c_y, "psfSigma", 0);
    emlrtHeapReferenceStackLeaveFcn();
    return y;
}

static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
    l_emlrt_marshallIn(emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

static void e_emlrt_marshallIn(const mxArray *fjac, const char_T *identifier, emxArray_real_T *y)
{
    emlrtMsgIdentifier thisId;
    thisId.fIdentifier = identifier;
    thisId.fParent = NULL;
    f_emlrt_marshallIn(emlrtAlias(fjac), &thisId, y);
    emlrtDestroyArray(&fjac);
}

static int32_T emlrt_marshallIn(const mxArray *mode, const char_T *identifier)
{
    int32_T y;
    emlrtMsgIdentifier thisId;
    thisId.fIdentifier = identifier;
    thisId.fParent = NULL;
    y = b_emlrt_marshallIn(emlrtAlias(mode), &thisId);
    emlrtDestroyArray(&mode);
    return y;
}

static const mxArray *emlrt_marshallOut(int32_T u)
{
    const mxArray *y;
    const mxArray *m0;
    y = NULL;
    m0 = mxCreateNumericMatrix(1, 1, mxINT32_CLASS, mxREAL);
    *(int32_T *)mxGetData(m0) = u;
    emlrtAssign(&y, m0);
    return y;
}

static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
    m_emlrt_marshallIn(emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const mxArray *user, const char_T *identifier, emxArray_real_T *y_index, real_T y_psfSigma[2])
{
    emlrtMsgIdentifier thisId;
    thisId.fIdentifier = identifier;
    thisId.fParent = NULL;
    h_emlrt_marshallIn(emlrtAlias(user), &thisId, y_index, y_psfSigma);
    emlrtDestroyArray(&user);
}

static void h_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y_index, real_T y_psfSigma[2])
{
    emlrtMsgIdentifier thisId;
    static const char * fieldNames[2] = { "index", "psfSigma" };
    thisId.fParent = parentId;
    emlrtCheckStructR2011a(parentId, u, 2, fieldNames, 0U, 0);
    thisId.fIdentifier = "index";
    i_emlrt_marshallIn(emlrtAlias(emlrtGetField(u, 0, "index")), &thisId, y_index);
    thisId.fIdentifier = "psfSigma";
    j_emlrt_marshallIn(emlrtAlias(emlrtGetField(u, 0, "psfSigma")), &thisId, y_psfSigma);
    emlrtDestroyArray(&u);
}

static void i_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
    n_emlrt_marshallIn(emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

static void j_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier *parentId, real_T y[2])
{
    o_emlrt_marshallIn(emlrtAlias(u), parentId, y);
    emlrtDestroyArray(&u);
}

static int32_T k_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId)
{
    int32_T ret;
    emlrtCheckBuiltInR2011a(msgId, src, "int32", FALSE, 0U, 0);
    ret = *(int32_T *)mxGetData(src);
    emlrtDestroyArray(&src);
    return ret;
}

static void l_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
    int32_T iv3[1];
    boolean_T bv0[1];
    int32_T i5;
    iv3[0] = -1;
    bv0[0] = TRUE;
    emlrtCheckVsBuiltInR2011a(msgId, src, "double", FALSE, 1U, iv3, bv0, ret->size);
    i5 = ret->size[0];
    emxEnsureCapacity((emxArray__common *)ret, i5, (int32_T)sizeof(real_T));
    emlrtImportArrayR2008b(src, ret->data, 8);
    emlrtDestroyArray(&src);
}

static void m_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
    int32_T i;
    int32_T iv4[2];
    boolean_T bv1[2];
    for (i = 0; i < 2; i++) {
        iv4[i] = -1;
        bv1[i] = TRUE;
    }
    emlrtCheckVsBuiltInR2011a(msgId, src, "double", FALSE, 2U, iv4, bv1, ret->size);
    i = ret->size[0] * ret->size[1];
    emxEnsureCapacity((emxArray__common *)ret, i, (int32_T)sizeof(real_T));
    emlrtImportArrayR2008b(src, ret->data, 8);
    emlrtDestroyArray(&src);
}

static void n_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
    int32_T i6;
    int32_T iv5[2];
    static const boolean_T bv2[2] = { TRUE, FALSE };
    boolean_T bv3[2];
    for (i6 = 0; i6 < 2; i6++) {
        iv5[i6] = (i6 << 2) - 1;
        bv3[i6] = bv2[i6];
    }
    emlrtCheckVsBuiltInR2011a(msgId, src, "double", FALSE, 2U, iv5, bv3, ret->size);
    i6 = ret->size[0] * ret->size[1];
    ret->size[1] = 3;
    emxEnsureCapacity((emxArray__common *)ret, i6, (int32_T)sizeof(real_T));
    emlrtImportArrayR2008b(src, ret->data, 8);
    emlrtDestroyArray(&src);
}

static void o_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *msgId, real_T ret[2])
{
    int32_T i7;
    int32_T iv6[2];
    for (i7 = 0; i7 < 2; i7++) {
        iv6[i7] = 1 + i7;
    }
    emlrtCheckBuiltInR2011a(msgId, src, "double", FALSE, 2U, iv6);
    for (i7 = 0; i7 < 2; i7++) {
        ret[i7] = (*(real_T (*)[2])mxGetData(src))[i7];
    }
    emlrtDestroyArray(&src);
}

void objFcn_api(const mxArray * const prhs[9], const mxArray *plhs[4])
{
    emxArray_real_T *x;
    emxArray_real_T *fjac;
    struct_T user;
    emxArray_real_T *f;
    int32_T mode;
    int32_T m;
    int32_T n;
    int32_T ldfj;
    int32_T needfi;
    int32_T nstate;
    emlrtHeapReferenceStackEnterFcn();
    emxInit_real_T(&x, 1, TRUE);
    b_emxInit_real_T(&fjac, 2, TRUE);
    emxInitStruct_struct_T(&user, TRUE);
    b_emxInit_real_T(&f, 2, TRUE);
    /* Marshall function inputs */
    mode = emlrt_marshallIn(emlrtAliasP(prhs[0]), "mode");
    m = emlrt_marshallIn(emlrtAliasP(prhs[1]), "m");
    n = emlrt_marshallIn(emlrtAliasP(prhs[2]), "n");
    ldfj = emlrt_marshallIn(emlrtAliasP(prhs[3]), "ldfj");
    needfi = emlrt_marshallIn(emlrtAliasP(prhs[4]), "needfi");
    c_emlrt_marshallIn(emlrtAliasP(prhs[5]), "x", x);
    e_emlrt_marshallIn(emlrtAliasP(prhs[6]), "fjac", fjac);
    nstate = emlrt_marshallIn(emlrtAliasP(prhs[7]), "nstate");
    g_emlrt_marshallIn(emlrtAliasP(prhs[8]), "user", user.index, user.psfSigma);
    /* Invoke the target function */
    objFcn(&mode, m, n, ldfj, needfi, x, fjac, nstate, &user, f);
    /* Marshall function outputs */
    plhs[0] = emlrt_marshallOut(mode);
    plhs[1] = b_emlrt_marshallOut(f);
    plhs[2] = b_emlrt_marshallOut(fjac);
    plhs[3] = c_emlrt_marshallOut(user.index, user.psfSigma);
    emxFree_real_T(&f);
    emxFree_real_T(&fjac);
    emxFree_real_T(&x);
    emxFreeStruct_struct_T(&user);
    emlrtHeapReferenceStackLeaveFcn();
}
/* End of code generation (objFcn_api.c) */
