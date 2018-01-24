/*
 * multiGaussListND_mexCode_api.c
 *
 * Code generation for function 'multiGaussListND_mexCode_api'
 *
 * C source code generated on: Sun Dec  2 23:59:22 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "multiGaussListND_mexCode.h"
#include "multiGaussListND_mexCode_api.h"
#include "multiGaussListND_mexCode_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static void c_emlrt_marshallIn(const mxArray *X, const char_T *identifier,
  real_T (**y_data)[1500], int32_T y_size[2]);
static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T (**y_data)[1500], int32_T y_size[2]);
static void e_emlrt_marshallIn(const mxArray *sigma, const char_T *identifier,
  real_T (**y_data)[3], int32_T y_size[2]);
static void emlrt_marshallIn(const mxArray *maskCoord, const char_T *identifier,
  emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(emxArray_real_T *u);
static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T (**y_data)[3], int32_T y_size[2]);
static void g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);
static void h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T (**ret_data)[1500], int32_T ret_size[2]);
static void i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T (**ret_data)[3], int32_T ret_size[2]);

/* Function Definitions */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  g_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void c_emlrt_marshallIn(const mxArray *X, const char_T *identifier,
  real_T (**y_data)[1500], int32_T y_size[2])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  d_emlrt_marshallIn(emlrtAlias(X), &thisId, y_data, y_size);
  emlrtDestroyArray(&X);
}

static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T (**y_data)[1500], int32_T y_size[2])
{
  h_emlrt_marshallIn(emlrtAlias(u), parentId, y_data, y_size);
  emlrtDestroyArray(&u);
}

static void e_emlrt_marshallIn(const mxArray *sigma, const char_T *identifier,
  real_T (**y_data)[3], int32_T y_size[2])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  f_emlrt_marshallIn(emlrtAlias(sigma), &thisId, y_data, y_size);
  emlrtDestroyArray(&sigma);
}

static void emlrt_marshallIn(const mxArray *maskCoord, const char_T *identifier,
  emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  b_emlrt_marshallIn(emlrtAlias(maskCoord), &thisId, y);
  emlrtDestroyArray(&maskCoord);
}

static const mxArray *emlrt_marshallOut(emxArray_real_T *u)
{
  const mxArray *y;
  static const int32_T iv0[2] = { 0, 0 };

  const mxArray *m0;
  y = NULL;
  m0 = mxCreateNumericArray(2, (int32_T *)&iv0, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m0, (void *)u->data);
  mxSetDimensions((mxArray *)m0, u->size, 2);
  emlrtAssign(&y, m0);
  return y;
}

static void f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, real_T (**y_data)[3], int32_T y_size[2])
{
  i_emlrt_marshallIn(emlrtAlias(u), parentId, y_data, y_size);
  emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv1[2];
  boolean_T bv0[2];
  int32_T i2;
  for (i2 = 0; i2 < 2; i2++) {
    iv1[i2] = 500000 + -499997 * i2;
    bv0[i2] = TRUE;
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv1, bv0, ret->size);
  ret->size[0] = ret->size[0];
  ret->size[1] = ret->size[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static void h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T (**ret_data)[1500], int32_T ret_size[2])
{
  int32_T iv2[2];
  boolean_T bv1[2];
  int32_T i3;
  for (i3 = 0; i3 < 2; i3++) {
    iv2[i3] = 500 + -497 * i3;
    bv1[i3] = TRUE;
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv2, bv1, ret_size);
  *ret_data = (real_T (*)[1500])mxGetData(src);
  emlrtDestroyArray(&src);
}

static void i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, real_T (**ret_data)[3], int32_T ret_size[2])
{
  int32_T iv3[2];
  boolean_T bv2[2];
  int32_T i4;
  static const boolean_T bv3[2] = { FALSE, TRUE };

  for (i4 = 0; i4 < 2; i4++) {
    iv3[i4] = 1 + (i4 << 1);
    bv2[i4] = bv3[i4];
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv3, bv2, ret_size);
  ret_size[0] = 1;
  *ret_data = (real_T (*)[3])mxGetData(src);
  emlrtDestroyArray(&src);
}

void multiGaussListND_mexCode_api(const mxArray * const prhs[3], const mxArray
  *plhs[1])
{
  emxArray_real_T *maskCoord;
  emxArray_real_T *gaussList;
  int32_T X_size[2];
  real_T (*X_data)[1500];
  int32_T sigma_size[2];
  real_T (*sigma_data)[3];
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  b_emxInit_real_T(&maskCoord, 2, TRUE);
  b_emxInit_real_T(&gaussList, 2, TRUE);

  /* Marshall function inputs */
  emlrt_marshallIn(emlrtAlias(prhs[0]), "maskCoord", maskCoord);
  c_emlrt_marshallIn(emlrtAlias(prhs[1]), "X", &X_data, X_size);
  e_emlrt_marshallIn(emlrtAlias(prhs[2]), "sigma", &sigma_data, sigma_size);

  /* Invoke the target function */
  multiGaussListND_mexCode(maskCoord, *X_data, X_size, *sigma_data, sigma_size,
    gaussList);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(gaussList);
  gaussList->canFreeData = FALSE;
  emxFree_real_T(&gaussList);
  maskCoord->canFreeData = FALSE;
  emxFree_real_T(&maskCoord);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (multiGaussListND_mexCode_api.c) */
