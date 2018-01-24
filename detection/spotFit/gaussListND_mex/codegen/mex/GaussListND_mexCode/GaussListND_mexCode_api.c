/*
 * GaussListND_mexCode_api.c
 *
 * Code generation for function 'GaussListND_mexCode_api'
 *
 * C source code generated on: Sun Dec  2 22:57:02 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "GaussListND_mexCode.h"
#include "GaussListND_mexCode_api.h"
#include "GaussListND_mexCode_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static void c_emlrt_marshallIn(const mxArray *sigma, const char_T *identifier,
  emxArray_real_T *y);
static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static void e_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);
static void emlrt_marshallIn(const mxArray *coordList, const char_T *identifier,
  emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(emxArray_real_T *u);
static void f_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);

/* Function Definitions */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  e_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void c_emlrt_marshallIn(const mxArray *sigma, const char_T *identifier,
  emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  d_emlrt_marshallIn(emlrtAlias(sigma), &thisId, y);
  emlrtDestroyArray(&sigma);
}

static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  f_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void e_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv1[2];
  boolean_T bv0[2];
  int32_T i;
  for (i = 0; i < 2; i++) {
    iv1[i] = -1;
    bv0[i] = TRUE;
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

static void emlrt_marshallIn(const mxArray *coordList, const char_T *identifier,
  emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  b_emlrt_marshallIn(emlrtAlias(coordList), &thisId, y);
  emlrtDestroyArray(&coordList);
}

static const mxArray *emlrt_marshallOut(emxArray_real_T *u)
{
  const mxArray *y;
  static const int32_T iv0[3] = { 0, 0, 0 };

  const mxArray *m0;
  y = NULL;
  m0 = mxCreateNumericArray(3, (int32_T *)&iv0, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m0, (void *)u->data);
  mxSetDimensions((mxArray *)m0, u->size, 3);
  emlrtAssign(&y, m0);
  return y;
}

static void f_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv2[2];
  boolean_T bv1[2];
  int32_T i3;
  static const boolean_T bv2[2] = { FALSE, TRUE };

  for (i3 = 0; i3 < 2; i3++) {
    iv2[i3] = 1 + -2 * i3;
    bv1[i3] = bv2[i3];
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv2, bv1, ret->size);
  ret->size[0] = 1;
  ret->size[1] = ret->size[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

void GaussListND_mexCode_api(const mxArray * const prhs[3], const mxArray *plhs
  [1])
{
  emxArray_real_T *coordList;
  emxArray_real_T *sigma;
  emxArray_real_T *center;
  emxArray_real_T *gaussList;
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  emxInit_real_T(&coordList, 2, TRUE);
  emxInit_real_T(&sigma, 2, TRUE);
  emxInit_real_T(&center, 2, TRUE);
  b_emxInit_real_T(&gaussList, 3, TRUE);

  /* Marshall function inputs */
  emlrt_marshallIn(emlrtAlias(prhs[0]), "coordList", coordList);
  c_emlrt_marshallIn(emlrtAlias(prhs[1]), "sigma", sigma);
  emlrt_marshallIn(emlrtAlias(prhs[2]), "center", center);

  /* Invoke the target function */
  GaussListND_mexCode(coordList, sigma, center, gaussList);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(gaussList);
  gaussList->canFreeData = FALSE;
  emxFree_real_T(&gaussList);
  center->canFreeData = FALSE;
  emxFree_real_T(&center);
  sigma->canFreeData = FALSE;
  emxFree_real_T(&sigma);
  coordList->canFreeData = FALSE;
  emxFree_real_T(&coordList);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (GaussListND_mexCode_api.c) */
