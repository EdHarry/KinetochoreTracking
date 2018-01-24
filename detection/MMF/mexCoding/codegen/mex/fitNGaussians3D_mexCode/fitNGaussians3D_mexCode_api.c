/*
 * fitNGaussians3D_mexCode_api.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode_api'
 *
 * C source code generated on: Tue Nov 19 11:16:20 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode.h"
#include "fitNGaussians3D_mexCode_api.h"
#include "fitNGaussians3D_mexCode_emxutil.h"

/* Function Declarations */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static void c_emlrt_marshallIn(const mxArray *b_index, const char_T *identifier,
  emxArray_real_T *y);
static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y);
static real_T (*e_emlrt_marshallIn(const mxArray *psfSigma, const char_T
  *identifier))[2];
static void emlrt_marshallIn(const mxArray *x0, const char_T *identifier,
  emxArray_real_T *y);
static const mxArray *emlrt_marshallOut(emxArray_real_T *u);
static real_T (*f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId))[2];
static void g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);
static void h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret);
static real_T (*i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *
  msgId))[2];

/* Function Definitions */
static void b_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  g_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void c_emlrt_marshallIn(const mxArray *b_index, const char_T *identifier,
  emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  d_emlrt_marshallIn(emlrtAlias(b_index), &thisId, y);
  emlrtDestroyArray(&b_index);
}

static void d_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId, emxArray_real_T *y)
{
  h_emlrt_marshallIn(emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static real_T (*e_emlrt_marshallIn(const mxArray *psfSigma, const char_T
  *identifier))[2]
{
  real_T (*y)[2];
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  y = f_emlrt_marshallIn(emlrtAlias(psfSigma), &thisId);
  emlrtDestroyArray(&psfSigma);
  return y;
}
  static void emlrt_marshallIn(const mxArray *x0, const char_T *identifier,
  emxArray_real_T *y)
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
  static const int32_T iv2[2] = { 0, 0 };

  const mxArray *m0;
  y = NULL;
  m0 = mxCreateNumericArray(2, (int32_T *)&iv2, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m0, (void *)u->data);
  mxSetDimensions((mxArray *)m0, u->size, 2);
  emlrtAssign(&y, m0);
  return y;
}

static real_T (*f_emlrt_marshallIn(const mxArray *u, const emlrtMsgIdentifier
  *parentId))[2]
{
  real_T (*y)[2];
  y = i_emlrt_marshallIn(emlrtAlias(u), parentId);
  emlrtDestroyArray(&u);
  return y;
}
  static void g_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv3[1];
  boolean_T bv0[1];
  int32_T iv4[1];
  iv3[0] = -1;
  bv0[0] = TRUE;
  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 1U,
    iv3, bv0, iv4);
  ret->size[0] = iv4[0];
  ret->allocatedSize = ret->size[0];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static void h_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier
  *msgId, emxArray_real_T *ret)
{
  int32_T iv5[2];
  boolean_T bv1[2];
  int32_T i4;
  static const boolean_T bv2[2] = { TRUE, FALSE };

  int32_T iv6[2];
  for (i4 = 0; i4 < 2; i4++) {
    iv5[i4] = (i4 << 2) - 1;
    bv1[i4] = bv2[i4];
  }

  emlrtCheckVsBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv5, bv1, iv6);
  ret->size[0] = iv6[0];
  ret->size[1] = iv6[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static real_T (*i_emlrt_marshallIn(const mxArray *src, const emlrtMsgIdentifier *
  msgId))[2]
{
  real_T (*ret)[2];
  int32_T iv7[2];
  int32_T i5;
  for (i5 = 0; i5 < 2; i5++) {
    iv7[i5] = 1 + i5;
  }

  emlrtCheckBuiltInR2012b(emlrtRootTLSGlobal, msgId, src, "double", FALSE, 2U,
    iv7);
  ret = (real_T (*)[2])mxGetData(src);
  emlrtDestroyArray(&src);
  return ret;
}
  void fitNGaussians3D_mexCode_api(const mxArray *prhs[4], const mxArray *plhs[2])
{
  emxArray_real_T *x0;
  emxArray_real_T *image;
  emxArray_real_T *b_index;
  emxArray_real_T *F;
  emxArray_real_T *J;
  real_T (*psfSigma)[2];
  emlrtHeapReferenceStackEnterFcnR2012b(emlrtRootTLSGlobal);
  emxInit_real_T(&x0, 1, TRUE);
  emxInit_real_T(&image, 1, TRUE);
  b_emxInit_real_T(&b_index, 2, TRUE);
  b_emxInit_real_T(&F, 2, TRUE);
  b_emxInit_real_T(&J, 2, TRUE);
  prhs[0] = emlrtProtectR2012b(prhs[0], 0, FALSE, -1);

  /* Marshall function inputs */
  emlrt_marshallIn(emlrtAlias(prhs[0]), "x0", x0);
  emlrt_marshallIn(emlrtAlias(prhs[1]), "image", image);
  c_emlrt_marshallIn(emlrtAlias(prhs[2]), "index", b_index);
  psfSigma = e_emlrt_marshallIn(emlrtAlias(prhs[3]), "psfSigma");

  /* Invoke the target function */
  fitNGaussians3D_mexCode(x0, image, b_index, *psfSigma, F, J);

  /* Marshall function outputs */
  plhs[0] = emlrt_marshallOut(F);
  plhs[1] = emlrt_marshallOut(J);
  J->canFreeData = FALSE;
  emxFree_real_T(&J);
  F->canFreeData = FALSE;
  emxFree_real_T(&F);
  b_index->canFreeData = FALSE;
  emxFree_real_T(&b_index);
  image->canFreeData = FALSE;
  emxFree_real_T(&image);
  x0->canFreeData = FALSE;
  emxFree_real_T(&x0);
  emlrtHeapReferenceStackLeaveFcnR2012b(emlrtRootTLSGlobal);
}

/* End of code generation (fitNGaussians3D_mexCode_api.c) */
