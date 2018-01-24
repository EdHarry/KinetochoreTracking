/*
 * mmfMex_api.c
 *
 * Code generation for function 'mmfMex_api'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "mmfMex_api.h"
#include "mmfMex_emxutil.h"
#include "eml_warning.h"
#include "mmfMex_mexutil.h"

/* Variable Definitions */
static emlrtRTEInfo lb_emlrtRTEI = { 1, 1, "mmfMex_api", "" };

/* Function Declarations */
static const mxArray *b_emlrt_marshallOut(emxArray_real_T *u);
static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *maskAmp,
  const char_T *identifier, emxArray_real_T *y);
static const mxArray *c_emlrt_marshallOut(emxArray_real_T *u);
static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *maskCoord,
  const char_T *identifier, emxArray_real_T *y);
static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y);
static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *sigma, const
  char_T *identifier, real_T (**y_data)[3], int32_T y_size[2]);
static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T (**y_data)[3], int32_T y_size[2]);
static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret);
static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T (**ret_data)[3], int32_T ret_size[2]);

/* Function Definitions */
static const mxArray *b_emlrt_marshallOut(emxArray_real_T *u)
{
  const mxArray *y;
  static const int32_T iv28[2] = { 0, 0 };

  const mxArray *m15;
  y = NULL;
  m15 = mxCreateNumericArray(2, (int32_T *)&iv28, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m15, (void *)u->data);
  mxSetDimensions((mxArray *)m15, u->size, 2);
  emlrtAssign(&y, m15);
  return y;
}

static void c_emlrt_marshallIn(const emlrtStack *sp, const mxArray *maskAmp,
  const char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  d_emlrt_marshallIn(sp, emlrtAlias(maskAmp), &thisId, y);
  emlrtDestroyArray(&maskAmp);
}

static const mxArray *c_emlrt_marshallOut(emxArray_real_T *u)
{
  const mxArray *y;
  static const int32_T iv29[1] = { 0 };

  const mxArray *m16;
  y = NULL;
  m16 = mxCreateNumericArray(1, (int32_T *)&iv29, mxDOUBLE_CLASS, mxREAL);
  mxSetData((mxArray *)m16, (void *)u->data);
  mxSetDimensions((mxArray *)m16, u->size, 1);
  emlrtAssign(&y, m16);
  return y;
}

static void d_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  j_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void e_emlrt_marshallIn(const emlrtStack *sp, const mxArray *maskCoord,
  const char_T *identifier, emxArray_real_T *y)
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  f_emlrt_marshallIn(sp, emlrtAlias(maskCoord), &thisId, y);
  emlrtDestroyArray(&maskCoord);
}

static void f_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, emxArray_real_T *y)
{
  k_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y);
  emlrtDestroyArray(&u);
}

static void g_emlrt_marshallIn(const emlrtStack *sp, const mxArray *sigma, const
  char_T *identifier, real_T (**y_data)[3], int32_T y_size[2])
{
  emlrtMsgIdentifier thisId;
  thisId.fIdentifier = identifier;
  thisId.fParent = NULL;
  h_emlrt_marshallIn(sp, emlrtAlias(sigma), &thisId, y_data, y_size);
  emlrtDestroyArray(&sigma);
}

static void h_emlrt_marshallIn(const emlrtStack *sp, const mxArray *u, const
  emlrtMsgIdentifier *parentId, real_T (**y_data)[3], int32_T y_size[2])
{
  l_emlrt_marshallIn(sp, emlrtAlias(u), parentId, y_data, y_size);
  emlrtDestroyArray(&u);
}

static void j_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  int32_T iv31[1];
  boolean_T bv0[1];
  int32_T iv32[1];
  iv31[0] = -1;
  bv0[0] = TRUE;
  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", FALSE, 1U, iv31, bv0, iv32);
  ret->size[0] = iv32[0];
  ret->allocatedSize = ret->size[0];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static void k_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, emxArray_real_T *ret)
{
  int32_T iv33[2];
  boolean_T bv1[2];
  int32_T i8;
  int32_T iv34[2];
  for (i8 = 0; i8 < 2; i8++) {
    iv33[i8] = (i8 << 2) - 1;
    bv1[i8] = TRUE;
  }

  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", FALSE, 2U, iv33, bv1, iv34);
  ret->size[0] = iv34[0];
  ret->size[1] = iv34[1];
  ret->allocatedSize = ret->size[0] * ret->size[1];
  ret->data = (real_T *)mxGetData(src);
  ret->canFreeData = FALSE;
  emlrtDestroyArray(&src);
}

static void l_emlrt_marshallIn(const emlrtStack *sp, const mxArray *src, const
  emlrtMsgIdentifier *msgId, real_T (**ret_data)[3], int32_T ret_size[2])
{
  int32_T iv35[2];
  boolean_T bv2[2];
  int32_T i9;
  static const boolean_T bv3[2] = { FALSE, TRUE };

  int32_T iv36[2];
  for (i9 = 0; i9 < 2; i9++) {
    iv35[i9] = 1 + (i9 << 1);
    bv2[i9] = bv3[i9];
  }

  emlrtCheckVsBuiltInR2012b(sp, msgId, src, "double", FALSE, 2U, iv35, bv2, iv36);
  ret_size[0] = iv36[0];
  ret_size[1] = iv36[1];
  *ret_data = (real_T (*)[3])mxGetData(src);
  emlrtDestroyArray(&src);
}

void mmfMex_api(emlrtStack *sp, const mxArray * const prhs[4], const mxArray
                *plhs[6])
{
  emxArray_real_T *maskAmp;
  emxArray_real_T *maskCoord;
  emxArray_real_T *X;
  emxArray_real_T *jacobian;
  emxArray_real_T *resi;
  emxArray_real_T *amp;
  int32_T sigma_size[2];
  real_T (*sigma_data)[3];
  real_T nDim;
  real_T nGauss;
  real_T bg;
  emlrtHeapReferenceStackEnterFcnR2012b(sp);
  c_emxInit_real_T(sp, &maskAmp, 1, &lb_emlrtRTEI, TRUE);
  emxInit_real_T(sp, &maskCoord, 2, &lb_emlrtRTEI, TRUE);
  emxInit_real_T(sp, &X, 2, &lb_emlrtRTEI, TRUE);
  emxInit_real_T(sp, &jacobian, 2, &lb_emlrtRTEI, TRUE);
  c_emxInit_real_T(sp, &resi, 1, &lb_emlrtRTEI, TRUE);
  c_emxInit_real_T(sp, &amp, 1, &lb_emlrtRTEI, TRUE);

  /* Marshall function inputs */
  c_emlrt_marshallIn(sp, emlrtAlias(prhs[0]), "maskAmp", maskAmp);
  e_emlrt_marshallIn(sp, emlrtAlias(prhs[1]), "maskCoord", maskCoord);
  e_emlrt_marshallIn(sp, emlrtAlias(prhs[2]), "X", X);
  g_emlrt_marshallIn(sp, emlrtAlias(prhs[3]), "sigma", &sigma_data, sigma_size);

  /* Invoke the target function */
  mmfMex(sp, maskAmp, maskCoord, X, *sigma_data, sigma_size, jacobian, resi, amp,
         &bg, &nGauss, &nDim);

  /* Marshall function outputs */
  plhs[0] = b_emlrt_marshallOut(jacobian);
  plhs[1] = c_emlrt_marshallOut(resi);
  plhs[2] = c_emlrt_marshallOut(amp);
  plhs[3] = emlrt_marshallOut(bg);
  plhs[4] = emlrt_marshallOut(nGauss);
  plhs[5] = emlrt_marshallOut(nDim);
  amp->canFreeData = FALSE;
  emxFree_real_T(&amp);
  resi->canFreeData = FALSE;
  emxFree_real_T(&resi);
  jacobian->canFreeData = FALSE;
  emxFree_real_T(&jacobian);
  X->canFreeData = FALSE;
  emxFree_real_T(&X);
  maskCoord->canFreeData = FALSE;
  emxFree_real_T(&maskCoord);
  maskAmp->canFreeData = FALSE;
  emxFree_real_T(&maskAmp);
  emlrtHeapReferenceStackLeaveFcnR2012b(sp);
}

/* End of code generation (mmfMex_api.c) */
