/*
 * multiGaussListND_mexCode_mex.c
 *
 * Code generation for function 'multiGaussListND_mexCode'
 *
 * C source code generated on: Sun Dec  2 23:59:22 2012
 *
 */

/* Include files */
#include "mex.h"
#include "multiGaussListND_mexCode_api.h"
#include "multiGaussListND_mexCode_initialize.h"
#include "multiGaussListND_mexCode_terminate.h"

/* Type Definitions */

/* Function Declarations */
static void multiGaussListND_mexCode_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);
MEXFUNCTION_LINKAGE mxArray *emlrtMexFcnProperties(void);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "multiGaussListND_mexCode", NULL, false, {2045744189U,2170104910U,2743257031U,4284093946U}, NULL };
emlrtCTX emlrtRootTLSGlobal = NULL;

/* Function Definitions */
static void multiGaussListND_mexCode_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *outputs[1];
  mxArray *inputs[3];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  int nInputs = nrhs;
  /* Check for proper number of arguments. */
  if(nrhs != 3) {
    mexErrMsgIdAndTxt("emlcoder:emlmex:WrongNumberOfInputs","3 inputs required for entry-point 'multiGaussListND_mexCode'.");
  } else if(nlhs > 1) {
    mexErrMsgIdAndTxt("emlcoder:emlmex:TooManyOutputArguments","Too many output arguments for entry-point 'multiGaussListND_mexCode'.");
  }
  /* Temporary copy for mex inputs. */
  for (n = 0; n < nInputs; ++n) {
    inputs[n] = (mxArray *)prhs[n];
  }
  /* Module initialization. */
  multiGaussListND_mexCode_initialize(&emlrtContextGlobal);
  /* Call the function. */
  multiGaussListND_mexCode_api((const mxArray**)inputs,(const mxArray**)outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  multiGaussListND_mexCode_terminate();
}

void multiGaussListND_mexCode_atexit_wrapper(void)
{
  multiGaussListND_mexCode_atexit();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(multiGaussListND_mexCode_atexit_wrapper);
  /* Dispatch the entry-point. */
  multiGaussListND_mexCode_mexFunction(nlhs, plhs, nrhs, prhs);
}

mxArray *emlrtMexFcnProperties(void)
{
    const char *mexProperties[] = {
        "Version",
        "EntryPoints"};
    const char *epProperties[] = {
        "Name",
        "NumberOfInputs",
        "NumberOfOutputs",
        "ConstantInputs"};
    mxArray *xResult = mxCreateStructMatrix(1,1,2,mexProperties);
    mxArray *xEntryPoints = mxCreateStructMatrix(1,1,4,epProperties);
    mxArray *xInputs = NULL;
    xInputs = mxCreateLogicalMatrix(1, 3);
    mxSetFieldByNumber(xEntryPoints, 0, 0, mxCreateString("multiGaussListND_mexCode"));
    mxSetFieldByNumber(xEntryPoints, 0, 1, mxCreateDoubleScalar(3));
    mxSetFieldByNumber(xEntryPoints, 0, 2, mxCreateDoubleScalar(1));
    mxSetFieldByNumber(xEntryPoints, 0, 3, xInputs);
    mxSetFieldByNumber(xResult, 0, 0, mxCreateString("8.0.0.783 (R2012b)"));
    mxSetFieldByNumber(xResult, 0, 1, xEntryPoints);

    return xResult;
}
/* End of code generation (multiGaussListND_mexCode_mex.c) */
