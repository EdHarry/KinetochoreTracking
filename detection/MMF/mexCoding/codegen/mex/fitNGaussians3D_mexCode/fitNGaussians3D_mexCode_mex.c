/*
 * fitNGaussians3D_mexCode_mex.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode'
 *
 * C source code generated on: Tue Nov 19 11:16:20 2013
 *
 */

/* Include files */
#include "mex.h"
#include "fitNGaussians3D_mexCode_api.h"
#include "fitNGaussians3D_mexCode_initialize.h"
#include "fitNGaussians3D_mexCode_terminate.h"

/* Function Declarations */
static void fitNGaussians3D_mexCode_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "fitNGaussians3D_mexCode", NULL, false, {2045744189U,2170104910U,2743257031U,4284093946U}, NULL };
void *emlrtRootTLSGlobal = NULL;

/* Function Definitions */
static void fitNGaussians3D_mexCode_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  mxArray *outputs[2];
  mxArray *inputs[4];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  int nInputs = nrhs;
  /* Module initialization. */
  fitNGaussians3D_mexCode_initialize(&emlrtContextGlobal);
  /* Check for proper number of arguments. */
  if (nrhs != 4) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:WrongNumberOfInputs", 5, mxINT32_CLASS, 4, mxCHAR_CLASS, 23, "fitNGaussians3D_mexCode");
  } else if (nlhs > 2) {
    emlrtErrMsgIdAndTxt(emlrtRootTLSGlobal, "EMLRT:runTime:TooManyOutputArguments", 3, mxCHAR_CLASS, 23, "fitNGaussians3D_mexCode");
  }
  /* Temporary copy for mex inputs. */
  for (n = 0; n < nInputs; ++n) {
    inputs[n] = (mxArray *)prhs[n];
  }
  /* Call the function. */
  fitNGaussians3D_mexCode_api((const mxArray**)inputs, (const mxArray**)outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  fitNGaussians3D_mexCode_terminate();
}

void fitNGaussians3D_mexCode_atexit_wrapper(void)
{
   fitNGaussians3D_mexCode_atexit();
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(fitNGaussians3D_mexCode_atexit_wrapper);
  /* Dispatch the entry-point. */
  fitNGaussians3D_mexCode_mexFunction(nlhs, plhs, nrhs, prhs);
}
/* End of code generation (fitNGaussians3D_mexCode_mex.c) */
