/*
 * fitNGaussians3D_mexCode_F_mex.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode_F'
 *
 * C source code generated on: Mon May 21 17:33:36 2012
 *
 */

/* Include files */
#include "mex.h"
#include "fitNGaussians3D_mexCode_F_api.h"
#include "fitNGaussians3D_mexCode_F_initialize.h"
#include "fitNGaussians3D_mexCode_F_terminate.h"

/* Type Definitions */

/* Function Declarations */
static void fitNGaussians3D_mexCode_F_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "fitNGaussians3D_mexCode_F", NULL, false, NULL };

/* Function Definitions */
static void fitNGaussians3D_mexCode_F_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Temporary copy for mex outputs. */
  mxArray *outputs[1];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  /* Check for proper number of arguments. */
  if(nrhs != 4) {
    mexErrMsgIdAndTxt("emlcoder:emlmex:WrongNumberOfInputs","4 inputs required for entry-point 'fitNGaussians3D_mexCode_F'.");
  } else if(nlhs > 1) {
    mexErrMsgIdAndTxt("emlcoder:emlmex:TooManyOutputArguments","Too many output arguments for entry-point 'fitNGaussians3D_mexCode_F'.");
  }
  /* Module initialization. */
  fitNGaussians3D_mexCode_F_initialize(&emlrtContextGlobal);
  /* Call the function. */
  fitNGaussians3D_mexCode_F_api(prhs,(const mxArray**)outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  fitNGaussians3D_mexCode_F_terminate();
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(fitNGaussians3D_mexCode_F_atexit);
  emlrtClearAllocCount(&emlrtContextGlobal, 0, 0, NULL);
  /* Dispatch the entry-point. */
  fitNGaussians3D_mexCode_F_mexFunction(nlhs, plhs, nrhs, prhs);
}
/* End of code generation (fitNGaussians3D_mexCode_F_mex.c) */
