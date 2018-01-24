/*
 * objFcn_mex.c
 *
 * Code generation for function 'objFcn'
 *
 * C source code generated on: Fri May 25 21:48:51 2012
 *
 */

/* Include files */
#include "mex.h"
#include "objFcn_api.h"
#include "objFcn_initialize.h"
#include "objFcn_terminate.h"

/* Type Definitions */

/* Function Declarations */
static void objFcn_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]);

/* Variable Definitions */
emlrtContext emlrtContextGlobal = { true, false, EMLRT_VERSION_INFO, NULL, "objFcn", NULL, false, NULL };

/* Function Definitions */
static void objFcn_mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Temporary copy for mex outputs. */
  mxArray *outputs[4];
  int n = 0;
  int nOutputs = (nlhs < 1 ? 1 : nlhs);
  /* Check for proper number of arguments. */
  if(nrhs != 9) {
    mexErrMsgIdAndTxt("emlcoder:emlmex:WrongNumberOfInputs","9 inputs required for entry-point 'objFcn'.");
  } else if(nlhs > 4) {
    mexErrMsgIdAndTxt("emlcoder:emlmex:TooManyOutputArguments","Too many output arguments for entry-point 'objFcn'.");
  }
  /* Module initialization. */
  objFcn_initialize(&emlrtContextGlobal);
  /* Call the function. */
  objFcn_api(prhs,(const mxArray**)outputs);
  /* Copy over outputs to the caller. */
  for (n = 0; n < nOutputs; ++n) {
    plhs[n] = emlrtReturnArrayR2009a(outputs[n]);
  }
  /* Module finalization. */
  objFcn_terminate();
}
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
  /* Initialize the memory manager. */
  mexAtExit(objFcn_atexit);
  emlrtClearAllocCount(&emlrtContextGlobal, 0, 0, NULL);
  /* Dispatch the entry-point. */
  objFcn_mexFunction(nlhs, plhs, nrhs, prhs);
}
/* End of code generation (objFcn_mex.c) */
