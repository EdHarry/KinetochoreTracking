/* MEX file for the Bresenham algorithm */

/* include to communicate with MatLab */
#include "mex.h"
#include "matrix.h"

/* other header files required */
#include <math.h>
#include <memory.h>
#include "tools.h"
#include "macros.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
/* mexFunction is the gateway routine for the MEX-file. 
   
   LHS
   ---
   plhs[0]: 2xn matrix with integer line points

   RHS
   ---
   prhs[0] : 2x1 matrix with integer start point coordinates
   prhs[1] : 2x1 matrix with integer end point coordinates
   */
{
  double *dPtr;
  int* x1;
  int* x2;
  int xStart[2];
  int xEnd[2];
  int i,n;

  dPtr = (double*)mxGetPr(prhs[0]);
  xStart[0] = DROUND(dPtr[0]);
  xStart[1] = DROUND(dPtr[1]);
  dPtr = (double*)mxGetPr(prhs[1]);
  xEnd[0] = DROUND(dPtr[0]);
  xEnd[1] = DROUND(dPtr[1]);
  
  n = tbx_getline(xStart[0],xStart[1],xEnd[0],xEnd[1],&x1,&x2);

  if(n>0){
    plhs[0] = mxCreateDoubleMatrix(2,n,mxREAL);
    dPtr = mxGetPr(plhs[0]);
    for(i=0; i<n; i++){
      dPtr[2*i] = (double)x1[i];
      dPtr[2*i+1] = (double)x2[i];
    }
    free(x1);
    free(x2);
  }
  return;
}
