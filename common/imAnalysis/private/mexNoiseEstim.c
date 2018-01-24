/* MEX file for noise estimation according to Vorhees and Poggio, 1987 */

/* include to communicate with MatLab */
#include "mex.h"
#include "matrix.h"

/* other header files required */
#include "mexUtils.h"
#include "tools.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
/* mexFunction is the gateway routine for the MEX-file. 
   
   LHS
   ---
   plhs[0]: scalar with noise estimate
   plhs[1]: scalar with threshold
   plhs[2]: nInGradBufx1 matrix with double valued gradient values

   RHS
   ---
   prhs[0] : nx x ny matrix with image
   prhs[1] : scalar with risk probability

   */
{
  int status;
  double *buf,*gradMagn = NULL;
  double c,noise,thresh;
  int     nx,ny,nVals;

  status = mat2dblMtx(&buf,&nx,&ny,prhs[0]);
  if(!status)
    return;

  status = mat2dblScalar(&c,prhs[1]);

  tbx_dNoiseEstim(buf,nx,ny,&gradMagn,&nVals,c,&noise,&thresh);

  plhs[2] = dblMtx2mat(gradMagn,nVals,1);

  if(gradMagn) 
    free(gradMagn);

  plhs[0]=dblScalar2mat(noise);
  plhs[1]=dblScalar2mat(thresh);

}
