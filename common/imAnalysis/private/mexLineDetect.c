/* MEX file for line filtering / detection */

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
   plhs[0]: "lineness" (filter response)
   plhs[1]: orientation [0,M_PI], perpendicular to the local axis direction
   plhs[2]: details, structure with the following fields:
            "status" status returned from tbx_MsLineFilter()
            "hessian" local orientation computed via the Hessian 
                      (sometimes obsolete)
	    "scaleSelect" uint8 buffer with selected scales, encoded as 
	                  scaleIndex+1;
	    "noise"   noise value based on which all the internal threshholds
	              are set.

	     
   RHS
   ---
   prhs[0] : image matrix (alraedy in "C-permutation")
   prhs[1] : options, command structure with the following fields:
             "scales" : vector with the scales for multi-scale filtering
	     "noise"  : estimate of the image noise; if <0 an internal value
	                is computed 
	     "linetype": definition of the linetype sought;
	                 is constructed as bitwise ored of the following types
                         1 positive line (bright line on dark background)
                         2 negative line (dark line on bright background)
			 4 wave line (line with bright and dark stripe on 
                           grey background)
			 therefore:
                         3 looks for positive AND negative lines
                         5 looks for positive AND wave lines
                         6 looks for negative AND wave lines
                         7 looks for ALL three types
             "nonMaxSupp" non maximum suppression yes = 1; no = 0;
    */
{
  const char* fieldNames[] = {"status","hessian","scaleSelect","noise","maxMap"};
  int     nFields = 5;
  int     status;
  double *img,*scales;
  double  noise;
  int     nx,ny,nScales,linetype,nonMaxSupp;
  double *outAngleCanny,*outAngleHessian,*outResp;
  unsigned char *outWhichScales,*maxMap = NULL;
  mxArray*   aux;

  /* unpack image */
  if(!mat2dblMtx(&img,&nx,&ny,prhs[0]))
    mexErrMsgTxt("Invalid input image");


  /* allocate the buffers needed */
  outAngleCanny = (double*)mxMalloc(nx*ny*sizeof(double));
  outAngleHessian = (double*)mxMalloc(nx*ny*sizeof(double));
  outResp = (double*)mxMalloc(nx*ny*sizeof(double));
  outWhichScales = (unsigned char*)mxMalloc(nx*ny*sizeof(unsigned char));

  /* unpack the option fields */
  if(!mat2dblVec(&scales,&nScales,mxGetField(prhs[1],0,"scales")))
    mexErrMsgTxt("Invalid vector with line scales");
  
  if(!mat2dblScalar(&noise,mxGetField(prhs[1],0,"noise")))
      mexErrMsgTxt("Invalid noise value");

  if(!mat2intScalar(&linetype,mxGetField(prhs[1],0,"linetype")))
    mexErrMsgTxt("Invalid linetype selection");
  
  if(!mat2intScalar(&nonMaxSupp,mxGetField(prhs[1],0,"nonMaxSupp")))
    mexErrMsgTxt("Invalid non-maximum suppression option");

  if(nonMaxSupp)
    maxMap = (unsigned char*)mxMalloc(nx*ny*sizeof(unsigned char));

  status = tbx_dMsLineFilter(img,NULL,0,scales,nScales,
			     outAngleCanny,
			     outAngleHessian,
			     outResp,
			     outWhichScales,
			     &noise,nx,ny,linetype,maxMap);


  /* fill result buffers */
  plhs[0] = dblMtx2mat(outResp,nx,ny);
  plhs[1] = dblMtx2mat(outAngleCanny,nx,ny);

  /* create return structure with all the details */
  plhs[2] = mxCreateStructMatrix(1,1,nFields,fieldNames);
  
  /* fill status */
  aux = intScalar2mat(status);
  mxSetField(plhs[2],0,"status",aux);
  
  /* fill Hessian */
  aux = dblMtx2mat(outAngleHessian,nx,ny);
  mxSetField(plhs[2],0,"hessian",aux);
  
  /* fill selected scales array */
  aux = uCharMtx2mat(outWhichScales,nx,ny);
  mxSetField(plhs[2],0,"scaleSelect",aux);

  /* fill noise value */
  aux = dblScalar2mat(noise);
  mxSetField(plhs[2],0,"noise",aux);
  
   /* fill maxMap array */
  aux = uCharMtx2mat(maxMap,nx,ny);
  mxSetField(plhs[2],0,"maxMap",aux);

  /* free data buffers */
  mxFree(img);
  mxFree(scales);
  mxFree(outAngleHessian);
  mxFree(outAngleCanny);
  mxFree(outResp);
  mxFree(outWhichScales);
  if(maxMap)
    mxFree(maxMap);

  return;
}
