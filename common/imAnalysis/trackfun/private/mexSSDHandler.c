/* MEX file to handle SSD the SSD tracking */

/* include to communicate with MatLab */
#include "mex.h"
#include "matrix.h"

/* other header files required */
#include <math.h>
#include <memory.h>
#include "macros.h"
#include "dLsm.h"
#include "mexSSDHandler.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
/* mexFunction is the gateway routine for the MEX-file. 
   Called by either
   
   LHS depends on task
   -------------------
   plhs[0]: scalar with  status: 1 successful completion; 0 error occurred;]
            this status considers only the interface functions;
            SSD errors must be checked in the return structure

   plhs[1]: return structure with the following fields:
  
   "status": DLSMMAINSTRUCT.status
   "pos":    final position of the patch in MatLab coordinates
   "lori":   scalar value with the local orientation of the template texture
   "shape":  2x2 matrix with the shape parameters [mx sx; sy my]
   "preci": 6x1 matrix with precision of the estimated parameters
   "goffit": 2x1 goodness-of-fit values [sigma0, cross_correl_coefficient]
   "patchStack": stack of patch series; only filled if prhs[0].opts(6) == 1
   "tImg" : updated (in case of prefiltering) version of template image;
   "sImg" : updated (in case of prefiltering) version of search image;

   RHS
   ---
   prhs[0] command structures with the following fields:
   
   "fname" filename with the default values
   "tImg" template image
   "tMsk" template mask
   "sImg" search image
   "pSze" patch size [width, height] (in MatLab definition; will be converted
                                      to # of pixels in horizontal and vertical
                                      direction)
   "tPos" position of template center; tPos(1) and tPos(2) in MatLab define 
          the center with respect to a left handed coordinate frame with 
          origin 1,1 in the upper left corner of the image.
   "pos0" initial guess of the patch position in the search image; 
          tPos(1) and tPos(2) in MatLab define 
          the center with respect to a left handed coordinate frame with 
          origin 1,1 in the upper left corner of the image.
   "shape0"  initial guess of the shape matrix (optional: if either the field 
             is missing or [], then the matching starts with an identity 
	     transformation [mx = 1.0 sx = 0.0; sy = 0.0 my = 1.0] as a 
             default 
   "maxEx" maximal excursion of the patch in either displacement or deformation
   "resol" resolution where resol(1) is the resolution for displacement and 
           resol(2) the multiplication factor to get the resolution in shape 
           measurements.
   "cParams" constraint parameters (obsolete if opts(3) is an unconstrained
             parameter set)
   "opts" options:
          (1): computation mode -> 0 = no statistics; 1 = statistics; 
                                   2 = diagnostics;
          (2): interpolation method -> 1 = nearest neighbour; 2 = bilinear;
          (3): parameter set -> 1 = all parameters
                                2 = shift only
                                3 = shift + scales
                                4 = congruent
				5 = shift + rotation 
				6 = all parameters , line constrained
				7 = shift , line constrained
				8 = shift + scales, line constrained
				9 = congruent, line constrained
			       10 = shift + rotation , line constrained
          (4): prefiltering of template image (1 = yes; 0 = no)
          (5): prefiltering of search image (1 = yes; 0 = no)
          (6): patch series stored and returned (1 = yes; 0 = no)
          (7): verbose (1 = yes; 0 = no)
         
   "filtParam": this field is read only if either opts(4) or opts(5) is set to 1.
                contains the sigma parameter for Gaussian prefiltering
   
   
*/ 
{
  double* pStatus;
  const char* fieldNames[] = {"status","pos","lori","shape","preci","goffit",
			      "patchStack","tImg","sImg"};
  int     nFields = 9;
  int     nCells,nTPos,nPos0,nOpts,nResol,nShape0,mShape0;
  char*   fname;
  int     bufSze;
  DLSMIMAGE   tImg,sImg;
  DLSMMAINSTRUCT* dLsmData = NULL;
  unsigned char* mask;
  int     nx,ny; 
  int     pDim[2];
  int     maxEx;
  double* tPos;
  double* pos0;
  double* shape0 = NULL;
  double* resol;
  double* cParams;
  int     nCparams;
  double  sigma = (double)0.0;
  int*    opts;
  int     compMode;
  int     paramSet[6]=FULL;
  int     constraints = CNONE;
  char*   msg;
  double  auxD;
#ifdef VERBOSE
  int     i;
#endif

  /* create array for the status report */
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  pStatus = (double *)mxGetPr(plhs[0]);
  *pStatus = (double)0;

  /* create the result structure */
  plhs[1] = mxCreateStructMatrix(1,1,nFields,fieldNames);
   
  /* some security checks on the command structure */
  nCells = mxGetN(prhs[0])*mxGetM(prhs[0]);
  if((nCells != 1) || !mxIsStruct(prhs[0]))
    mexErrMsgTxt("Invalid input command structure");

  /* convert the filename with the defaults */
  string2charBuf(&fname,&bufSze,mxGetField(prhs[0],0,"fname"));

  /* convert images */
  if(!mat2dlsmImg(&tImg,mxGetField(prhs[0],0,"tImg")))
    mexErrMsgTxt("Invalid template image");
  if(!mat2dlsmImg(&sImg,mxGetField(prhs[0],0,"sImg")))
    mexErrMsgTxt("Invalid search image");

  /* convert patch dimensions */
  if(!wh2nm(pDim,mxGetField(prhs[0],0,"pSze")))
    mexErrMsgTxt("Invalid patch dimensions");


  /* extract mask if available */
  if(!mat2uCharMtx(&mask,&nx,&ny,mxGetField(prhs[0],0,"tMsk")))
    mask = NULL;
  else{
    if((nx!=pDim[0]) || (ny!=pDim[1])){
      free(mask);
      mask = NULL;
    }
  }

  /* convert border size */
  if(!mat2intScalar(&maxEx,mxGetField(prhs[0],0,"maxEx")))
    mexErrMsgTxt("Invalid excursion limit option");
#ifdef VERBOSE  
  mexPrintf("Tolerable excursion: %d\n",maxEx);
#endif

  /* convert the resolution vector */
  if(!mat2dblVec(&resol,&nResol,mxGetField(prhs[0],0,"resol")))
    mexErrMsgTxt("Invalid pos0 field");
  if(nResol != 2 )  mexErrMsgTxt("Invalid resol field");
#ifdef VERBOSE
  mexPrintf("Computational resolution : %8.3lf %8.3lf\n",resol[0],resol[1]);
#endif

  /* convert constraints' parameter vector */
  if(!mat2dblVec(&cParams,&nCparams,mxGetField(prhs[0],0,"cParams")))
    cParams = NULL;

  /* convert positional vectors and the inital shape, if available */
  if(!mat2dblVec(&pos0,&nPos0,mxGetField(prhs[0],0,"pos0")))
    mexErrMsgTxt("Invalid pos0 field");
  if(nPos0 != 2 )  mexErrMsgTxt("Invalid pos0 field");

  if(mat2dblMtx(&shape0,&nShape0,&mShape0,mxGetField(prhs[0],0,"shape0"))){
    /* transpose the shape matrix to a C-style matrix */
    if(shape0){
      if((nShape0 == 2) && (mShape0 == 2)){
	/* swap the elements shape0[1] and shape0[2] */
	auxD = shape0[1];
	shape0[1] = shape0[2];
	shape0[2] = auxD;
      }
    }
  }

  if(!mat2dblVec(&tPos,&nTPos,mxGetField(prhs[0],0,"tPos")))
    mexErrMsgTxt("Invalid tPos field");
  if(nTPos != 2 )  mexErrMsgTxt("Invalid tPos field");
#ifdef VERBOSE
  mexPrintf("Template position : %8.3lf %8.3lf\n",tPos[0],tPos[1]);
  mexPrintf("Initial position  : %8.3lf %8.3lf\n",pos0[0],pos0[1]); 
#endif

  /* convert option flags */
  if(!mat2intVec(&opts,&nOpts,mxGetField(prhs[0],0,"opts")))
    mexErrMsgTxt("Invalid opts field");
  if(nOpts != 7 )  mexErrMsgTxt("Invalid opts field");
#ifdef VERBOSE
  mexPrintf("Options set: ");
  for(i=0; i<nOpts; i++) mexPrintf("%d ",opts[i]);
  mexPrintf("\n");
#endif
  switch(opts[0]){
  case 0: {compMode = NOSTATISTICS; break;}
  case 1: {compMode = STATISTICS; break;}
  case 2: {compMode = DIAGNOSTICS; break;}
  }
  switch(opts[2]){
  case 2: {memset(paramSet+2,0,4*sizeof(int)); constraints = CNONE; break;}
  case 3: {memset(paramSet+3,0,2*sizeof(int)); constraints = CNONE; break;}
  case 4: {constraints = CCONGRUENT; break;}  
  case 5: {constraints = CCONGRUENT | CROT; break;} 
  case 6: {constraints = CLINE; break;}
  case 7: {memset(paramSet+2,0,4*sizeof(int)); constraints = CLINE; break;}
  case 8: {memset(paramSet+3,0,2*sizeof(int)); constraints = CLINE; break;}
  case 9: {constraints = CCONGRUENT | CLINE; break;}  
  case 10:{constraints = CCONGRUENT | CLINE | CROT; break;} 
  }

  /* get filter parameter if necessary */
  if(opts[3] || opts[4] ){
    if(!mat2dblScalar(&sigma,mxGetField(prhs[0],0,"filtParam")))
      mexErrMsgTxt("Invalid filter parameter option");
  }

  /* start the matching */
  if(!dLsmInitStruct(&dLsmData,fname,&tImg,&sImg,pDim[0],pDim[1],maxEx,
		     tPos[0]-1,tPos[1]-1,mask,resol[0],resol[1],
		     opts[1]|(SAVEPATCHES*opts[5]))){
    packResult(plhs[1],dLsmData);
    return;
  }     

  /* prefilter if requested */
  if(!dLsmPrefilter(dLsmData,sigma,
		    NOIMAGE|(opts[3]*TEMPLATEIMAGE)|(opts[4]*PATCHIMAGE),
		    GAUSSIAN)){
    packResult(plhs[1],dLsmData);
    return;
  }

  /* initialize the parameters */
  if(!dLsmInitParams(dLsmData,paramSet,constraints,cParams,compMode)){
    packResult(plhs[1],dLsmData);
    return;
  } 

  /* initialize the patch */
  if(!dLsmInitPatch(dLsmData,pos0[0]-1,pos0[1]-1,shape0)){
    packResult(plhs[1],dLsmData);
    return;
  } 

  /* compute the solution */
  if(!dLsmCompute(dLsmData)){
    packResult(plhs[1],dLsmData);
    return;
  }

  /* print results in verbose mode */
  if(opts[6]){
    msg = (char*)mxCalloc(RESULTBUFNCHAR,sizeof(char));
    dLsmPrintResults(dLsmData,msg,RESULTBUFNCHAR*sizeof(char));
    mexPrintf("%s",msg);
    if(dLsmData->statistics.exist){
      dLsmPrintStatistics(dLsmData,msg,RESULTBUFNCHAR*sizeof(char));
      mexPrintf("%s",msg);
    }
    mxFree(msg);
  }

  /* pack the results */
  packResult(plhs[1],dLsmData);

  /* free the stuff which has been allocated during field conversion */
  mxFree(opts);
  if(fname) 
    mxFree(fname);
  if(mask)
    mxFree(mask);
  if(shape0)
    mxFree(shape0);
  mxFree(tPos);
  mxFree(pos0);
  mxFree(resol);
  if(cParams)
    mxFree(cParams);
  mxFree(sImg.data);
  mxFree(tImg.data);
}


int mat2dlsmImg(DLSMIMAGE* img,mxArray* mat)
{
  int status;

  status = mat2dblMtx(&(img->data),&(img->dim_x),&(img->dim_y),mat);

  return(status);
  /*
  int  nDims;
  const int* dims;
  double*  dPtr;
  
  img->dim_x = img->dim_y = 0;

  if(mxIsEmpty(mat))
    return(0);
  nDims = mxGetNumberOfDimensions(mat);
  if(nDims != 2)
    return(0);
  dims = mxGetDimensions(mat);
  if( (dims[0] <= 0) || (dims[1] <=0))
    return(0);
  if(!mxIsDouble(mat))
    return(0);

  img->dim_x = dims[0];
  img->dim_y = dims[1];
  img->data = (double*)mxMalloc(dims[0]*dims[1]*sizeof(double));
  dPtr = mxGetPr(mat);
  memcpy(img->data,dPtr,dims[0]*dims[1]*sizeof(double));

#ifdef VERBOSE
  mexPrintf("%d x %d image converted\n",img->dim_x,img->dim_y);
#endif
  return(1);

  */
}

mxArray* dlsmImg2mat(DLSMIMAGE* img)
{

  mxArray* mat;

  mat = dblMtx2mat(img->data,img->dim_x,img->dim_y);

  return(mat);
}


int wh2nm(int* nm,mxArray* mat)
{
  double* pData;

  nm[0] = nm[1] =0;

  if(mxIsEmpty(mat))
    return(0);
  if((mxGetN(mat)*mxGetM(mat)) != 2)
    return(0);
  
  pData = mxGetPr(mat);

  nm[0] = DROUND(pData[0])+1;
  nm[1] = DROUND(pData[1])+1;

#ifdef VERBOSE
  mexPrintf("%d x %d entered\n",nm[0],nm[1]);
#endif
  return(1);
}


void packResult(mxArray* rep, DLSMMAINSTRUCT* data)
     /* packs the results and frees the structure 
	rep is an 1x1 array of structures with field allocated in the 
	main mex routine
	*/
{
  mxArray* aux;
  double*  dPtr;
  double*  inStack;
  double*  inPreci;
  int      nPatchEl;
  int      i;
  int nDims[3];

  if(!data) return;

  /* NOTE all the aux, which are pointers to mxArray structures are created
     with mxCreate....(). MATLAB has memory infrastructure which guarantees
     proper freeing of memory when leaving the function */

  /* set status field */
  aux=intScalar2mat(data->status);
  mxSetField(rep,0,"status",aux);

  /* return updated versions of template and patch image */
  aux=dlsmImg2mat(&(data->tmpl_image));
  mxSetField(rep,0,"tImg",aux);
  aux = dlsmImg2mat(&(data->patch_image));
  mxSetField(rep,0,"sImg",aux);

  /* fill position field */
  aux=mxCreateDoubleMatrix(1,2,mxREAL);
  dPtr = mxGetPr(aux);
  dLsmGetResults(data,dPtr,mxGetNumberOfElements(aux)*mxGetElementSize(aux));
  /* set the position back to MatLab defined coordinates */
  dPtr[0] += (double)1;
  dPtr[1] += (double)1;
  mxSetField(rep,0,"pos",aux);

  /* fill the local orientation field, if computed */
  if(data->tmpl_lOri)
    aux=dblScalar2mat(*(data->tmpl_lOri));
  else
    /* create empty matrix */
    aux=dblMtx2mat(NULL,0,0);
  mxSetField(rep,0,"lori",aux);
  
  /* fill the shape field 
     it is directly transposed here */
  if((data->status>0)||(data->status<=LSM_NOPARAMSTRUCT)){
    aux=mxCreateDoubleMatrix(2,2,mxREAL);
    dPtr = mxGetPr(aux);
    dPtr[0] = data->result.paramsup[2]; /* mx */
    dPtr[2] = data->result.paramsup[3]; /* sx */
    dPtr[1] = data->result.paramsup[4]; /* sy */
    dPtr[3] = data->result.paramsup[5]; /* my */
    mxSetField(rep,0,"shape",aux);
  }

  /* fill precision */
  if(data->statistics.exist){
    aux = mxCreateDoubleMatrix(1,6,mxREAL);
    dPtr = mxGetPr(aux);
    memset(dPtr,0,6*sizeof(double));
    inPreci = data->statistics.sparams;
    for(i=0; i<6; i++){
      if(data->result.paramset[i]){
	dPtr[i]=*(inPreci++);
      }
    }
    mxSetField(rep,0,"preci",aux);
  }

  /* fill "goffit field */
  if(data->statistics.exist){
    aux = mxCreateDoubleMatrix(1,2,mxREAL);
    dPtr = mxGetPr(aux);
    dPtr[0] = data->statistics.s0;
    dPtr[1] = data->statistics.c;
    mxSetField(rep,0,"goffit",aux);
  }

  /* fill patch stack */
  if(data->savepatches && data->patchseries){
    if(data->patchseries[0].data){
      nDims[0] = data->patchseries[0].dim_x;
      nDims[1] = data->patchseries[0].dim_y;
      nDims[2] = data->nIter?data->nIter:1;
      aux = mxCreateNumericArray(3,nDims,mxDOUBLE_CLASS,mxREAL); 
      dPtr = mxGetPr(aux);
      memset(dPtr,0,mxGetNumberOfElements(aux)*mxGetElementSize(aux));
      inStack = dPtr;
      nPatchEl = data->patchseries[0].dim_x*data->patchseries[0].dim_y;
      for(i = 0; i<data->nIter; i++){
	if(data->patchseries[i].data){
	  memcpy(inStack,data->patchseries[i].data,
		 nPatchEl*sizeof(double));
	  inStack += nPatchEl;
	}
      }
      mxSetField(rep,0,"patchStack",aux);
    }
  }

  dLsmFreeStruct(&data);
  return;
}




