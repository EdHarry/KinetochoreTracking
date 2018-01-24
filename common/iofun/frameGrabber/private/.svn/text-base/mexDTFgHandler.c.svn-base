/* MEX file to handle DATA TRANSLATION framegrabber */

#include "dtFg.h"

/* #include "membox.h" (not used at the moment; Win32 memory managment 
 applied instead */

/* header files required for the Matlab connection */
#include "mex.h"

/* other header files required */
#include "macros.h"
#include <math.h>
#include <string.h>


/* global variable for this dll */
static DEVINFO dev;

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
/* mexFunction is the gateway routine for the MEX-file. 
   Called by either dtFgOpen.m, dtFgClose.m, dtFgPrint.m,
   dtFgGrab.m ;
   
   LHS depends on task
   -------------------
   task == ("open" || "close" || "print")
   nrhs = 2
   plhs[0] return scalar: 1 successful completion; 0 error occurred.
   plhs[1] return scalar: 1 if device is open, 0 if closed.

   task == ("grab")
   nrhs = 3
   plhs[0] return scalar: 1 successful completion; 0 error occurred.
   plhs[1] return scalar: 1 if device is open, 0 if closed.
   plhs[2] image with 4 bytes/pixel (RGBDummy) in the structure 4*dim_Hori*dim_Verti
   plhs[3] roi [ulpos1,ulpos2,width,height] 
           with respect to the total grabbing area
           width and height are measured as distance NOT # of pixels

   RHS depends on task
   -------------------
   nrhs = 1 : task string

   task == ("grab")
   optional argument prhs[1] defintion of roi [ulpos1,ulpos2,width,height]
   
*/ 
{
  char msg[OLC_MAX_STATUS_MESSAGE_SIZE];
  char *task;
  int buflen;
  int status, dtStatusFlag;
  double *pStatus,*pOpen,*pdRoi;
  unsigned char* pLhsData;
  ULNG imgBufSize;
  int nDims;
  int* dims;
  int roi[4];

  /* create array for the status report */
  plhs[0] = mxCreateDoubleMatrix(1,1,mxREAL);
  pStatus = (double *)mxGetPr(plhs[0]);
  *pStatus = (double)0;
  /* create array for the open report */
  plhs[1] = mxCreateDoubleMatrix(1,1,mxREAL);
  pOpen = (double *)mxGetPr(plhs[1]);
  *pOpen = (double)0;


  /* get the task string */
  if ( mxIsChar(prhs[0])){
    buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0]) * sizeof(mxChar)) + 1;
    task = mxCalloc(buflen, sizeof(char));
    mxGetString(prhs[0], task, buflen);
#ifdef VERBOSE
    mexPrintf("The request is %s\n",task);
#endif
  }
  else
    return;

  /* interpretation of the task string */
  if (!strcmp(task,"grab")){
    /* create output array */
    nDims=3;
    dims = (int*) mxMalloc (nDims*sizeof(int));
    if(nrhs == 2){
      pdRoi = (double*)mxGetPr(prhs[1]);
      /* convert the coordinates to C standard coordinates */
      roi[0] = DROUND(pdRoi[0]-1);
      roi[1] = DROUND(pdRoi[1]-1);
      /* NOTE: the MatLab roi defintion is given by width and height 
	 while the roi defintion of the framgrabber is 
	 NUMBER of PIXELS in horizontal and vertical direction; therefore
	 the MatLab roi width and height have to be extended by 1 (GD 07-13-98)
      */
      roi[2] = DROUND(pdRoi[2]+1);
      roi[3] = DROUND(pdRoi[3]+1);
      /* clip to correct dimensions */
      roi[0] = roi[0]<(int)0?(int)0:roi[0];
      roi[1] = roi[1]<(int)0?(int)0:roi[1];
      roi[2] = (roi[2]+roi[0])>(int)dev.fmWidth?dev.fmWidth-roi[0]:roi[2];
      roi[3] = (roi[3]+roi[1])>(int)dev.fmHeight?dev.fmHeight-roi[1]:roi[3];
      /* set dimensions */
      dims[0]=dev.pixDepth;
      dims[1]=roi[2];
      dims[2]=roi[3];
    }
    else{
      roi[0]=roi[1]= 0;
      dims[0]=dev.pixDepth;
      dims[1]=roi[2]=dev.fmWidth;
      dims[2]=roi[3]=dev.fmHeight;
    }
    plhs[2]=mxCreateNumericArray(nDims,dims,mxUINT8_CLASS,mxREAL);
    mxFree(dims);
    plhs[3]=mxCreateDoubleMatrix(4,1,mxREAL);
    if(plhs[2] && plhs[3]){
      imgBufSize = roi[2]*roi[3]*dev.pixDepth;
      pLhsData = (unsigned char*)mxGetPr(plhs[2]);
      /* grab the image */
      status = grabFromDTFg((unsigned char*)pLhsData,roi,
			    imgBufSize*mxGetElementSize(plhs[2]),&dtStatusFlag);
      pdRoi = (double*)mxGetPr(plhs[3]);
      /* convert back to MatLab roi definition */
      pdRoi[0] = (double)(roi[0]+1);
      pdRoi[1] = (double)(roi[1]+1);
      pdRoi[2] = (double)(roi[2]-1);
      pdRoi[3] = (double)(roi[3]-1);
    }
    else{
      status = ALLOCATIONERROR;
    }
  }
  else{
    if( !strcmp(task,"open")){
      /* open the framegrabber device */
      status = openDTFg(&dtStatusFlag);
    }
    else{
      if ( !strcmp(task,"close")){
	status = closeDTFg(&dtStatusFlag);
      }
      else {
	if ( !strcmp(task,"print")){
	  printCtrlDTFg();
	  status = 1;
	}
      }
    }
  }
  *pOpen = (double)dev.isOpen;

  mxFree(task);
  if( status < 0 ){
    if( status == DTSPECIAL ){
      mexPrintf("DT SPECIAL error occurred:\n");
      OlImgGetStatusMessage(dtStatusFlag,msg,
			    OLC_MAX_STATUS_MESSAGE_SIZE*sizeof(char));
    }
    else{
      sprintf(msg,"DT error occurred: %d",status);
    }
    mexErrMsgTxt(msg);
  }
  
  *pStatus = (double)1;
}

int openDTFg(int *status)
/*
  INPUT 
  OUTPUT    status code which is evaluated in the above mexFunction
  RETURN  1 successful completion
         <0 error codes (see dtFg.h)
*/

{
  int nImgDev=0;
  int iImgDev;
  int slctImgDev;
  LPOLT_IMGDEVINFO devInfoList;

  /* get number of imaging devices on the system */
  if( (*status = OlImgGetDeviceCount(&nImgDev)) != OLC_STS_NORMAL)
    return(DTSPECIAL);

#ifdef VERBOSE
  mexPrintf("%d device(s) found\n",nImgDev);
#endif
  if ( nImgDev == 0 ){
    return(NODTDEVICEFOUND);
  }

  /* allocate the space for info */
  devInfoList = (OLT_IMGDEVINFO*)mxMalloc(nImgDev * sizeof(OLT_IMGDEVINFO));
  if ( !devInfoList){
    return(ALLOCATIONERROR); 
  }

  /* get device info */
  for(iImgDev = 0; iImgDev< nImgDev; iImgDev++)
    devInfoList[iImgDev].StructSize =  sizeof(OLT_IMGDEVINFO);
  if ((*status = OlImgGetDeviceInfo(devInfoList, 
				   nImgDev * sizeof(OLT_IMGDEVINFO))) 
      != OLC_STS_NORMAL){    
    mxFree(devInfoList);
    return(DTSPECIAL);
  }
#ifdef VERBOSE
  for(iImgDev = 0; iImgDev < nImgDev; iImgDev++){
    printImgDevInfoDTFg(&(devInfoList[iImgDev]));
  }
#endif

  /* if there are several devices attached to the system ... */
  if ( nImgDev > 1 ){
    /* ... select one manually */
    slctImgDev = selectDTFg(devInfoList, nImgDev);
  }
  else{
    slctImgDev = 0;
  }
#ifdef VERBOSE
  mexPrintf("%d selected\n",slctImgDev);
#endif
  
  /* physically open the selected device */
  if ((*status = OlImgOpenDevice(devInfoList[slctImgDev].Alias,
				 &(dev.devId)))
      != OLC_STS_NORMAL){
    dev.isOpen = 0;
    mxFree(devInfoList);
    return(DTSPECIAL);
  }
  
  dev.isOpen = 1;
  mxFree(devInfoList);

  /* get the framegrabber controls */
  if(!(getCtrlDTFg(status) > 0)){
    return(DTSPECIAL);
  }

  /* allocate the frame ID needed for subsequent grabbing */
  if((*status =OlFgAllocateBuiltInFrame ( dev.devId,
					  OLC_FG_DEV_MEM_VOLATILE,
					  OLC_FG_NEXT_FRAME,
					  &(dev.fmId)))
     !=  OLC_STS_NORMAL){ 
    return(DTSPECIAL);
  }

  return(1);
}

int closeDTFg(int *status)
/*
  INPUT 
  OUTPUT    status code which is evaluated in the above mexFunction
  RETURN  1 successful completion
  <0 error codes (see dtFg.h)
  */
{
  /* destroy the device frame frame */
  if((*status =OlFgDestroyFrame(dev.devId,
				dev.fmId))
     !=  OLC_STS_NORMAL){ 
    return(DTSPECIAL);
  }
  
  /* close the device */
  if(dev.isOpen){
    if ((*status = OlImgCloseDevice(dev.devId))
	!= OLC_STS_NORMAL){
      return(DTSPECIAL);
    }
    else{
      dev.isOpen = 0;
      return(1);
    }
  }
}

int selectDTFg(LPOLT_IMGDEVINFO pDevInfoList, int nImgDev)
{
  int iImgDev;
  double slct;
  mxArray* inputArray[1];
  mxArray* outputArray[1];

  /* create the inputArray holding the prompt string */
  inputArray[0] = 
    mxCreateString("Select framegrabber among the above devices : ");
  
  do{
    for(iImgDev = 0; iImgDev < nImgDev; iImgDev++){
      mexPrintf("%3d: %s %s\n",iImgDev+1,pDevInfoList[iImgDev].Alias,
		pDevInfoList[iImgDev].DeviceName);
    }
    /* call the matlab "input" function */ 
    mexCallMATLAB(1,outputArray,1,inputArray,"input");
    slct = mxGetScalar(outputArray[0]);
    mxDestroyArray(outputArray[0]);
  }while((slct<(double)1) || (slct>(double)(nImgDev)));
  mxDestroyArray(inputArray[0]);
  return((int)floor(slct-1));
}


int grabFromDTFg(unsigned char* buf,int roi[4],ULNG bufSize,int *status)
/* NOTE: buf must be allocated by the caller */
{

  if((*status =OlFgAcquireFrameToDevice ( dev.devId, dev.fmId ))
     !=  OLC_STS_NORMAL){ 
    return(DTSPECIAL);
  }

  if((*status =OlFgReadFrameRect( dev.devId, dev.fmId,
				  roi[0],roi[1],roi[2],roi[3],
				  buf,bufSize ))
     !=  OLC_STS_NORMAL){ 
    return(DTSPECIAL);
  }

  return(1);
}

int getCtrlDTFg(int *status)
{
  ULNG nInpSrc;
  char devName[OLC_MAX_DEVICE_NAME_STR_SIZE + 1];

  /* get device name */
  if((*status = OlImgQueryDeviceCaps(dev.devId,
				     OLC_IMG_DC_DEVICE_NAME,
				     dev.devName,
				     (OLC_MAX_DEVICE_NAME_STR_SIZE + 1)*sizeof(char)))
     != OLC_STS_NORMAL){ 
    return(DTSPECIAL);
  }
  
  /* get input source info */
  if((*status = OlFgQueryInputCaps(dev.devId,
				   OLC_FG_IC_INPUT_SOURCE_COUNT,
				   &(dev.nInpSrc),
				   sizeof(ULNG)))
     !=  OLC_STS_NORMAL){ 
    return(DTSPECIAL);
  }
  if ((*status = OlFgQueryInputVideoSource(dev.devId, 
					   &(dev.inpSrc)))
      !=  OLC_STS_NORMAL){ 
    return(DTSPECIAL);
  }

  /* get input frame type info */
  if((*status = OlFgQueryInputCaps(dev.devId,
				   OLC_FG_IC_FRAME_TYPE_LIMITS,
				   &(dev.fmTypeCaps),
				   sizeof(ULNG)))
     != OLC_STS_NORMAL){ 
    return(DTSPECIAL);
  }
  if ((*status = OlFgQueryInputControlValue(dev.devId,
					    dev.inpSrc,
					    OLC_FG_CTL_FRAME_TYPE,
					    &(dev.fmType)))
      != OLC_STS_NORMAL){ 
    return(DTSPECIAL);
  }
  
  /* get pixel depth info */
  if((*status = OlFgQueryInputCaps(dev.devId,
				   OLC_FG_IC_PIXEL_DEPTH,
				   &(dev.pixDepth),
				   sizeof(ULNG)))
     != OLC_STS_NORMAL){ 
    return(DTSPECIAL);
  }
 
  /* get frame dim info */
  if ((*status = OlFgQueryInputControlValue(dev.devId,
					    dev.inpSrc,
					    OLC_FG_CTL_FRAME_HEIGHT,
					    &(dev.fmHeight)))
      !=  OLC_STS_NORMAL){ 
    return(DTSPECIAL);
  }
  if ((*status = OlFgQueryInputControlValue(dev.devId,
					    dev.inpSrc,
					    OLC_FG_CTL_FRAME_WIDTH,
					    &(dev.fmWidth)))
      != OLC_STS_NORMAL){ 
    return(DTSPECIAL);
  }
  return(1);
}

void printCtrlDTFg(void)
{
  mexPrintf("-------------------------------------\n");
  mexPrintf("Name (open):  %s (%d)\n",dev.devName,dev.isOpen);
  mexPrintf("Nof Video Sources (current): %3d(%d)\n",
	    dev.nInpSrc,dev.inpSrc);
  mexPrintf("Frame Type Caps (set):  FRAME EVEN %1d (%1d)\n",
	    (dev.fmTypeCaps & OLC_FG_FRM_IL_FRAME_EVEN) > 0,
	    (dev.fmType & OLC_FG_FRM_IL_FRAME_EVEN) > 0);
  mexPrintf("                        FRAME ODD  %1d (%1d)\n",
	    (dev.fmTypeCaps & OLC_FG_FRM_IL_FRAME_ODD) > 0,
	    (dev.fmType & OLC_FG_FRM_IL_FRAME_ODD) > 0 );
  mexPrintf("                        FRAME NEXT %1d (%1d)\n",
	    (dev.fmTypeCaps & OLC_FG_FRM_IL_FRAME_NEXT) > 0,
	    (dev.fmType & OLC_FG_FRM_IL_FRAME_NEXT) > 0);
  mexPrintf("                        FIELD EVEN %1d (%1d)\n",
	    (dev.fmTypeCaps & OLC_FG_FRM_FIELD_EVEN) > 0,
	    (dev.fmType & OLC_FG_FRM_FIELD_EVEN) > 0 );
  mexPrintf("                        FIELD ODD  %1d (%1d)\n",
	    (dev.fmTypeCaps & OLC_FG_FRM_FIELD_ODD) > 0 ,
	    (dev.fmType & OLC_FG_FRM_FIELD_ODD) > 0 );
  mexPrintf("                        FIELD NEXT %1d (%1d)\n",
	    (dev.fmTypeCaps & OLC_FG_FRM_FIELD_NEXT) > 0,
	    (dev.fmType & OLC_FG_FRM_FIELD_NEXT) > 0);
  mexPrintf("                        NON INTERL %1d (%1d)\n",
	    (dev.fmTypeCaps & OLC_FG_FRM_NON_INTERLACED) > 0 ,
	    (dev.fmType & OLC_FG_FRM_NON_INTERLACED) > 0);
  mexPrintf("Pixel Depth [bytes]: %3d\n",dev.pixDepth);
  mexPrintf("Frame Width x Height: %3d x %3d\n",
	    dev.fmWidth,dev.fmHeight);
  mexPrintf("-------------------------------------\n");
}

void printImgDevInfoDTFg(LPOLT_IMGDEVINFO data)
{
  mexPrintf("--OLT_IMGDEVINFO-START---------------\n");
  mexPrintf("Struct Size [bytes]: %3d\n",data->StructSize);
  mexPrintf("Device Type: ");
  if(data->DeviceType & OLC_IMG_DEV_BOGUS)
    mexPrintf("BOGUS\n");
  if(data->DeviceType & OLC_IMG_DEV_MONO_FRAME_GRABBER)
    mexPrintf("MONO_FRAME_GRABBER\n");
  if(data->DeviceType & OLC_IMG_DEV_COLOR_FRAME_GRABBER)
    mexPrintf("COLOR_FRAME_GRABBER\n");
  mexPrintf("Device Id: %3d\n",data->DeviceId);
  mexPrintf("Device Alias: %s\n",data->Alias);
  mexPrintf("Device Name: %s\n",data->DeviceName);
  mexPrintf("Host Dev Mem Size [bytes]: %8d\n",data->HostDeviceMemSize);
  mexPrintf("Host Lin Mem Size [bytes]: %8d\n",data->HostLinearMemSize);
  mexPrintf("--OLT_IMGDEVINFO-END-----------------\n");
}


