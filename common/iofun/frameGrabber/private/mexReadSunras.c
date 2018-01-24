/* MEX file to read a SUN raster image */

#include "sunras.h"
#include "mex.h"

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[] )
/* mexFunction is the gateway routine for the MEX-file. 
   Called by ReadSunras() MATLAB function;
   input arguments required (tested in ReadSunras()):
   nlhs = 1 : outgoing uInt8 image
   nrhs = 1 : incoming filename string
*/ 
{
  char* fname;
  int   nChar,nX,nY;
  unsigned char* imBuf;

  nChar = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
  
  /* Allocate enough memory to hold the converted string. */ 
  if(!(fname = (char*)mxCalloc(nChar, sizeof(char)))){
    mexErrMsgTxt("Allocation failure.");
  }
  /* convert the String array to the char* buffer */
  if( mxGetString(prhs[0], fname, nChar)){
    mexErrMsgTxt("Could not convert string data.");
  }
  
  mexPrintf("The converted string is \n%s.\n", fname);
  
  /* read the image from file */
  readImage(&imBuf,fname,&nX,&nY);

  /* create the output matrix */
  
  /* copy the data from the image buffer to the array */

  /* free the memory allocated in readSunras() */
}

int readImage(unsigned char **image,const char *fname,int *nx,int *ny)
/*
  INPUT   fname : name of file to be read
  OUTPUT  image : buffer with image data (allocated in this procedure)
          nx,ny : size of image
  RETURN  1 : if successful completion
         <0 : in case of any reading error (see error codes in sunras.h) 
*/

{
  FILE  *fp;
  int header[8],i,i_size,ret =0;
  unsigned char *colormap=NULL,*rgb=NULL,tmp;
    
  if((fp = fopen(fname,"r"))== NULL){
    *image=NULL;
    *nx=*ny=0;   
    return(FILEOPENFAILED);
  }
  
  fread(header,sizeof(int),8,fp);
  *nx=header[1];
  *ny=header[2];
  i_size=(*nx)*(*ny);
  
  if (header[5] == RT_STANDARD){
    *image=(unsigned char *) malloc(i_size*sizeof(unsigned char));
    switch (header[3]){
    case 8:
      if ((header[6] == RMT_EQUAL_RGB)&&(header[7])){
	colormap=(unsigned char *) malloc(header[7]*sizeof(char));
	fread(colormap,1,header[7],fp);
      }
      for (i=0;i<*ny;i++){
	fread((*image+i*(*nx)),1,(*nx),fp);
	if ((*nx)%2) fread(&tmp,1,1,fp);
      }
      if (colormap) 
	for (i=0;i<i_size;i++) (*image)[i]=colormap[(*image)[i]];
      ret=1;
      break;
    case 24:
      if (header[6] == RMT_NONE){
	rgb=(unsigned char *) malloc(header[4]*sizeof(unsigned char));
	fread(rgb,1,header[4],fp);
	for (i=0;i<i_size;i++) 
	  (*image)[i]=(unsigned char)(rgb[i*3]*0.299+
				      rgb[i*3+1]*0.587+
				      rgb[i*3+2]*0.114);
      }
      ret=1;
      break;
    default: 
      ret = UNSUPPORTEDDEPTH;
    }
  }
  else{
    ret = NOTSTANDARDSUNRAS;
  }
  if(colormap) free(colormap);
  if(rgb) free(rgb);
  return(ret);
}
