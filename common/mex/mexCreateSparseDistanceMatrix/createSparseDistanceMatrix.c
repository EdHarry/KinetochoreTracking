
#include "matrix.h"
#include "mex.h"
#include <math.h>



#if defined(NAN_EQUALS_ZERO)
#define IsNonZero(d) ((d) != 0.0 || mxIsNaN(d))
#else
#define IsNonZero(d) ((d) != 0.0)
#endif


double percent_sparse=0.01;
double percent_step=0.01;

void calcDistMatrix1D(mxArray *plhs[], double *M, double *N, mwSize Mrows, mwSize Nrows, mwSize nzmax, 
					  double threshold, double eps)
{
    double dist;
	mwIndex  i,j,k;

	double *sr;
	mwSize *irs, *jcs;

	
	sr=mxGetPr(plhs[0]);
	irs=mxGetIr(plhs[0]);
	jcs=mxGetJc(plhs[0]);

    
	k=0;
	for (i=0;i<Nrows;i++) 
	{
        
		jcs[i]=k;

        for (j=0;j<Mrows;j++)
		{        
            
            dist=*(N+i)-*(M+j);            

			
			if (abs((int)dist)<=threshold)
			{
				if (dist==0) 
					dist=eps;

				
				if(k >= nzmax) 
				{
					
					int oldnzmax = nzmax;

					
					percent_sparse += percent_step;
					nzmax=(mwSize)ceil((double)Mrows*(double)Nrows*percent_sparse);

					
					if(nzmax <= oldnzmax)
						nzmax=oldnzmax+1;

					
					mxSetNzmax(plhs[0], nzmax);

					
					mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
					mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(mwIndex)));

					
					sr = mxGetPr(plhs[0]);
					irs = mxGetIr(plhs[0]);

				}
					
				
				sr[k]=dist;
				irs[k]=j;
				k++;
			}

        }	
	}
	jcs[Nrows]=k;

	
    

	if (k>0)
		nzmax=k;
	else
		nzmax=1;

	
	mxSetNzmax(plhs[0], nzmax);

	
	mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
	mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(mwIndex)));

	
	sr=mxGetPr(plhs[0]);
	irs=mxGetIr(plhs[0]);
	jcs=mxGetJc(plhs[0]);
}

void calcDistMatrix2D(mxArray *plhs[], double *M, double *N, mwSize Mrows, mwSize Nrows, mwSize nzmax, 
					  double threshold, double eps)
{
    double	mX=0, mY=0;
    double	nX=0, nY=0;
    double  *posM, *posN;
    double dist;
	mwIndex  i,j,k;

	
	double *sr;
	mwSize *irs, *jcs;
	

	sr=mxGetPr(plhs[0]);
	irs=mxGetIr(plhs[0]);
	jcs=mxGetJc(plhs[0]);

    
	k=0;
	for (i=0;i<Nrows;i++) 
	{
       
        posN=N+i;
        nX=*posN;
        nY=*(posN+Nrows);
        
		
		jcs[i]=k;

		for (j=0;j<Mrows;j++) 
		{        
           
            posM=M+j;
            mX=*posM;
			mY=*(posM+Mrows);
            
			 
            dist=sqrt((nY-mY)*(nY-mY)+(nX-mX)*(nX-mX));           

			
			if (dist<=threshold)
			{
				if (dist==0) 
					dist=eps;

				
				if(k >= nzmax) 
				{	

					
					int oldnzmax = nzmax;

					
					percent_sparse += percent_step;
					nzmax=(mwSize)ceil((double)Mrows*(double)Nrows*percent_sparse);

					
					if(nzmax <= oldnzmax)
						nzmax=oldnzmax+1;

					
					mxSetNzmax(plhs[0], nzmax);

					
					mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
					mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(mwIndex)));

				
					sr = mxGetPr(plhs[0]);
					irs = mxGetIr(plhs[0]);

				}
				
				
				sr[k]=dist;
				irs[k]=j;
				k++;
			}

        }	
	}
	jcs[Nrows]=k;

	
	

	if (k>0)
		nzmax=k;
	else
		nzmax=1;

	
	mxSetNzmax(plhs[0], nzmax);

	
	mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
	mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(mwIndex)));

	
	sr=mxGetPr(plhs[0]);
	irs=mxGetIr(plhs[0]);
	jcs=mxGetJc(plhs[0]);
}


void calcDistMatrix3D(mxArray *plhs[], double *M, double *N, mwSize Mrows, mwSize Nrows, mwSize nzmax, 
					  double threshold, double eps)
{
    double	mX=0, mY=0, mZ=0;
    double	nX=0, nY=0, nZ=0;
    double  *posM, *posN;
    mwIndex  i,j,k;
    double dist;

	double *sr;
	mwIndex *irs, *jcs;


	sr=mxGetPr(plhs[0]);
	irs=mxGetIr(plhs[0]);
	jcs=mxGetJc(plhs[0]);


	
    k=0;
	for (i=0;i<Nrows;i++) 
	{    
       
        posN=N+i;
        nX=*posN;
        nY=*(posN+Nrows);
		nZ=*(posN+2*Nrows);

	
		jcs[i]=k;
        
        for (j=0;j<Mrows;j++) 
		{       
          
            posM=M+j;
            mX=*posM;
			mY=*(posM+Mrows);
			mZ=*(posM+2*Mrows);
            
            
            dist=sqrt((nY-mY)*(nY-mY)+(nX-mX)*(nX-mX)+(nZ-mZ)*(nZ-mZ));

			
			if (dist<=threshold)
			{
				if (dist==0) 
					dist=eps;

				
				if(k >= nzmax)
				{
					
					int oldnzmax = nzmax;

					
					percent_sparse += percent_step;
					nzmax=(mwSize)ceil((double)Mrows*(double)Nrows*percent_sparse);

					
					if(nzmax <= oldnzmax)
						nzmax=oldnzmax+1;

					
					mxSetNzmax(plhs[0], nzmax);

					
					mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
					mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(mwIndex)));

					
					sr = mxGetPr(plhs[0]);
					irs = mxGetIr(plhs[0]);

				}
				
				
				sr[k]=dist;
				irs[k]=j;
				k++;
			}

        }
	}
	jcs[Nrows]=k;

	

	if (k>0)
		nzmax=k;
	else
		nzmax=1;

	
	mxSetNzmax(plhs[0], nzmax);

	
	mxSetPr(plhs[0], mxRealloc(sr, nzmax*sizeof(double)));
	mxSetIr(plhs[0], mxRealloc(irs, nzmax*sizeof(mwIndex)));

	
	sr=mxGetPr(plhs[0]);
	irs=mxGetIr(plhs[0]);
	jcs=mxGetJc(plhs[0]);
}


void mexFunction(int nlhs, mxArray *plhs[],int nrhs, const mxArray *prhs[])
{
    
	mwSize nzmax; 
	int cmplx; 
	double threshold;
	double eps;

    
    double *M, *N;
	

	mwSize Mrows,Mcols;
	mwSize Nrows,Ncols;

    
	if((nrhs<3) || (nrhs>4))
		mexErrMsgTxt("Three input parameters required.");
	if(nlhs > 1)
		mexErrMsgTxt("One output parameter required.");	


	if(!(mxIsDouble(prhs[0]))) 
		mexErrMsgTxt("Input argument must be of type double.");
	

	switch (nrhs)
	{
	case 3: 
		{
			eps=1e-10;	
			break;
		}
	default:
		{
			eps=mxGetScalar(prhs[3]);
			break;
		}
	
	}

	
	Mrows=mxGetM(prhs[0]);
	Mcols=mxGetN(prhs[0]);
	Nrows=mxGetM(prhs[1]);
	Ncols=mxGetN(prhs[1]);
	

	threshold=mxGetScalar(prhs[2]);

	
	if ((Mcols>3) || (Ncols>3))
		mexErrMsgTxt("Point coordinates in more than 3 dimensions are not supported.");
	
	if (Mcols!=Ncols)
		mexErrMsgTxt("The points in the coordinate matrices have different number of dimensions.");

    
	M=mxGetPr(prhs[0]);
	N=mxGetPr(prhs[1]);	
	
	
    nzmax=(int)ceil((double)Mrows*(double)Nrows*percent_sparse);
	if (nzmax==0)
		nzmax=1;

	
	cmplx=0;

	
	plhs[0]=mxCreateSparse(Mrows,Nrows,nzmax,cmplx);

   	
	if (Mcols==1) { calcDistMatrix1D(plhs,M,N,Mrows,Nrows,nzmax,threshold,eps); }
	if (Mcols==2) { calcDistMatrix2D(plhs,M,N,Mrows,Nrows,nzmax,threshold,eps); }
	if (Mcols==3) { calcDistMatrix3D(plhs,M,N,Mrows,Nrows,nzmax,threshold,eps); }

}
