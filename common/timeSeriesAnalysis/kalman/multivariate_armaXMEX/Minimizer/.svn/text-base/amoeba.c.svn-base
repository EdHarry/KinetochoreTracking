/* 
 * File: amoeba.c
 * --------------
 *
 */


#include <math.h>
#include "amoeba.h"
#include "nrutil.h"
#include "prob.h"

#define NRANSI
#define TINY 1.0e-10
#define NMAX 10000
#define GET_PSUM \
					for (j=1;j<=ndim;j++) {\
					for (sum=0.0,i=1;i<=mpts;i++) sum += p[i][j];\
					psum[j]=sum;}

#define SWAP(a,b) {swap=(a);(a)=(b);(b)=swap;}



/* Private function prototype */

static double amotry(double **p, double y[], double psum[], int ndim, objFn fn, int ihi, double fac, void *prob);



/*
 * Function: amoeba
 * ----------------
 * p      Simplex
 * y      Function values at simplex vertices
 * ftol   Tolerance of decrease in function value in terminating step
 * fn     Objective function
 * nfunk  Number of function evaluations
 * prob   Optional client data that defines the problem
 */

int amoeba(double **p, double y[], int ndim, double ftol, objFn fn, int *nfunk, void *prob)
{
	int i,ihi,ilo,inhi,j,mpts=ndim+1;

	double rtol,sum,swap,ysave,ytry,*psum;

	psum=vector(1,ndim);
	*nfunk=0;
	GET_PSUM
	while (1) {
		ilo=1;
		ihi = y[1]>y[2] ? (inhi=2,1) : (inhi=1,2);
		for (i=1;i<=mpts;i++) {
			if (y[i] <= y[ilo]) ilo=i;
			if (y[i] > y[ihi]) {
				inhi=ihi;
				ihi=i;
			} else if (y[i] > y[inhi] && i != ihi) inhi=i;
		}
		rtol=2.0*fabs(y[ihi]-y[ilo])/(fabs(y[ihi])+fabs(y[ilo])+TINY);

		if (rtol < ftol) {
			SWAP(y[1],y[ilo])
			for (i=1;i<=ndim;i++) SWAP(p[1][i],p[ilo][i])
			break;
		}

		if (*nfunk >= NMAX) return (1); /* nrerror("NMAX exceeded"); */


		*nfunk += 2;

		ytry=amotry(p,y,psum,ndim,fn,ihi,-1.0, prob);

		if (ytry <= y[ilo]){

		  ytry=amotry(p,y,psum,ndim,fn,ihi,2.0, prob);

		} else if (ytry >= y[inhi]) {
			ysave=y[ihi];
			
			ytry=amotry(p,y,psum,ndim,fn,ihi,0.5, prob);

			if (ytry >= ysave) {
				for (i=1;i<=mpts;i++) {
					if (i != ilo) {
						for (j=1;j<=ndim;j++)
							p[i][j]=psum[j]=0.5*(p[i][j]+p[ilo][j]);

						y[i]=(*fn)(psum, prob);
					}
				}
				*nfunk += ndim;
				GET_PSUM
			}
		} else --(*nfunk);
	}

	free_vector(psum,1,ndim);
	return (0);
}



/*
 * Helper function: amotry
 * -----------------------
 * INPUT
 *     p      Simplex
 *     y      Function values at simplex vertices
 *     fn     Objective function
 *     ihi    Index of the high point
 *     fac    Factor by which simplex is reflected from the high point
 *     prob   Optional client data that defines the problem
 */

static double amotry(double **p, double y[], double psum[], int ndim, objFn fn, int ihi, double fac, void *prob)
{
  int j;
	double fac1,fac2,ytry,*ptry;

	ptry=vector(1,ndim);
	fac1=(1.0-fac)/ndim;
	fac2=fac1-fac;
	for (j=1;j<=ndim;j++) ptry[j]=psum[j]*fac1-p[ihi][j]*fac2;

	ytry = (*fn)(ptry, prob);    /* function value at the trial point */

	if (ytry < y[ihi]) {
		y[ihi]=ytry;
		for (j=1;j<=ndim;j++) {
			psum[j] += ptry[j]-p[ihi][j];
			p[ihi][j]=ptry[j];
		}
	}

	free_vector(ptry,1,ndim);
	return ytry;
}



#undef SWAP
#undef GET_PSUM
#undef NMAX
#undef NRANSI
