/*
 * File: leftdiv.c
 * ---------------
 * This file combines the gaussj.c and part of nrutil.c (Numerical Recipes utility)
 * to implement the leftdiv.h (left division) interface.
 *
 * Original code:
 * Press, W.H., et al. Numerical Recipes: the Art of Scientific Computing.
 * 3rd Ed. Cambridge University Press. New York, NY (2007)
 *
 * Modification:
 * Pei-hsin Hsu, March 2008
 */

/*
 * Implementation note
 * -------------------
 * Be aware that in Numerical Recipes (NR), indices of matrices and vectors all start
 * from 1, just like MATLAB, despite the fact that the code was written in ANSI C. As
 * a consequence NR has to allocate one additional row and column for a given matrix
 * and never uses the row 0 and column 0. The global variable NR_END is defined as 1 
 * for this off-by-one problem.
 *
 * In this implementation of leftdiv.h interface, the indices have been shifted from original [1..n]
 * to canonical C version [0..n-1], and memories for the additional row and column are not allocated.
 */


#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include "leftdiv.h"


#define NRANSI
#define SWAP(a,b) {temp=(a);(a)=(b);(b)=temp;}
#define NR_END 0 /* Modified, please see the implementation note for details */
#define FREE_ARG char*


/* Private function declarations */

static void nrerror(char error_text[]);
static int *ivector(long nl, long nh);
static void free_ivector(int *v, long nl, long nh);


/*
 * Function: gaussj
 * ----------------
 * Modifications:
 * (1) The type of entries of matrices A and B is defined as double precision
 *     floating-point numbers (double), whereas the original gaussj.c from
 *     Numerical Recipes uses float.
 * (2) The indices have been shifted from original [1..n] to canonical C version [0..n-1].
 *     Please see implementation note for details.
 * (3) gaussj will return 0 if matrix A is singular. (NR error handling will be invoked in original design.)
 */


void gaussj(double **a, int n, double **b, int m)
{
	int *indxc,*indxr,*ipiv;
	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv,temp;

	indxc=ivector(0,n-1);
	indxr=ivector(0,n-1);
	ipiv=ivector(0,n-1);
	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) {
		big=0.0;
		for (j=0;j<n;j++)
			if (ipiv[j] != 1)
				for (k=0;k<n;k++) {
					if (ipiv[k] == 0) {
						if (fabs(a[j][k]) >= big) {
							big=fabs(a[j][k]);
							irow=j;
							icol=k;
						}
					}
				}
		++(ipiv[icol]);
		if (irow != icol) {
		  for (l=0;l<n;l++) SWAP(a[irow][l],a[icol][l]) /* modified */
			for (l=0;l<m;l++) SWAP(b[irow][l],b[icol][l])
		}
		indxr[i]=irow;
		indxc[i]=icol;

		if (a[icol][icol] == 0.0){
		  nrerror("gaussj: Singular Matrix");
		}

		pivinv=1.0/a[icol][icol];
		a[icol][icol]=1.0;
		for (l=0;l<n;l++) a[icol][l] *= pivinv;
		for (l=0;l<m;l++) b[icol][l] *= pivinv;
		for (ll=0;ll<n;ll++)
			if (ll != icol) {
				dum=a[ll][icol];
				a[ll][icol]=0.0;
				for (l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
				for (l=0;l<m;l++) b[ll][l] -= b[icol][l]*dum;
			}
	}
	for (l=n-1;l>=0;l--) {
		if (indxr[l] != indxc[l])
			for (k=0;k<n;k++)
				SWAP(a[k][indxr[l]],a[k][indxc[l]]);
	}
	free_ivector(ipiv,0,n-1);
	free_ivector(indxr,0,n-1);
	free_ivector(indxc,0,n-1);
}






/* Numerical Recipes standard error handler */

static void nrerror(char error_text[])
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}


/* Allocate an integer vector with subscript range v[nl..nh] */

static int *ivector(long nl, long nh)
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}


/* Free an integer vector allocated with ivector() */

static void free_ivector(int *v, long nl, long nh)
{
	free((FREE_ARG) (v+nl-NR_END));
}


#undef SWAP
#undef NRANSI
