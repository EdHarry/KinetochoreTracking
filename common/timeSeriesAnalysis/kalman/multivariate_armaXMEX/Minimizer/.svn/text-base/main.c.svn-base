/*
 * File: main.c
 * ------------
 */

#include <stdio.h>
#include <math.h>

#include "amoeba.h"
#include "nrutil.h"
#include "prob.h"

#define NRANSI
#define MP 4
#define NP 3
#define FTOL 1.0e-6


/*
 * Function: fn
 * ------------
 * The type space of objective functions is defined in amoeba.h.  
 */

double fn(double x[], void *prob)
{

  return (pow((x[1] - 16), 2) - pow((x[2] - 16), 2) + pow((x[3] - 16), 2));
}



/*
 * Function: main
 * --------------
 */

int main(void)
{
	int i,nfunc,j,ndim=NP;
	double *x,*y,**p;
	probADT prob;

	x=vector(1,NP);
	y=vector(1,MP);
	p=matrix(1,MP,1,NP);

	prob = NewProb(1, 1);

	for (i=1;i<=MP;i++) {
		for (j=1;j<=NP;j++)
			x[j]=p[i][j]=(i == (j+1) ? 2.0 : 0.0);
		y[i]=fn(x, prob);
	}

	if (amoeba(p, y, ndim, FTOL, fn, &nfunc, prob)) {
	  printf("\nSuccess.\n");
	} else {
	  printf("\nNMAX is exceeded.\n\n");
	  return 0;
	};


	printf("\nNumber of function evaluations: %3d\n",nfunc);
	printf("Vertices of final 3-d simplex and\n");
	printf("function values at the vertices:\n\n");
	printf("%3s %10s %12s %12s %14s\n\n",
		"i","x[i]","y[i]","z[i]","function");

	for (i=1;i<=MP;i++) {
		printf("%3d ",i);
		for (j=1;j<=NP;j++) printf("%12.6f ",p[i][j]);
		printf("%12.6f\n",y[i]);
	}

	printf("\nTrue minimum is at ()\n");

	FreeProb(prob);

	free_matrix(p,1,MP,1,NP);
	free_vector(y,1,MP);
	free_vector(x,1,NP);
	return 0;
}

#undef NRANSI
