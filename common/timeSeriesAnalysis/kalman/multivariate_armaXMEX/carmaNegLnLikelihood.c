/*
 * File: carmaNegLnLikelihood.c
 * ----------------------------
 * The objective function carmaNegLnLikelihood has been modified to
 * accommodate multiple nodes, each contains multiple time series.
 *
 *   
 * Hunter Elliott
 * Modified by Pei-hsin Hsu
 */

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "prob.h"
#include "prob.c"
#include "carmaFitModel.h"

#ifndef MAX
#define	MAX(A, B)	((A) > (B) ? (A) : (B))
#endif



double carmaNegLnLikelihood(double *paramNR, void *d)
{

  probADT prob = (probADT) d;

  /* Convert Numerical Recipes array to C array */
  double *paramV = paramNR + 1; 

  double *TOPO, *TOPOp;

  double *maPARAMS, *maPARAMSp;

  double *innovations, *innovationVars, *wnV;

  double likelihood, sum1, sum2;

  int i, j, k, n;

  int trajLen, nMissing, totLen, totMissing;
  
  double *traj;

  int nMovies = prob->nMovies;
  int nNodes = prob->nNodes;
  int nParams = prob->nParams;
  int arOrderMax = prob->arOrderMax;
  int maOrderMax = prob->maOrderMax;
  int nLags = arOrderMax + 1;
  int *topoBIN = prob->topoBIN;
  int *maBIN = prob->maBIN;

  double wnVariance = prob->wnVariance;



  /* Allocate memory for full and partial parameters */

  TOPO = (double *) malloc(sizeof(double) * nNodes * nNodes * nLags);
  TOPOp = (double *) malloc(sizeof(double) * nNodes * nNodes * nLags);
  maPARAMS = (double *) malloc(sizeof(double) * nNodes * maOrderMax);
  maPARAMSp = (double *) malloc(sizeof(double) * nNodes * maOrderMax);


  /* Retrieve partial parameters from parameter vector */

  paramsFromVector(&TOPOp, nNodes, arOrderMax,
		   &maPARAMSp, maOrderMax,
		   topoBIN, maBIN,
		   paramV, nParams);



  /*
   * Convert partial TOPO parameters to full parameters
   * --------------------------------------------------
   * (1) Along the diagonal are AR parameters.
   * (2) At lag 0 there are no AR parameters.
   * (3) Partial AR parameters are converted to full for innovations calculation.
   * (4) X-parameters are directly copied.
   */

  for (i = 0; i < nNodes; i++){
    for (j = 0; j < nNodes; j++){

      if (i == j) {      /* (1) */

	*(TOPO + i * nLags + j * nNodes * nLags) = 0;      /* (2) */

	levinsonDurbinExpoAR((TOPOp + i * nLags + j * nNodes * nLags + 1),
			     arOrderMax,
			     (TOPO + i * nLags + j * nNodes * nLags + 1));   /* (3) */
      } else {      /* (4) */

	for (k = 0; k < nLags; k++)
	  *(TOPO + i * nLags + j * nNodes * nLags + k) = 
	    *(TOPOp + i * nLags + j * nNodes * nLags + k);

      }
	
    }
  }



  /* Convert partial MA parameters to full parameters */

  for (i = 0; i < nNodes; i++)
    levinsonDurbinExpoMA(maPARAMSp + i * maOrderMax,
			 maOrderMax,
			 maPARAMS + i * maOrderMax);






  /*
   * Calculate -2 ln likelihood
   * --------------------------
   * The actual -2 ln likelihood is declared as likelihood here.
   * For each node comprising several time series (multiple movies),
   * the calculation is based on combination of these references.
   * 
   * Jaqaman K et al. (2006) Biophys. J. 91: 2312-2325. Eq (4)
   * Jones RH (1980) Technometrics 22: 389-395. Eq (3.13) to (3.15)
   */

  likelihood = 0;

  for (i = 0; i < nNodes; i++) {

    /* 
     * For a concatenated movies in a given node,
     * totLen: total length
     * totMissing: total number of missing observations
     */

    totLen = totMissing = 0;

    sum1 = sum2 = 0;      /* Two terms in Eq 3.15 of Jones paper */

    for (n = 0; n < nMovies; n++) {

      traj = prob->data[n].traj;
      trajLen = prob->data[n].trajLen;
      nMissing = 0;

      innovations = (double *) malloc(sizeof(double) * trajLen);
      innovationVars = (double *) malloc(sizeof(double) * trajLen);
      wnV = (double *) malloc(sizeof(double) * trajLen);

      carmaCalcKalmanInnov(traj, trajLen,
			   TOPO, topoBIN,
			   arOrderMax, nNodes,
			   maPARAMS, maOrderMax,
			   wnVariance, i,
			   &innovations, &innovationVars, &wnV,
			   &sum1, &sum2, &nMissing);
      
      totLen += trajLen;
      totMissing += nMissing;
      prob->data[n].nMissing = nMissing;
      /* mexPrintf("c: [%d] trajLen = %d, nMissing = %d\n", n, trajLen, nMissing); */

      free(innovations);
      free(innovationVars);
      free(wnV);
    }

    likelihood += sum1 + (totLen - totMissing) * log(sum2);

  }

  /* mexPrintf("c: sum1 = %g, sum2 = %g\n", sum1, sum2); */

  free(TOPO);
  free(TOPOp);
  free(maPARAMS);
  free(maPARAMSp);

  return likelihood;
}
