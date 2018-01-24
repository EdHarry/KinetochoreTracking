#include <stdlib.h>
#include "carmaUtils.h"

/*
 * File: paramsFromVector.c
 * ------------------------
 * The function paramsFromVector extracts CARMA parameters
 * from vector used by the minimizer.
 */



paramsFromVector(double **TOPO, int nNodes, int arOrderMax, 
		 double **maPARAMS, int maOrderMax, 
		 int *topoBIN, int *maBIN, 
		 double *vec, int numParams)
{
 int i, j, k;

 int start, ctr, nLags;

 nLags = arOrderMax + 1;
 ctr = numParams - 1;

 /* Retrieve the MA parameters */

 for (i = nNodes - 1; i >= 0; i--){
   for (j = maOrderMax - 1; j >= 0; j--){
     if (*(maBIN + j + i * maOrderMax)){
	 *(*maPARAMS + j + i * maOrderMax) = *(vec + ctr);
	 ctr--;
     } else {
	 *(*maPARAMS + j + i * maOrderMax) = 0;
     }
   }
 }

 /* Retrieve the AR and X params */ 

 for (i = nNodes - 1; i >= 0; i--){
   for (j = nNodes - 1; j >= 0; j--){

     /* Get the params for node / connection, excluding NaNs */

     if (i == j){
       start = 1;	
       *(*TOPO + j * nLags + i * nNodes * nLags) = 0;
     } else {
       start = 0;
     }

     for (k = arOrderMax; k >= start; k--){
       if (*(topoBIN + k + j * nLags + i * nNodes * nLags)){
	 *(*TOPO + k + j * nLags + i * nNodes * nLags) = *(vec + ctr);
	 ctr--;
       } else{
	 *(*TOPO + k + j * nLags + i * nNodes * nLags) = 0;
       }
     }
   }
 }
}
