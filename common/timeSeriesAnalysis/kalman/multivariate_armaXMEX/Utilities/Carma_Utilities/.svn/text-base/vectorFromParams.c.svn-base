#include <stdlib.h>
#include "carmaUtils.h"
#include "../../prob.h"
#include "../../prob.c"

/*
 * File: vectorFromParams.c
 * ------------------------
 * The function converts CARMA parameters to a vector for the minimizer,
 * using the connections allowed by adjacency matrix tryCONN.
 */



int vectorFromParams(double *TOPO, double *maPARAMS, void *problem, double **vec)
{

  int i, j, k, ctr, start;
  int nLags, nNodes, arOrderMax, maOrderMax, maxNumParams;

  probADT prob = (probADT) problem;

  nLags = prob->arOrderMax + 1;
  nNodes = prob->nNodes;

  arOrderMax = prob->arOrderMax;
  maOrderMax = prob->maOrderMax;

  maxNumParams = nLags * nNodes * nNodes + nNodes * maOrderMax;

  *vec = (double *) malloc((sizeof(double) * maxNumParams));

  /* Set the first parameter to zero for numerical recipes */

  *vec[0] = 0;

  ctr = 1;

  /* First add AR and X params */ 

  for (i = 0; i < nNodes; i++){
    for (j = 0; j < nNodes; j++){	  
	
	/* if it is a node instead of a connection, skip lag 0 */

      (i == j) ? (start = 1) : (start = 0);
	
      for (k = start; k <= arOrderMax; k++){
	  
	  /* Get the params for node / connection, if a parameter is present as indicated by topoBIN */
  
	if( *(prob->topoBIN + k + j * nLags + i * nNodes * nLags) ){	    
	  *(*vec + ctr) = *(TOPO + k + j * nLags + i * nNodes * nLags);
	  ctr++;
	}
      }	
    }  
  }
  

  
  /* Now add the moving average parameters */

  for (i = 0; i < nNodes; i++){
    
    /* Retrieve the Moving Average (MA) parameters at this node, if maBIN indicates a parameter present */

    for (j = 0; j < maOrderMax; j++){	
      if (*(prob->maBIN + j + i * maOrderMax)){
	*(*vec + ctr) = *(maPARAMS + j + i * maOrderMax);
	ctr++;
      }
    }
  }
  
  
  /* Reallocate parameter vector for actual number of params */

  *vec = (double *) realloc(*vec, sizeof(double) * ctr);
  
  /* Return number of parameters, which is one fewer than ctr due to zero at beginning for NR*/

  prob->nParams = ctr - 1;
}
