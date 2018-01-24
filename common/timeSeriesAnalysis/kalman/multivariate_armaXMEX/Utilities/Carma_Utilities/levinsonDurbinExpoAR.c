#include <math.h>
#include <stdio.h>
#include "carmaUtils.h"

/* Calculates Autoregressive parameters given input partial autoregressive parameters */
/* 	for one ARMA process */

levinsonDurbinExpoAR(double *arParamP, int arOrder, double *arParam)
{

int i,j,k;
double tempVec[arOrder];
double tempMat[2][arOrder];

for (i = 0; i < arOrder; i++){
	tempVec[i] = (1 - exp( *(arParamP+i) ) ) / (1 + exp( *(arParamP+i) ) );
}

for (i = 0; i < arOrder; i++){

	tempMat[0][i] = 0;
	tempMat[1][i] = tempVec[i];

	for (j = 0; j < i; j++){
	tempMat[1][j] = tempMat[0][j] - ( tempVec[i] * tempMat[0][i-j-1] );
	}

	for (k = 0; k <= i; k++){
	tempMat[0][k] = tempMat[1][k];
	}
}

for (k = 0; k < arOrder; k++){
	*(arParam+k) = tempMat[1][k];
}

}
