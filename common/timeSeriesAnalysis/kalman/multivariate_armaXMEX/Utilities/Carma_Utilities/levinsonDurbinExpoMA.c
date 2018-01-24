#include <math.h>
#include <stdio.h>
#include "carmaUtils.h"

/* Calculates moving average coeficients given partial moving average coeficients */
/* 	for one ARMA process */

levinsonDurbinExpoMA(double *maParamP, int maOrder, double *maParam)
{

int i,j,k;
double tempVec[maOrder];
double tempMat[2][maOrder];

double *ptMat = &tempMat[0][0];

for (i = 0; i < maOrder; i++){
	tempVec[i] = (1 - exp( *(maParamP+i) ) ) / (1 + exp( *(maParamP+i) ) );
}

for (i = 0; i < maOrder; i++){

	tempMat[0][i] = 0;
	tempMat[1][i] = tempVec[i];

	for (j = 0; j < i; j++){
	tempMat[1][j] = tempMat[0][j] + ( tempVec[i] * tempMat[0][i-j-1] );
	}

	for (k = 0; k <= i; k++){
	tempMat[0][k] = tempMat[1][k];
	}
}

for (k = 0; k < maOrder; k++){
	*(maParam+k) = tempMat[1][k];
}

}
