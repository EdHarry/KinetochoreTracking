#include <math.h>
#include "carmaUtils.h"

/* Calculates partial moving average parameters given input moving average */
/* 	parameters for one ARMA process */

inverseLevinsonDurbinExpoMA(double *maParam, int maOrder, double *maParamP)
{

int i,j,k;
double tempVec[maOrder];
double tempVec2[maOrder];
double multFact;


for (k = 0; k < maOrder; k++){
	tempVec[k] = *(maParam+k);
}

for (i = (maOrder-1); i >= 0; i--){

	*(maParamP+i) = tempVec[i];
	multFact = 1 / ( pow( *(maParamP+i) , 2) - 1 );	     


	for (j = 0; j < i; j++){
		tempVec2[j] = multFact * ( *(maParamP+i) * tempVec[i-j-1] - tempVec[j]);
	}

	for (k = 0; k < i; k++){
		tempVec[k] = tempVec2[k];
	}
}

for (i = 0; i < maOrder; i++){
	*(maParamP+i) = log( (1 - *(maParamP+i) ) / (1 + *(maParamP+i) ) );
}


}
