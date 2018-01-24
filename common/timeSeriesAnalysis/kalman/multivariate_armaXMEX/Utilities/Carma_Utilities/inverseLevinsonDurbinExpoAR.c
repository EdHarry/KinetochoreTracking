#include <math.h>
#include "carmaUtils.h"

/* Calculates partial autoregressive parameters given input autoregressive parameters */
/* 	for one ARMA process */

inverseLevinsonDurbinExpoAR(double *arParam, int arOrder, double *arParamP)
{

int i,j,k;
double tempVec[arOrder];
double tempVec2[arOrder];
double multFact;


for (k = 0; k < arOrder; k++){
	tempVec[k] = *(arParam+k);
}


for (i = (arOrder-1); i >= 0; i--){

	*(arParamP+i) = tempVec[i];
	multFact = 1 / ( 1 - pow( *(arParamP+i) , 2.00) );

	for (j = 0; j < i; j++){
	tempVec2[j] = multFact * ( *(arParamP+i) * tempVec[i-j-1] + tempVec[j]);
	}

	for (k = 0; k < i; k++){
		tempVec[k] = tempVec2[k];
	}
}

for (i = 0; i < arOrder; i++){
	*(arParamP+i) = log( (1 - *(arParamP+i) ) / (1 + *(arParamP+i) ) );
}


}
