/*******************************************************************************
 *	Header File For utility subroutines for ARMAX / CARMA model fitting    *
 *									       *
 *******************************************************************************/

/***  Public Variable Descriptions  ***/

/***  Public Subroutine Descriptions  ***/



/* Public Variable External Declarations */

/* Public Subroutine Prototypes */

extern int levinsonDurbinExpoAR(double *arParamP, int arOrder, double *arParam);

extern int levinsonDurbinExpoMA(double *maParamP, int maOrder, double *maParam);

extern int inverseLevinsonDurbinExpoAR(double *arParam, int arOrder, double *arParamP);

extern int inverseLevinsonDurbinExpoMA(double *maParam, int maOrder, double *maParamP);

extern int vectorFromParams(double *TOPO, double *maPARAMS, void *problem, double **vec);

extern int paramsFromVector(double **TOPO, int nNodes, int arOrderMax, double **maPARAMS, int maOrderMax, 
			    int *topoBIN, int *maBIN, double *vec, int numParams);

extern int autoCorr(double *TRAJ, int trajLength, int iNode, int maxLag, double **aCorrs);

