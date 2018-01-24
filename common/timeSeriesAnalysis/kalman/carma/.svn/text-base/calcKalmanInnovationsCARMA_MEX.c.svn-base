
#include "mex.h"
#include "matrix.h"
#include <math.h>

#if !defined(MAX)

#define	MAX(A, B)	((A) > (B) ? (A) : (B))

#endif

void matrixMultiply(double *pMatA, int mA, int nA,
                    double *pMatB, int mB, int nB,
                    double *pMatC, int mC, int nC)
{
    int m, n, p;
    double tmp = 0;
    if (nA != mB){mxErrMsgTxt("matrix dimension mismatch!");}
    
    for (m = 0; m < mC; m++){
        for (n = 0; n < nC; n++){
            for (p = 0; p < nA; p++){
                tmp += *(pMatA + nA*m + p) * *(pMatB + nB*p + n);
            }
            *(pMatC + nC*m + n) = tmp;
            tmp = 0;
        }
    }    
    
}
void matrixAdd(double *pMatA, int mA, int nA,
               double *pMatB, int mB, int nB,
               double *pMatC, int mC, int nC)
{
    int m, n;
    
    if (nA != nB){mxErrMsgTxt("matrix dimension mismatch!");}
    if (mA != mB){mxErrMsgTxt("matrix dimension mismatch!");}
    
    for (m = 0; m < mA; m++){
        for (n = 0; n < nA; n++){    
            *(pMatC + nC*m + n) = *(pMatA + nA*m + n) + *(pMatB + nB*m + n);
        }
    }    
    
}
void matrixSubtract(double *pMatA, int mA, int nA,
                    double *pMatB, int mB, int nB,
                    double *pMatC, int mC, int nC)
{
    int m, n;
    
    if (nA != nB){mxErrMsgTxt("matrix dimension mismatch!");}
    if (mA != mB){mxErrMsgTxt("matrix dimension mismatch!");}
    
    for (m = 0; m < mA; m++){
        for (n = 0; n < nA; n++){    
            *(pMatC + nC*m + n) = *(pMatA + nA*m + n) - *(pMatB + nB*m + n);
        }
    }    
    
}
void matrixMultiplyConst(double *pInMat, int mInMat, int nInMat,
                         double constFactor,
                         double *pOutMat)
{
    int m, n;    
    
    for (m = 0; m < mInMat; m++){
        for (n = 0; n < nInMat; n++){    
            *(pOutMat + nInMat*m + n) = *(pInMat + nInMat*m + n) * constFactor;
        }
    }    
    
}
void matrixTranspose(double *pInMat, int mInMat, int nInMat, double *pOutMat){
    int m, n;
    for (m = 0; m < mInMat; m++){
        for (n = 0; n < nInMat; n++){
            *(pOutMat + n*mInMat + m) = *(pInMat + m*nInMat + n);
        }
    }
    
}
void calcInnovations(
                    double *TRAJ,
                    double *params,
                    double *innovations,
                    double *innovationVars,
                    double *whiteNoise,
                    int connFrom[],
                    int connTo[],  
                    int arOrder,
                    int maOrder,
                    int xOrder,
                    int nExoIn,
                    int currNode,
                    mwSize numParam,
                    mwSize trajLength,
                    mwSize nConn,
                    int nNodes,
                    double *wnVariance
                    )                    
{        
    int j, k, l, m, n;
    int maxOrder = MAX(arOrder,maOrder+1);
    int iConn[nConn];
    int nConnTo = 0;
    double arParamMod[maxOrder], maParamMod[maxOrder];
    double *arPtr, *arParam, *maPtr, *maParam, *G;
    double stateVecT1_T[maxOrder], stateVecT_T[maxOrder];
    double stateCovMatT_T[maxOrder][maxOrder];
    double stateCovMatT1_T[maxOrder][maxOrder];
    double observationVec[maxOrder];
    double tmp = 0;  
    double *matProbe1, *matProbe2;

    /* initialize variables for call to covKalmanInit */
    
    mxArray *plhs1[1], *prhs1[6];
    double *tmp1;
    
    prhs1[0] = mxCreateDoubleMatrix(1,arOrder,mxREAL);
    prhs1[1] = mxCreateDoubleMatrix(1,maOrder,mxREAL);
    prhs1[2] = mxCreateDoubleMatrix(maxOrder,1,mxREAL);
    arParam  = mxGetPr(prhs1[0]);
    maParam  = mxGetPr(prhs1[1]);
    G        = mxGetPr(prhs1[2]);
    prhs1[3] = mxCreateDoubleScalar(arOrder);
    prhs1[4] = mxCreateDoubleScalar(maOrder);
    prhs1[5] = mxCreateDoubleScalar(maxOrder);
    
    /* Initialize state vector, ARMA parameters, G etc */    

    arPtr = (params + (currNode-1) * (arOrder + maOrder));
    maPtr = (params + (currNode-1) * (arOrder + maOrder) + arOrder );
    
    for (k = 0; k < maxOrder; k++){
        
        stateVecT_T[k] = 0;        
        stateVecT1_T[k] = 0;
           
        if (k < arOrder){
            arParam[k] = *(arPtr+k);
            arParamMod[k] = arParam[k];
        }
        else{arParamMod[k] = 0;}
        
        if (k < maOrder){
            maParam[k] = *(maPtr+k);
            maParamMod[k] = maParam[k];
        }
        else{maParamMod[k] = 0;}
        
        /*Initialize vector G*/
        if (k < 1){
            G[k] = 1;
        }
        else {
            tmp = maParamMod[k-1];
            for (j = 1; j <= k; j++){
                tmp += arParamMod[j-1]*G[k-j];
            }
            G[k] = tmp;
        }
        
        /* build "observation vector" */
        
        if (k == 0){
            observationVec[k] = 1;
        }
        else{
            observationVec[k] = 0;
        }

    }
    
    /* determine pertinent connections for current node */
    
    for (k = 0; k < nConn; k++){
        if (connTo[k] == currNode + nExoIn){
            nConnTo++;
            iConn[nConnTo-1] = k;
        }
    }
    
    mexCallMATLAB(1,plhs1,6,prhs1,"covKalmanInit");
    tmp1 = mxGetPr(plhs1[0]);    

    /*multiply initial covariance matrix by white noise variance */    
    matrixMultiplyConst(tmp1,maxOrder,maxOrder,*wnVariance,&stateCovMatT_T[0][0]);
    
    double matGGprime[maxOrder][maxOrder];
    double matGGprimeWN[maxOrder][maxOrder];
    double *GGprimeWN = &matGGprimeWN[0][0];    
    double *GGprime = &matGGprime[0][0];        
    matrixMultiply(G,maxOrder,1,G,1,maxOrder,GGprime,maxOrder,maxOrder);
    matrixMultiplyConst(GGprime,maxOrder,maxOrder,*wnVariance, GGprimeWN);
    
    /* construct transition matrix F */    
    
    double F[maxOrder][maxOrder];
    double Fprime[maxOrder][maxOrder];
    
    for (k = 0; k < maxOrder; k++){
        for (l = 0; l < maxOrder; l++){
            if (l == k + 1){
                F[k][l] = 1.0;
                Fprime[l][k] = 1.0;
            }
            else if (k == maxOrder - 1){
                F[k][l] = arParamMod[maxOrder-(l+1)];
                Fprime[l][k] = arParamMod[maxOrder-(l+1)];
            }
            else {
                F[k][l] = 0;
                Fprime[l][k] = 0;
            }

        }
    }

    /* Construct input coef matrix B */
    
    double B[nConnTo][maxOrder][xOrder+1];
    
    for (m = 0; m < nConnTo ; m++){        
        for (k = 0; k < maxOrder; k++){
            for (l = 0; l < xOrder+1 ; l ++){
                if (k == maxOrder-1){
                    B[m][k][l] = params[ (arOrder+maOrder)*nNodes + iConn[m]*(xOrder+1) - l + xOrder];
                }
                else{
                    B[m][k][l] = 0;
                }
            }
        }
    }

    /* Initialize variables for innovation calculation */
    double innovationVar[trajLength];
    double delta[maxOrder];
    double *inputVectors[nConnTo];
    int t2;
    double tmpMat[maxOrder][maxOrder];
    double tmpMat2[maxOrder][maxOrder];
    double tmpVec[maxOrder];
    bool nanInInput = 0;
    bool nanPresent = 0;
    
    /* Predict and filter timeseries */
    /* This is the meat! */
                
    
    matProbe1 = &stateCovMatT1_T[0][0];
    matProbe2 = &stateCovMatT_T[0][0];
    
    
    for (j = MAX(xOrder,0); j < trajLength; j++){

        /* get input timeseries */
        
        
        t2 = j - 1 + maxOrder - xOrder;
        
        for (m = 0 ; m < nConnTo; m++){
            inputVectors[m] = TRAJ + (trajLength*(connFrom[iConn[m]]-1)) + t2;
        } 
        
        /* check for NaNs in input vectors */        
        for (m = 0; m < nConnTo; m++){                        
            for (k = 0; k < xOrder+1; k++){
                
                if ( mxIsNaN( *(inputVectors[m] + k) ) ){
                    nanInInput = 1;
                }                
            }
        }        
        
        if ( !nanInInput){
            
            /* If no NaN, predict the state at t+1 */

            /* apply transition matrix*/
            matrixMultiply(&F[0][0],maxOrder,maxOrder,&stateVecT_T[0],maxOrder,1,&stateVecT1_T[0],maxOrder,1);

            for (m = 0; m < nConnTo; m++){
                /*multiply input matrix B by input trajectory */
                matrixMultiply(&B[m][0][0],maxOrder,xOrder+1,inputVectors[m],xOrder+1,1,&tmpVec[0],maxOrder,1);
                /* add resulting vector to stateVec */
                for (k = 0; k < maxOrder; k++){
                    stateVecT1_T[k] += tmpVec[k];
                }
            }

            /* Predict the covariance matrix of the state */

            /* calculate F * stateCovT_T */
            matrixMultiply(&F[0][0],maxOrder,maxOrder,&stateCovMatT_T[0][0],maxOrder,maxOrder,&tmpMat[0][0],maxOrder,maxOrder);
            /* above times F' */
            matrixMultiply(&tmpMat[0][0],maxOrder,maxOrder,&Fprime[0][0],maxOrder,maxOrder,&tmpMat2[0][0],maxOrder,maxOrder);
            /* add G*G' */
            matrixAdd(&tmpMat2[0][0],maxOrder,maxOrder,GGprimeWN,maxOrder,maxOrder,&stateCovMatT1_T[0][0],maxOrder,maxOrder);
        }
        /*check if there is an observation at current timepoint */
        
        
        if ( mxIsNaN( *(TRAJ + (trajLength * (currNode+nExoIn-1)) + j) ) | nanInInput ){
            nanPresent = 1;
            nanInInput = 0;
        }
        
        if (nanPresent){
            for (k = 0; k < maxOrder; k++){

                stateVecT_T[k] = stateVecT1_T[k];
                
                for (m = 0; m < maxOrder; m++){
                    stateCovMatT_T[k][m] = stateCovMatT1_T[k][m];                                        
                }                               
                
            }
            
             /* insure that covariance matrix is symmetric */
            matrixTranspose(&stateCovMatT_T[0][0],maxOrder,maxOrder,&tmpMat[0][0]);
            matrixAdd(&stateCovMatT_T[0][0],maxOrder,maxOrder,&tmpMat[0][0],maxOrder,maxOrder,&tmpMat2[0][0],maxOrder,maxOrder);
            matrixMultiplyConst(&tmpMat2[0][0],maxOrder,maxOrder,.5,&stateCovMatT_T[0][0]);                                                            
            
            nanPresent = 0;

        }        
        else{
            /* If all observations are present */

            
            /* Calculate innovations */
            *(innovations + j + ((currNode-1) * trajLength)) = *(TRAJ + (trajLength * (currNode+nExoIn-1)) + j) - stateVecT1_T[0];
            /* And variance in innovations */
            *(innovationVar + j + ((currNode-1) * trajLength)) = stateCovMatT1_T[0][0] + ( *(TRAJ + (nNodes+currNode+2*nExoIn-1)*trajLength + j ) * *(TRAJ + (nNodes+currNode+2*nExoIn-1)*trajLength + j)); 

            /* calculate Delta */
            matrixMultiplyConst(&stateCovMatT1_T[0][0],maxOrder,1,( 1 / *(innovationVar + j + (currNode-1) * trajLength) ),&delta[0]);
            /* delta times the next innovation */
            matrixMultiplyConst(&delta[0],maxOrder,1,*(innovations + j + ((currNode-1) * trajLength)),&tmpVec[0]);        
            /* add this to previous state prediction */
            matrixAdd(&stateVecT1_T[0],maxOrder,1,&tmpVec[0],maxOrder,1,&stateVecT_T[0],maxOrder,1);

            /* update covatiance matrix */
            matrixMultiply(&delta[0],maxOrder,1,&observationVec[0],1,maxOrder,&tmpMat[0][0],maxOrder,maxOrder);
            matrixMultiply(&tmpMat[0][0],maxOrder,maxOrder,&stateCovMatT1_T[0][0],maxOrder,maxOrder,&tmpMat2[0][0],maxOrder,maxOrder);
            matrixSubtract(&stateCovMatT1_T[0][0],maxOrder,maxOrder,&tmpMat2[0][0],maxOrder,maxOrder,&stateCovMatT_T[0][0],maxOrder,maxOrder);

            /* insure that covariance matrix is symmetric */
            matrixTranspose(&stateCovMatT_T[0][0],maxOrder,maxOrder,&tmpMat[0][0]);
            matrixAdd(&stateCovMatT_T[0][0],maxOrder,maxOrder,&tmpMat[0][0],maxOrder,maxOrder,&tmpMat2[0][0],maxOrder,maxOrder);
            matrixMultiplyConst(&tmpMat2[0][0],maxOrder,maxOrder,.5,&stateCovMatT_T[0][0]);                                                                        
            
            /* Calculate white noise. Same as innovation without obs error */
            *(whiteNoise + j + ((currNode-1) * trajLength) ) = stateVecT_T[0] - stateVecT1_T[0];

        }
    }
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{  
    
    double *params, *TRAJ;
    double *innovations, *innovationVar, *whiteNoise;
    double *connFromTMP, *connToTMP;
    double *wnVariance;
    int arOrder, maOrder, xOrder;
    int nExoIn, nNodes;
    mxArray *tmp;
    double *tmp2;    
    
    mwSize numParam, trajLength, sizebuf, nConn;
    
    int j, currNode;    
    
    /*Determine number of input parameters */
    params = mxGetPr(prhs[0]);
    numParam = mxGetM(prhs[0]);        
    
    /* Get pointers for input matrices */
                
    tmp = mxGetField(prhs[1],0,"arOrder");
    tmp2 = mxGetPr(tmp);
    arOrder = (int) *tmp2;	
    
    tmp = mxGetField(prhs[1],0,"maOrder");
    tmp2 = mxGetPr(tmp);
    maOrder = (int) *tmp2;
    
    tmp = mxGetField(prhs[1],0,"xOrder");
    tmp2 = mxGetPr(tmp);
    xOrder = (int) *tmp2;

    tmp = mxGetField(prhs[1],0,"TRAJ");
    trajLength = mxGetM(tmp);
    TRAJ = mxGetPr(tmp);
    
    tmp = mxGetField(prhs[1],0,"connFrom");
    nConn = mxGetM(tmp);
    connFromTMP = mxGetPr(tmp);

    tmp = mxGetField(prhs[1],0,"connTo");
    connToTMP = mxGetPr(tmp);
    
    int connTo[nConn], connFrom[nConn];
    for (j = 0; j < nConn; j++){
        connTo[j] = (int) connToTMP[j];
        connFrom[j] = (int) connFromTMP[j];
    }
    
    tmp = mxGetField(prhs[1],0,"nExoIn");
    tmp2 = mxGetPr(tmp);
    nExoIn = (int) *tmp2;
    
    tmp = mxGetField(prhs[1],0,"nNodes");
    tmp2 = mxGetPr(tmp);
    nNodes = (int) *tmp2;    
    
    tmp = mxGetField(prhs[1],0,"wnVariance");
    wnVariance = mxGetPr(tmp);
    
    /* Get pointers for output matrices */
    
    plhs[0] = mxCreateDoubleMatrix(trajLength,nNodes,mxREAL);
    innovations = mxGetPr(plhs[0]);

    plhs[1] = mxCreateDoubleMatrix(trajLength,nNodes,mxREAL);
    innovationVar = mxGetPr(plhs[1]);
    
    plhs[2] = mxCreateDoubleMatrix(trajLength,nNodes,mxREAL);
    whiteNoise = mxGetPr(plhs[2]);
    
    
    for (currNode = 1; currNode <= nNodes; currNode++){


        calcInnovations(TRAJ,params,innovations,innovationVar,whiteNoise,
                        connFrom,connTo,arOrder,maOrder,xOrder,nExoIn,
                        currNode,numParam,trajLength,nConn,nNodes,wnVariance);


    }    

}
