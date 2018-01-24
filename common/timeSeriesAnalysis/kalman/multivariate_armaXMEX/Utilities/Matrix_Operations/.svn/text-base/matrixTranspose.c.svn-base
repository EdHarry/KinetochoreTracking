#include <stdio.h>
#include "matOps.h"

void matrixTranspose(double *pInMat, int mInMat, int nInMat, double *pOutMat){
    int m, n;
    for (m = 0; m < mInMat; m++){
        for (n = 0; n < nInMat; n++){
            *(pOutMat + n*mInMat + m) = *(pInMat + m*nInMat + n);
        }
    }

}
