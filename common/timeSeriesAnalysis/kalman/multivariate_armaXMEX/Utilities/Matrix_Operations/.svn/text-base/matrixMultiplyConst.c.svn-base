#include <stdio.h>
#include "matOps.h"

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
