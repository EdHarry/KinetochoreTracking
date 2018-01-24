#include <stdio.h>
#include "matOps.h"

/* Calculates A x B = C, where A, B and C are matrices. */

void matrixMultiply(double *pMatA, int mA, int nA,
                    double *pMatB, int mB, int nB,
                    double *pMatC)
{
double tmp = 0.0;
int m, n, p;
int mC = mA;
int nC = nB;
if (nA != mB){ printf("matrix dimension mismatch!\n");}

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
