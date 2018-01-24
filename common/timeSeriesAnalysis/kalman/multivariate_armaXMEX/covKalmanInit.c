/*
 * File: covKalmanInit.c
 * ---------------------
 * This file implements the covKalmanInit.h interface.
 * Please see the reference for equations used in this file.
 *
 * Note: A few portions (e.g., printing out entries of a matrix)
 *       of the function are suppressed in the final version, but
 *       they will be helpful in future debugging.
 *
 * Reference: Jones, R.H. (1980) Technometrics 22(3): 389-395.
 * Implementation: Pei-hsin Hsu, March 2008
 */


#include <stdio.h>
#include <stdlib.h>
#include "covKalmanInit.h"
#include "./Utilities/Matrix_Operations/leftdiv.h"


/* Private function prototypes */

static double **NewMat(int row, int col);
static void FreeMat(double **mat, int row);
static void PrintMat(double **mat, int row, int col);



/*
 * Function: covKalmanInit
 * -----------------------
 * The implementation of the only one function exported by covKalmanInit.h interface.
 *
 * INPUT
 *     arPar        : Autoregressive coefficients array
 *     maPar        : Moving average coefficients array
 *     procErrCov   : Vector G (2.16)
 *     arOrd        : Proposed AR order
 *     maOrd        : Proposed MA order
 *     maxOrd       : Max(arOrd, maOrd + 1) (2.3)
 *     initCov      : Initial state covariance matrix
 *                    Note: Memory must be allocated BEFORE function call.
 *
 * OUTPUT
 *     initCov      : Initial state covariance matrix, passed by reference
 *
 * LOCAL VARIABLES
 *     nNonZero     : Number of nonzero terms of the rhs vector (4.11)
 *     lenC         : Length of covariance vector C (4.11, and 4.9 if maOrd > arOrd)
 *     rhsVec       : Right hand side vector containing sums of products of beta and gamma terms.  (4.11)
 *     lhsMat       : Left hand side matrix containing covariances and coefficients (alpha terms). (4.11)
 *     maParMod     : Modified moving average vector, in which 1 is added to maPar as the first component.
 *     lhsMatMod    : Modified left hand side matrix. Ignoring that C(-i) = C(i).
 *     sum          : Temporary storage for local summation.
 */

void covKalmanInit(double *arPar, double *maPar, double *procErrCov, int arOrd, int maOrd, int maxOrd, double **initCov)
{
  int i, j, k;
  int nNonZero, nRow, nCol, lenC;
  double **rhsVec, **lhsMat, *C;
  double maParMod[maOrd + 1], lhsMatMod[arOrd + 1][2 * arOrd + 1];
  double sum;

  
  /*
   * Allocate memories for rhsVec and lhsMat. 
   * Note: Vectors have to be defined as one-column matrices
   * to call the function gaussj exported by leftdiv.h
   */
 

  rhsVec = NewMat(arOrd + 1, 1);            
  lhsMat = NewMat(arOrd + 1, arOrd + 1);


  /* Construct modified moving average parameter array */

  maParMod[0] = 1;

  for (i = 1; i < maOrd + 1; i++) maParMod[i] = maPar[i-1];


  /* Construct rhs vector of (4.11) */

  nNonZero = (arOrd < maOrd) ? (arOrd + 1) : (maOrd + 1); /* Number of non zero terms */

  for (i = 0; i < nNonZero; i++) {
    sum = 0;

    for (j = i; j < maOrd + 1; j++) {
      sum += maParMod[j] * procErrCov[j-i];
      /* printf("beta[%d] = %g gamma[%d] = %g\n", j, maParMod[j], j-i, procErrCov[j-i]); */
    }

    rhsVec[i][0] = sum;
    /* printf("rhsVec[%d][0] = %g\n", i, rhsVec[i][0]); */
  }

  for (i = nNonZero; i < arOrd + 1; i++){
    rhsVec[i][0] = 0;
    /* printf("rhsVec[%d][0] = %g\n", i, rhsVec[i][0]); */
  }

  /* printf("\n"); */


  /* Construct modified lhsMat of (4.11) */

  nRow = arOrd + 1;
  nCol = 2 * arOrd + 1;

  for (i = 0; i < nRow; i++) {
    for (j = 0; j < nCol; j ++) lhsMatMod[i][j] = 0;
  }

  for (j = 0; j < arOrd; j++) lhsMatMod[0][j] = -1 * arPar[arOrd - j - 1];
  lhsMatMod[0][arOrd] = 1;

  for (i = 1; i < nRow; i++) {
    for (j = i; j < nCol; j ++) lhsMatMod[i][j]  = lhsMatMod[i - 1][j - 1];
  }


  /* Construct lhsMat of (4.11) */

  nRow = arOrd + 1;
  nCol = arOrd + 1;

  for (i = 0; i < nRow; i++) {
    for (j = 0; j < nCol; j ++) lhsMat[i][j] = 0;
    lhsMat[i][0] = lhsMatMod[i][arOrd];
  }

  for (i = 0; i < nRow; i++) {
    for (j = 1; j < nCol; j++) {
      lhsMat[i][j] = lhsMatMod[i][arOrd - j] + lhsMatMod[i][arOrd + j];
    }
  }

  /*
  for (i = 0; i < nRow; i++) {
    for (j = 0; j < nCol; j ++) printf("lhsMat[%d][%d] = %g\t", i, j, lhsMat[i][j]);
    printf("\n");
  }
  */


  /* Solve linear equations (4.11) */

  gaussj(lhsMat, arOrd + 1, rhsVec, 1);


  /* Retrieve covariances */

  lenC = (arOrd < maOrd) ? (maOrd + 1) : (arOrd + 1); /* Length of covariance vector C */

  C = (double *) malloc(sizeof(double) * lenC);

  
  for (i = 0; i < (arOrd + 1); i++) {
    C[i] = rhsVec[i][0];
    /* printf("C(%d) = %g\n", i, C[i]); */
  }
  


  /* Calculate the remaining covariances, if exist. (4.9) */

  if (arOrd < maOrd) {

    for (i = (arOrd + 1); i < maxOrd; i++) {

      C[i] = 0;

      if (arOrd > 0) {
	for (j = 0; j < arOrd; j++) C[i] += arPar[j] * C[i - j - 1];
      }

      if (maOrd > 0) {
	for (j = (i - 1); j < maOrd; j++) C[i] += maPar[j] * procErrCov[j - i + 1];
      }

      /* printf("C(%d) = %g\n", i, C[i]); */

    }

  }


  /* Construct initStateCov matrix (4.12) */

  for (j = 0; j < maxOrd; j++) initCov[0][j] = C[j];

  for (i = 1; i < maxOrd; i++) {
    for (j = i; j < maxOrd; j++) {
      sum = 0;
      for (k = 0; k < i; k++) sum += procErrCov[k] * procErrCov[k + j - i];
      initCov[i][j]  = C[j - i] - sum;
    }
  }

  for (i = 1; i < maxOrd; i++) {
    for (j = 0; j < i; j++) initCov[i][j] = initCov[j][i];
  }
 

  /* Free memories */

  FreeMat(rhsVec, arOrd + 1);
  FreeMat(lhsMat, arOrd + 1);
  free(C);
} 
  





/* Private function implementation */


/* 
 * Function: NewMat
 * ----------------
 * Allocate memory for a new matrix whose entries are double-presion
 * floating point numbers.
 */

static double **NewMat(int row, int col)
{
  int i;
  double **mat;

  mat = (double **) malloc(sizeof(double *) * row);
  for (i = 0; i < row; i++) mat[i]  = (double *) malloc(sizeof(double) * col);
  return mat;
}


/* 
 * Function: FreeMat
 * -----------------
 * Free the memory associated with each row ptr,
 * as well as the matrix ptr
 */

static void FreeMat(double **mat, int row)
{
  int i;

  for (i = 0; i < row; i++) free(mat[i]);
  free(mat);
}


/*
 * Function: PrintMat
 * ------------------
 * This private function prints out the entries of a given matrix mat,
 * and therefore is helpful in debugging.
 */

static void PrintMat(double **mat, int row, int col)
{
  int i, j;

  printf("\n");
  for (i = 0; i < row; i++){
    for (j = 0; j < col; j++) printf("%-12.6g", mat[i][j]);
    printf("\n");
  }
  printf("\n");
}
