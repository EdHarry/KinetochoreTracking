/*
 * File: covKalmanInit.h
 * ---------------------
 * This interface exports one function, covKalmanInit, which
 * calculates the initial state covariance matrix used in
 * Kalman recursion. This function is a C version of the 
 * covKalmanInit.m implemented by Khuloud Jaqaman in 2004.
 * Reference: Jones, R.H. (1980) Technometrics 22(3): 389-395. 
 * Implementation: Pei-hsin Hsu, March, 2008
 */


/* 
 * Cautionary note:
 * ----------------
 * The major difference between the C implementation and MATLAT is that
 * the initial state covariance matrix has to be passed by reference in C.
 */


#ifndef _covKalmanInit_h
#define _covKalmanInit_h


/* 
 * Function: covKalmanInit
 * -----------------------
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
 */

void covKalmanInit(double *arPar, double *maPar, double *procErrCov, int arOrd, int maOrd, int maxOrd, double **initCov);


#endif
