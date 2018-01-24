/*
 * File: leftdiv.h
 * ---------------
 * This interface exports one procedure, gaussj, which solves linear
 * equations by Gauss-Jordan elimination, and therefore is a C version
 * of the left division operator (A\B) defined in MATLAB.
 *
 * Original code:
 * Press, W.H., et al. Numerical Recipes: the Art of Scientific Computing.
 * 3rd Ed. Cambridge University Press. New York, NY, USA. (2007)
 *
 * Modification:
 * Pei-hsin Hsu, March 2008
 */


#ifndef _leftdiv_h
#define _leftdiv_h


#include <math.h>


/*
 * Function: gaussj
 * ----------------
 * This procedure solves linear equations AX = B by Gauss-Jordan elimination.
 *
 * INPUT
 *     A    : Left-hand side square matrix (n x n)
 *     B    : Right-hand side matrix (n x m) containing m right-hand side vectors.
 *            If only one right-hand side vector exists, B has to be defined as
 *            an n x 1 matrix (n x 1 2D arrays), instead of a 1D array.
 *
 * OUTPUT (pass by reference)
 *     A    : Inverse matrix of the input A 
 *     B    : Solution matrix X whose columns are solutions vectors
 *            If input B comprises only one right-hand side vector,
 *            the solution matrix X is an n x 1 2D arrays.
 *
 * Note: The type of entries of matrices A and B is defined as double precision
 *       floating-point numbers (double), whereas the original gaussj.c from
 *       Numerical Recipes uses float.
 */

void gaussj(double **A, int n, double **B, int m);


#endif
