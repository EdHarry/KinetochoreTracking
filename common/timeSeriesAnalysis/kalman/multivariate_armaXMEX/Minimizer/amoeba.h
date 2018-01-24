/*
 * File: amoeba.h
 * --------------
 */



#ifndef _amoeba_h
#define _amoeba_h

#include "prob.h"



/*
 * Type: objFn
 * -----------
 * (1) The type space of objective functions is defined. 
 * (2) The second argument is an arbitrary data that defines the problem,
 *     so that the TOMLAB counterpart of prob is a statement like:
 *     prob = conAssign(...). 
 * (3) If the second argument is not used in function evaluation, as
 *     originally defined in NR, NULL can be passed.
 */

typedef double (*objFn)(double [], void *prob);



/*
 * Function: amoeba
 * ----------------
 * p      Current simplex
 * y      Function values at simplex vertices
 * ftol   Tolerance of decrease in function value in terminating step
 * fn     Objective function
 * nfunk  Number of function evaluations
 * prob   Generic client data that defines the problem
 */

int amoeba(double **p, double y[], int ndim, double ftol, objFn fn, int *nfunk, void *prob);



#endif




