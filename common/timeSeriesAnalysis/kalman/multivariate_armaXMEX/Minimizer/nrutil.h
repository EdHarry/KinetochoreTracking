/*
 * File: nrutil.h
 * --------------
 */

#ifndef _nrutil_h
#define _nrutil_h


/* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
 * utility file nrutil.c.  Do not confuse this file with the same-named
 * file nrutil.c that is supplied in the same subdirectory or archive
 * as the header file nrutil.h.  *That* file contains both ANSI and
 * traditional K&R versions, along with #ifdef macros to select the
 * correct version.  *This* file contains only ANSI C.
 */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>


void nrerror(char error_text[]);
/* Numerical Recipes standard error handler */



double *vector(long nl, long nh);
/* Allocate a double vector with subscript range v[nl..nh] */
/* Modified */



int *ivector(long nl, long nh);
/* allocate an int vector with subscript range v[nl..nh] */


unsigned char *cvector(long nl, long nh);
/* allocate an unsigned char vector with subscript range v[nl..nh] */


unsigned long *lvector(long nl, long nh);
/* allocate an unsigned long vector with subscript range v[nl..nh] */


double *dvector(long nl, long nh);
/* allocate a double vector with subscript range v[nl..nh] */



double **matrix(long nrl, long nrh, long ncl, long nch);
/* Allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
/* Modified */



double **dmatrix(long nrl, long nrh, long ncl, long nch);
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */


int **imatrix(long nrl, long nrh, long ncl, long nch);
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */



double **submatrix(double **a, long oldrl, long oldrh, long oldcl, long oldch,
		  long newrl, long newcl);
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
/* Modified */



double **convert_matrix(double *a, long nrl, long nrh, long ncl, long nch);
/* Allocate a double matrix m[nrl..nrh][ncl..nch] that points to the matrix
 * declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
 * and ncol=nch-ncl+1. The routine should be called with the address
 * &a[0][0] as the first argument.
 * Modified.
 */



double ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh);
/* Allocate a double 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
/* Modified */



void free_vector(double *v, long nl, long nh);
/* Free a double vector allocated with vector() */
/* Modified */


void free_ivector(int *v, long nl, long nh);
/* free an int vector allocated with ivector() */


void free_cvector(unsigned char *v, long nl, long nh);
/* free an unsigned char vector allocated with cvector() */


void free_lvector(unsigned long *v, long nl, long nh);
/* free an unsigned long vector allocated with lvector() */

void free_dvector(double *v, long nl, long nh);
/* free a double vector allocated with dvector() */



void free_matrix(double **m, long nrl, long nrh, long ncl, long nch);
/* Free a double matrix allocated by matrix() */
/* Modified */



void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
/* free a double matrix allocated by dmatrix() */


void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch);
/* free an int matrix allocated by imatrix() */



void free_submatrix(double **b, long nrl, long nrh, long ncl, long nch);
/* Free a submatrix allocated by submatrix() */
/* Modified */



void free_convert_matrix(double **b, long nrl, long nrh, long ncl, long nch);
/* Free a matrix allocated by convert_matrix() */
/* Modified */



void free_f3tensor(double ***t, long nrl, long nrh, long ncl, long nch,
		   long ndl, long ndh);
/* Free a double f3tensor allocated by f3tensor() */
/* Modified */


#endif
