/*
 * round.c
 *
 * Code generation for function 'round'
 *
 * C source code generated on: Tue Nov 19 14:12:19 2013
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "mmfMex.h"
#include "round.h"

/* Function Definitions */
void b_round(real_T *x)
{
  *x = muDoubleScalarRound(*x);
}

/* End of code generation (round.c) */
