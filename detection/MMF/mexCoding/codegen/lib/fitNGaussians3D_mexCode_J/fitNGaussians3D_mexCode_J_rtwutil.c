/*
 * fitNGaussians3D_mexCode_J_rtwutil.c
 *
 * Code generation for function 'fitNGaussians3D_mexCode_J_rtwutil'
 *
 * C source code generated on: Mon May 28 14:49:14 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "fitNGaussians3D_mexCode_J.h"
#include "fitNGaussians3D_mexCode_J_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
real_T rt_powd_snf(real_T u0, real_T u1)
{
  real_T y;
  real_T d1;
  real_T d2;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d1 = fabs(u0);
    d2 = fabs(u1);
    if (rtIsInf(u1)) {
      if (d1 == 1.0) {
        y = rtNaN;
      } else if (d1 > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d2 == 0.0) {
      y = 1.0;
    } else if (d2 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

real_T rt_roundd_snf(real_T u)
{
  real_T y;
  if (fabs(u) < 4.503599627370496E+15) {
    if (u >= 0.5) {
      y = floor(u + 0.5);
    } else if (u > -0.5) {
      y = u * 0.0;
    } else {
      y = ceil(u - 0.5);
    }
  } else {
    y = u;
  }

  return y;
}

/* End of code generation (fitNGaussians3D_mexCode_J_rtwutil.c) */
