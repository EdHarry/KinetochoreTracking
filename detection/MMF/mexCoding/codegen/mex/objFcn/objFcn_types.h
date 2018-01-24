/*
 * objFcn_types.h
 *
 * Code generation for function 'objFcn'
 *
 * C source code generated on: Fri May 25 21:48:50 2012
 *
 */

#ifndef __OBJFCN_TYPES_H__
#define __OBJFCN_TYPES_H__

/* Type Definitions */
#ifndef struct_emxArray__common
#define struct_emxArray__common
typedef struct emxArray__common
{
    void *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
} emxArray__common;
#endif
#ifndef struct_emxArray_real_T
#define struct_emxArray_real_T
typedef struct emxArray_real_T
{
    real_T *data;
    int32_T *size;
    int32_T allocatedSize;
    int32_T numDimensions;
    boolean_T canFreeData;
} emxArray_real_T;
#endif
typedef struct
{
    emxArray_real_T *index;
    real_T psfSigma[2];
} struct_T;

#endif
/* End of code generation (objFcn_types.h) */
