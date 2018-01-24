/*
 * objFcn_emxutil.c
 *
 * Code generation for function 'objFcn_emxutil'
 *
 * C source code generated on: Fri May 25 21:48:51 2012
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "objFcn.h"
#include "objFcn_emxutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */

void b_emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush)
{
    emxArray_real_T *emxArray;
    int32_T loop_ub;
    int32_T i;
    *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
    if (doPush) {
        emlrtPushHeapReferenceStack((void *)pEmxArray, (void (*)(void *, boolean_T))emxFree_real_T);
    }
    emxArray = *pEmxArray;
    emxArray->data = (real_T *)NULL;
    emxArray->numDimensions = numDimensions;
    emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
    emxArray->allocatedSize = 0;
    emxArray->canFreeData = TRUE;
    loop_ub = numDimensions - 1;
    for (i = 0; i <= loop_ub; i++) {
        emxArray->size[i] = 0;
    }
}

void c_emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush)
{
    emxArray_real_T *emxArray;
    int32_T loop_ub;
    int32_T i;
    *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
    if (doPush) {
        emlrtPushHeapReferenceStack((void *)pEmxArray, (void (*)(void *, boolean_T))emxFree_real_T);
    }
    emxArray = *pEmxArray;
    emxArray->data = (real_T *)NULL;
    emxArray->numDimensions = numDimensions;
    emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
    emxArray->allocatedSize = 0;
    emxArray->canFreeData = TRUE;
    loop_ub = numDimensions - 1;
    for (i = 0; i <= loop_ub; i++) {
        emxArray->size[i] = 0;
    }
}

void emxEnsureCapacity(emxArray__common *emxArray, int32_T oldNumel, int32_T elementSize)
{
    int32_T newNumel;
    int32_T loop_ub;
    int32_T i;
    void *newData;
    newNumel = 1;
    loop_ub = emxArray->numDimensions - 1;
    for (i = 0; i <= loop_ub; i++) {
        newNumel *= emxArray->size[i];
    }
    if (newNumel > emxArray->allocatedSize) {
        loop_ub = emxArray->allocatedSize;
        if (loop_ub < 16) {
            loop_ub = 16;
        }
        while (loop_ub < newNumel) {
            loop_ub <<= 1;
        }
        newData = calloc((uint32_T)loop_ub, (uint32_T)elementSize);
        if (emxArray->data != NULL) {
            memcpy(newData, emxArray->data, (uint32_T)(elementSize * oldNumel));
            if (emxArray->canFreeData) {
                free(emxArray->data);
            }
        }
        emxArray->data = newData;
        emxArray->allocatedSize = loop_ub;
        emxArray->canFreeData = TRUE;
    }
}

void emxFreeStruct_struct_T(struct_T *pStruct)
{
    emxFree_real_T(&pStruct->index);
}

void emxFree_real_T(emxArray_real_T **pEmxArray)
{
    if (*pEmxArray != (emxArray_real_T *)NULL) {
        if ((*pEmxArray)->canFreeData) {
            free((void *)(*pEmxArray)->data);
        }
        free((void *)(*pEmxArray)->size);
        free((void *)*pEmxArray);
        *pEmxArray = (emxArray_real_T *)NULL;
    }
}

void emxInitStruct_struct_T(struct_T *pStruct, boolean_T doPush)
{
    b_emxInit_real_T(&pStruct->index, 2, doPush);
}

void emxInit_real_T(emxArray_real_T **pEmxArray, int32_T numDimensions, boolean_T doPush)
{
    emxArray_real_T *emxArray;
    int32_T loop_ub;
    int32_T i;
    *pEmxArray = (emxArray_real_T *)malloc(sizeof(emxArray_real_T));
    if (doPush) {
        emlrtPushHeapReferenceStack((void *)pEmxArray, (void (*)(void *, boolean_T))emxFree_real_T);
    }
    emxArray = *pEmxArray;
    emxArray->data = (real_T *)NULL;
    emxArray->numDimensions = numDimensions;
    emxArray->size = (int32_T *)malloc((uint32_T)(sizeof(int32_T) * numDimensions));
    emxArray->allocatedSize = 0;
    emxArray->canFreeData = TRUE;
    loop_ub = numDimensions - 1;
    for (i = 0; i <= loop_ub; i++) {
        emxArray->size[i] = 0;
    }
}
/* End of code generation (objFcn_emxutil.c) */
