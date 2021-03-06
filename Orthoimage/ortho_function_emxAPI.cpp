#include "rt_nonfinite.h"
#include "ortho_function.h"
#include "ortho_function_emxAPI.h"
#include "ortho_function_emxutil.h"

emxArray_real_T *emxCreateND_real_T(int32_T numDimensions, int32_T *size)
{
  emxArray_real_T *emx;
  int32_T numEl;
  int32_T loop_ub;
  int32_T i;
  b_emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  loop_ub = numDimensions - 1;

  for (i = 0; i <= loop_ub; i++)
  {
		numEl *= size[i];
		emx->size[i] = size[i];
  }

  emx->data = (real_T *)calloc((uint32_T)numEl, sizeof(real_T));
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  return emx;
}


emxArray_real_T *emxCreateWrapperND_real_T(real_T *data, int32_T numDimensions,int32_T *size)
{
  emxArray_real_T *emx;
  int32_T numEl;
  int32_T loop_ub;
  int32_T i;
  b_emxInit_real_T(&emx, numDimensions);
  numEl = 1;
  loop_ub = numDimensions - 1;

  for (i = 0; i <= loop_ub; i++) 
  {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = numDimensions;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}


emxArray_real_T *emxCreateWrapper_real_T(real_T *data, int32_T rows, int32_T cols)
{
  emxArray_real_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  b_emxInit_real_T(&emx, 2);
  numEl = 1;

  for (i = 0; i < 2; i++) 
  {
	numEl *= size[i];
	emx->size[i] = size[i];
  }

  emx->data = data;
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  emx->canFreeData = FALSE;
  return emx;
}


emxArray_real_T *emxCreate_real_T(int32_T rows, int32_T cols)
{
  emxArray_real_T *emx;
  int32_T size[2];
  int32_T numEl;
  int32_T i;
  size[0] = rows;
  size[1] = cols;
  b_emxInit_real_T(&emx, 2);
  numEl = 1;

  for (i = 0; i < 2; i++) 
  {
    numEl *= size[i];
    emx->size[i] = size[i];
  }

  emx->data = (real_T *)calloc((uint32_T)numEl, sizeof(real_T));
  emx->numDimensions = 2;
  emx->allocatedSize = numEl;
  return emx;
}


void emxDestroyArray_real_T(emxArray_real_T *emxArray)
{
  emxFree_real_T(&emxArray);
}
