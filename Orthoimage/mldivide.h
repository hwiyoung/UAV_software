#ifndef __MLDIVIDE_H__
#define __MLDIVIDE_H__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"

#include "rtwtypes.h"
#include "ortho_function_types.h"

extern void mldivide(const real_T A[9], const real_T B[3], real_T Y[3]);
#endif
