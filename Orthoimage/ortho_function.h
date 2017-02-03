#ifndef __ORTHO_FUNCTION_H__
#define __ORTHO_FUNCTION_H__

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"

#include "rtwtypes.h"

extern void ortho_function(const real_T EO[6], const double em_r[6000*4000], const double
  em_g[6000*4000], const double em_b[6000*4000]);
#endif


void orthoimage_start();