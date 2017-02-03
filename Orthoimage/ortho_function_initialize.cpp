#include "rt_nonfinite.h"
#include "ortho_function.h"
#include "ortho_function_initialize.h"


void ortho_function_initialize(void)
{
  rt_InitInfAndNaN(8U);
}
