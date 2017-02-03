#include "rt_nonfinite.h"
#include "ortho_function.h"
#include "min.h"

real_T b_min(const real_T varargin_1[4])
{
	real_T minval;
	int32_T itmp;
	int32_T ix;
	boolean_T guard1 = FALSE;
	boolean_T searchingForNonNaN;
	int32_T k;
	boolean_T exitg1;
	minval = varargin_1[0];
	itmp = 1;
	ix = 0;
	guard1 = FALSE;
	if (rtIsNaN(varargin_1[0]))
	{
		searchingForNonNaN = TRUE;
		k = 2;
		exitg1 = 0U;
		while ((exitg1 == 0U) && (k < 5)) 
		{
			ix++;
			if (!rtIsNaN(varargin_1[ix]))
			{
				minval = varargin_1[ix];
				itmp = k;
				searchingForNonNaN = FALSE;
				exitg1 = 1U;
			}
			else
			{
				k++;
			}
		}

		if (searchingForNonNaN)
		{
		}
		else
		{
			guard1 = TRUE;
		}
	}
	else 
	{
		guard1 = TRUE;
	}

	if (guard1 == TRUE) 
	{
		while (itmp + 1 < 5)
		{
			if (varargin_1[itmp] < minval) 
			{
				minval = varargin_1[itmp];
			}

			itmp++;
		}
	}

	return minval;
}

