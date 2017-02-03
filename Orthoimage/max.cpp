#include "rt_nonfinite.h"
#include "ortho_function.h"
#include "max.h"

real_T b_max(const real_T varargin_1[4])
{
	real_T maxval;
	int32_T itmp;
	int32_T ix;
	boolean_T guard1 = FALSE;
	boolean_T searchingForNonNaN;
	int32_T k;
	boolean_T exitg1;
	maxval = varargin_1[0];
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
				maxval = varargin_1[ix];
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
			if (varargin_1[itmp] > maxval)
			{
				maxval = varargin_1[itmp];
			}

			itmp++;
		}
	}

  return maxval;
}
