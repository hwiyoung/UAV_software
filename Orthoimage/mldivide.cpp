#include "rt_nonfinite.h"
#include "ortho_function.h"
#include "mldivide.h"

void mldivide(const real_T A[9], const real_T B[3], real_T Y[3])
{
	real_T b_A[9];
	int32_T r1;
	int32_T r2;
	int32_T r3;
	real_T maxval;
	real_T a21;
	int32_T rtemp;
	memcpy((void *)&b_A[0], (void *)&A[0], 9U * sizeof(real_T));
	r1 = 0;
	r2 = 1;
	r3 = 2;
	maxval = fabs(A[0]);
	a21 = fabs(A[1]);

	if (a21 > maxval) 
	{
		maxval = a21;
		r1 = 1;
		r2 = 0;
	}

	if (fabs(A[2]) > maxval) 
	{
		r1 = 2;
		r2 = 1;
		r3 = 0;
	}

	b_A[r2] = A[r2] / A[r1];
	b_A[r3] /= b_A[r1];
	b_A[3 + r2] -= b_A[r2] * b_A[3 + r1];
	b_A[3 + r3] -= b_A[r3] * b_A[3 + r1];
	b_A[6 + r2] -= b_A[r2] * b_A[6 + r1];
	b_A[6 + r3] -= b_A[r3] * b_A[6 + r1];
	if (fabs(b_A[3 + r3]) > fabs(b_A[3 + r2]))
	{
		rtemp = r2;
		r2 = r3;
		r3 = rtemp;
	}

	b_A[3 + r3] /= b_A[3 + r2];
	b_A[6 + r3] -= b_A[3 + r3] * b_A[6 + r2];
	Y[0] = B[r1];
	Y[1] = B[r2] - Y[0] * b_A[r2];
	Y[2] = (B[r3] - Y[0] * b_A[r3]) - Y[1] * b_A[3 + r3];
	Y[2] /= b_A[6 + r3];
	Y[0] -= Y[2] * b_A[6 + r1];
	Y[1] -= Y[2] * b_A[6 + r2];
	Y[1] /= b_A[3 + r2];
	Y[0] -= Y[1] * b_A[3 + r1];
	Y[0] /= b_A[r1];
}

/* End of code generation (mldivide.cpp) */
