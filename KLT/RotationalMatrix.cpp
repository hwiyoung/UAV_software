/********************************************************
 * Project:			Automated AT
 * Last updated:	22 August 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "RotationalMatrix.h"
#include <math.h>

/********************************************************************
 * Description:	Construction
 ********************************************************************/
CRotationalMatrix::CRotationalMatrix()
{
	// No initialization required.
}

/********************************************************************
 * Description:	Compute 3D rotational matrix for all EOs in the 
 *				passed-in vector.
 ********************************************************************/
void CRotationalMatrix::Compute3DRotMatrix(VEC_EO const &vec_EO_m, int nimages, int ROT_ORDER, int ROT_TYPE)
{
	EO		eo;
	Matrix	R(3,3);			// 3D Rotational Matrix
	Matrix	Rx(3,3), Ry(3,3), Rz(3,3);
	int		i;	
	double	x;				// Angle		
	double	cos_x, sin_x;	// cos(Angle), sin(Angle)

	// Construct the indexing vector.
	ConstructIndexVector(vec_EO_m);

	// Assign the number of images into the member variable.
	m_nImages = nimages;

	// Compute the rotational matrix for each image.
	for (i = 0; i < nimages; i++)
	{
		eo = vec_EO_m[i];

		if (ROT_TYPE == 1)
		{
			//		|	1		0		0		|
			// Rx =	|	0	  cos(Om) -sin(Om)	|
			//		|	0	  sin(Om) cos(Om)	|		
			
			x		=	eo.dOmega;	// Omega
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			Rx.Set(0, 0, 1.0);
			Rx.Set(0, 1, 0.0);
			Rx.Set(0, 2, 0.0);
			Rx.Set(1, 0, 0.0);
			Rx.Set(1, 1, cos_x);
			Rx.Set(1, 2, -sin_x);
			Rx.Set(2, 0, 0.0);
			Rx.Set(2, 1, sin_x);
			Rx.Set(2, 2, cos_x);
			
			//		|	cos(Ph)		0		sin(Ph)	|
			// Ry =	|	  0			1		   0	|
			//		|	-sin(Ph)	0		cos(Ph)	|
			
			x		=	eo.dPhi;	// Phi
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			Ry.Set(0, 0, cos_x);
			Ry.Set(0, 1, 0.0);
			Ry.Set(0, 2, sin_x);
			Ry.Set(1, 0, 0.0);
			Ry.Set(1, 1, 1.0);
			Ry.Set(1, 2, 0.0);
			Ry.Set(2, 0, -sin_x);
			Ry.Set(2, 1, 0.0);
			Ry.Set(2, 2, cos_x);

			//		|	cos(Ka)		-sin(Ka)	0	|
			// Rz =	|	sin(Ka)		cos(Ka)		0	|
			//		|	  0			   0		1	|
			
			x		=	eo.dKappa;	// Kappa
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			Rz.Set(0, 0, cos_x);
			Rz.Set(0, 1, -sin_x);
			Rz.Set(0, 2, 0.0);
			Rz.Set(1, 0, sin_x);
			Rz.Set(1, 1, cos_x);
			Rz.Set(1, 2, 0.0);
			Rz.Set(2, 0, 0.0);
			Rz.Set(2, 1, 0.0);
			Rz.Set(2, 2, 1.0);
		}
		else // ROT_TYPE = 2
		{
			//		|	1			0		0		|
			// Rx =	|	0		 cos(Om)	sin(Om)	|
			//		|	0		 -sin(Om)	cos(Om)	|		
			
			x		=	eo.dOmega;	// Omega
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			Rx.Set(0, 0, 1.0);
			Rx.Set(0, 1, 0.0);
			Rx.Set(0, 2, 0.0);
			Rx.Set(1, 0, 0.0);
			Rx.Set(1, 1, cos_x);
			Rx.Set(1, 2, sin_x);
			Rx.Set(2, 0, 0.0);
			Rx.Set(2, 1, -sin_x);
			Rx.Set(2, 2, cos_x);
			
			//		|	cos(Ph)		0		-sin(Ph)|
			// Ry =	|	  0			1		   0	|
			//		|	sin(Ph)		0		cos(Ph)	|
			
			x		=	eo.dPhi;	// Phi
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			Ry.Set(0, 0, cos_x);
			Ry.Set(0, 1, 0.0);
			Ry.Set(0, 2, -sin_x);
			Ry.Set(1, 0, 0.0);
			Ry.Set(1, 1, 1.0);
			Ry.Set(1, 2, 0.0);
			Ry.Set(2, 0, sin_x);
			Ry.Set(2, 1, 0.0);
			Ry.Set(2, 2, cos_x);

			//		|	cos(Ka)		sin(Ka)		0	|
			// Rz =	|	-sin(Ka)	cos(Ka)		0	|
			//		|	  0			   0		1	|
			
			x		=	eo.dKappa;	// Kappa
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			Rz.Set(0, 0, cos_x);
			Rz.Set(0, 1, sin_x);
			Rz.Set(0, 2, 0.0);
			Rz.Set(1, 0, -sin_x);
			Rz.Set(1, 1, cos_x);
			Rz.Set(1, 2, 0.0);
			Rz.Set(2, 0, 0.0);
			Rz.Set(2, 1, 0.0);
			Rz.Set(2, 2, 1.0);		
		}

		// Rotational matrix
		if (ROT_ORDER == 1)
			R = Rx * Ry * Rz;
		else
			R = Rz * Ry * Rx;

		// Store the rotational matrix into the vector.
		m_vec_RotMatrix.push_back(R);
	}
}

/********************************************************************
 * Description:	Compute the derivative of 3D rotational matrix for  
 *				all EOs in the passed-in vector.
 ********************************************************************/
void CRotationalMatrix::ComputeDerivative3DRotMatrix(VEC_EO const &vec_EO_m, int nimages, int ROT_ORDER, int ROT_TYPE)
{
	EO		eo;
	Matrix	dR(3,3);							
	Matrix	Rx(3,3), Ry(3,3), Rz(3,3);
	Matrix	dRx_om(3,3), dRy_ph(3,3), dRz_kp(3,3);
	int		i;
	double	x;				// Angle		
	double	cos_x, sin_x;	// cos(Angle), sin(Angle)

	// Construct the indexing vector.
	ConstructIndexVector(vec_EO_m);

	// Assign the number of images into the member variable.
	m_nImages = nimages;

	// Compute the derivative rotational matrix for each image.
	for (i = 0; i < nimages; i++)
	{
		eo = vec_EO_m[i];

		if (ROT_TYPE == 1)
		{
			// --- Omega ---
			
			x		=	eo.dOmega;	
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			// Rotational matrix 
			Rx.Set(0, 0, 1.0);			Rx.Set(0, 1, 0.0);			Rx.Set(0, 2, 0.0);
			Rx.Set(1, 0, 0.0);			Rx.Set(1, 1, cos_x);		Rx.Set(1, 2, -sin_x);
			Rx.Set(2, 0, 0.0);			Rx.Set(2, 1, sin_x);		Rx.Set(2, 2, cos_x);

			// Derivative rotational matrix 
			dRx_om.Set(0, 0, 0.0);		dRx_om.Set(0, 1, 0.0);		dRx_om.Set(0, 2, 0.0);
			dRx_om.Set(1, 0, 0.0);		dRx_om.Set(1, 1, -sin_x);	dRx_om.Set(1, 2, -cos_x);
			dRx_om.Set(2, 0, 0.0);		dRx_om.Set(2, 1, cos_x);	dRx_om.Set(2, 2, -sin_x);

			// --- Phi ---
			
			x		=	eo.dPhi;
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			// Rotational matrix 
			Ry.Set(0, 0, cos_x);		Ry.Set(0, 1, 0.0);			Ry.Set(0, 2, sin_x);
			Ry.Set(1, 0, 0.0);			Ry.Set(1, 1, 1.0);			Ry.Set(1, 2, 0.0);
			Ry.Set(2, 0, -sin_x);		Ry.Set(2, 1, 0.0);			Ry.Set(2, 2, cos_x);

			// Derivative rotational matrix 
			dRy_ph.Set(0, 0, -sin_x);	dRy_ph.Set(0, 1, 0.0);		dRy_ph.Set(0, 2, cos_x);
			dRy_ph.Set(1, 0, 0.0);		dRy_ph.Set(1, 1, 0.0);		dRy_ph.Set(1, 2, 0.0);
			dRy_ph.Set(2, 0, -cos_x);	dRy_ph.Set(2, 1, 0.0);		dRy_ph.Set(2, 2, -sin_x);

			// --- Kappa ---
			
			x		=	eo.dKappa;
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			// Rotational matrix 
			Rz.Set(0, 0, cos_x);		Rz.Set(0, 1, -sin_x);		Rz.Set(0, 2, 0.0);
			Rz.Set(1, 0, sin_x);		Rz.Set(1, 1, cos_x);		Rz.Set(1, 2, 0.0);
			Rz.Set(2, 0, 0.0);			Rz.Set(2, 1, 0.0);			Rz.Set(2, 2, 1.0);

			// Derivative rotational matrix 
			dRz_kp.Set(0, 0, -sin_x);	dRz_kp.Set(0, 1, -cos_x);	dRz_kp.Set(0, 2, 0.0);
			dRz_kp.Set(1, 0, cos_x);	dRz_kp.Set(1, 1, -sin_x);	dRz_kp.Set(1, 2, 0.0);
			dRz_kp.Set(2, 0, 0.0);		dRz_kp.Set(2, 1, 0.0);		dRz_kp.Set(2, 2, 0.0);
		}
		else // ROT_TYPE = 2
		{
			// --- Omega ---
			
			x		=	eo.dOmega;	
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			// Rotational matrix 
			Rx.Set(0, 0, 1.0);			Rx.Set(0, 1, 0.0);			Rx.Set(0, 2, 0.0);
			Rx.Set(1, 0, 0.0);			Rx.Set(1, 1, cos_x);		Rx.Set(1, 2, sin_x);
			Rx.Set(2, 0, 0.0);			Rx.Set(2, 1, -sin_x);		Rx.Set(2, 2, cos_x);

			// Derivative rotational matrix 
			dRx_om.Set(0, 0, 0.0);		dRx_om.Set(0, 1, 0.0);		dRx_om.Set(0, 2, 0.0);
			dRx_om.Set(1, 0, 0.0);		dRx_om.Set(1, 1, -sin_x);	dRx_om.Set(1, 2, cos_x);
			dRx_om.Set(2, 0, 0.0);		dRx_om.Set(2, 1, -cos_x);	dRx_om.Set(2, 2, -sin_x);

			// --- Phi ---
			
			x		=	eo.dPhi;
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			// Rotational matrix 
			Ry.Set(0, 0, cos_x);		Ry.Set(0, 1, 0.0);			Ry.Set(0, 2, -sin_x);
			Ry.Set(1, 0, 0.0);			Ry.Set(1, 1, 1.0);			Ry.Set(1, 2, 0.0);
			Ry.Set(2, 0, sin_x);		Ry.Set(2, 1, 0.0);			Ry.Set(2, 2, cos_x);

			// Derivative rotational matrix 
			dRy_ph.Set(0, 0, -sin_x);	dRy_ph.Set(0, 1, 0.0);		dRy_ph.Set(0, 2, -cos_x);
			dRy_ph.Set(1, 0, 0.0);		dRy_ph.Set(1, 1, 0.0);		dRy_ph.Set(1, 2, 0.0);
			dRy_ph.Set(2, 0, cos_x);	dRy_ph.Set(2, 1, 0.0);		dRy_ph.Set(2, 2, -sin_x);

			// --- Kappa ---
			
			x		=	eo.dKappa;
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			// Rotational matrix 
			Rz.Set(0, 0, cos_x);		Rz.Set(0, 1, sin_x);		Rz.Set(0, 2, 0.0);
			Rz.Set(1, 0, -sin_x);		Rz.Set(1, 1, cos_x);		Rz.Set(1, 2, 0.0);
			Rz.Set(2, 0, 0.0);			Rz.Set(2, 1, 0.0);			Rz.Set(2, 2, 1.0);

			// Derivative rotational matrix 
			dRz_kp.Set(0, 0, -sin_x);	dRz_kp.Set(0, 1, cos_x);	dRz_kp.Set(0, 2, 0.0);
			dRz_kp.Set(1, 0, -cos_x);	dRz_kp.Set(1, 1, -sin_x);	dRz_kp.Set(1, 2, 0.0);
			dRz_kp.Set(2, 0, 0.0);		dRz_kp.Set(2, 1, 0.0);		dRz_kp.Set(2, 2, 0.0);		
		}

		// Rotational matrix
		if (ROT_ORDER == 1) 
		{
			// R = Rx * Ry * Rz;

			// Derivative rotational matrix wrt Omega
			dR = dRx_om * Ry * Rz;
			m_vec_om_dRotMatrix.push_back(dR);

			// Derivative rotational matrix wrt Phi
			dR = Rx * dRy_ph * Rz;
			m_vec_ph_dRotMatrix.push_back(dR);

			// Derivative rotational matrix wrt Kappa
			dR = Rx * Ry * dRz_kp;
			m_vec_kp_dRotMatrix.push_back(dR);
		}
		else 
		{
			// R = Rz * Ry * Rx;

			// Derivative rotational matrix wrt Omega
			dR = Rz * Ry * dRx_om;
			m_vec_om_dRotMatrix.push_back(dR);

			// Derivative rotational matrix wrt Phi
			dR = Rz * dRy_ph * Rx;
			m_vec_ph_dRotMatrix.push_back(dR);

			// Derivative rotational matrix wrt Kappa
			dR = dRz_kp * Ry * Rx;
			m_vec_kp_dRotMatrix.push_back(dR);
		}
	}
}

/********************************************************************
 * Description:	Construct the indexing vector for m_vec_* wrt. Image ID.
 ********************************************************************/
void CRotationalMatrix::ConstructIndexVector(VEC_EO const &vec_EO_m)
{
	int		i, j;
	int		imageID;
	int		iCurrID;
	EO		eo;

	if (m_vec_idx.size() > 0)
		return;	// The indexing vector has been constructed.

	// Get the number of images in the vector.
	int n_images = vec_EO_m.size();

	// Set current image ID to 0.
	iCurrID = 0;	

	// Iterate through vec_EO_m to fill up the m_vec_idx.
	for (i = 0; i < n_images; i++)
	{
		eo = vec_EO_m[i];

		// Get the image ID
		imageID = eo.iImg;

		// Fill blank ImageID with -1.
		for (j = iCurrID; j < imageID; j++)
		{
			m_vec_idx.push_back(-1);
		}

		m_vec_idx.push_back(i);
		iCurrID = imageID + 1;
	}
}

/********************************************************************
 * Description:	Return the index of the input Image ID.
 ********************************************************************/
int CRotationalMatrix::GetIndexFromImageID(int imageID) const
{
	return m_vec_idx[imageID];
}

/********************************************************************
 * Description:	Return the rotational matrix of the input Image ID.
 ********************************************************************/
Matrix CRotationalMatrix::GetRotMatrix(int imageID) const
{
	int idx = GetIndexFromImageID(imageID);

	return m_vec_RotMatrix[idx];
}

/********************************************************************
 * Description:	Return the derivative rot. matrix for Omega of the input Image ID.
 ********************************************************************/
Matrix CRotationalMatrix::GetDOmegaMatrix(int imageID)
{
	int idx = GetIndexFromImageID(imageID);

	return m_vec_om_dRotMatrix[idx];
}

/********************************************************************
 * Description:	Return the derivative rot. matrix for Phi of the input Image ID.
 ********************************************************************/
Matrix CRotationalMatrix::GetDPhiMatrix(int imageID)
{
	int idx = GetIndexFromImageID(imageID);

	return m_vec_ph_dRotMatrix[idx];
}

/********************************************************************
 * Description:	Return the derivative rot. matrix for Kappa of the input Image ID.
 ********************************************************************/
Matrix CRotationalMatrix::GetDKappaMatrix(int imageID)
{
	int idx = GetIndexFromImageID(imageID);

	return m_vec_kp_dRotMatrix[idx];
}

/********************************************************************
 * Description:	Debug 3D rotational matrix for the specified image ID.
 ********************************************************************/
void CRotationalMatrix::Debug3DRotMatrix(int image_id) const
{
	int r,c;

	// Get index of the input Image ID.
	int idx = GetIndexFromImageID(image_id);

	// Get the 3D rotational matrix of the specified index.
	Matrix R = m_vec_RotMatrix[idx];

	printf("3D Rotational Matrix of The Image Index %d:\n", image_id);

	for (r = 0; r < R.NumRow(); r++)
	{
		for (c = 0; c < R.NumCol(); c++)
		{
			printf("[%d,%d] %lf\t", r, c, R.Get(r,c));
		}
	
		printf("\n");
	}
}

/********************************************************************
 * Description:	Debug derivative rotational matrix for the specified image ID.
 ********************************************************************/
void CRotationalMatrix::DebugDerivativeRotMatrix(int image_id) const
{
	int r,c;
	Matrix R;

	printf("Derivative Rotational Matrix of The Image Index %d:\n", image_id);

	// Get index of the input Image ID.
	int idx = GetIndexFromImageID(image_id);

	// Get the derivative matrix for Omega.
	R = m_vec_om_dRotMatrix[idx];

	printf("Omega\n");

	for (r = 0; r < R.NumRow(); r++)
	{
		for (c = 0; c < R.NumCol(); c++)
		{
			printf("[%d,%d] %lf\t", r, c, R.Get(r,c));
		}

		printf("\n");
	}

	printf("\n");

	// Get the derivative matrix for Phi.
	R = m_vec_ph_dRotMatrix[idx];

	printf("Phi\n");

	for (r = 0; r < R.NumRow(); r++)
	{
		for (c = 0; c < R.NumCol(); c++)
		{
			printf("[%d,%d] %lf\t", r, c, R.Get(r,c));
		}

		printf("\n");
	}

	printf("\n");

	// Get the derivative matrix for Kappa.
	R = m_vec_kp_dRotMatrix[idx];

	printf("Kappa\n");

	for (r = 0; r < R.NumRow(); r++)
	{
		for (c = 0; c < R.NumCol(); c++)
		{
			printf("[%d,%d] %lf\t", r, c, R.Get(r,c));
		}

		printf("\n");
	}
}
