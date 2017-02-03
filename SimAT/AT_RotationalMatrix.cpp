/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "AT_RotationalMatrix.h"
#include <math.h>

#include "stdio.h"
#include <iostream>
#include <direct.h>   
#include <fstream>		
#include <strstream>
/********************************************************************
 * Description:	Construction
 ********************************************************************/
AT_CRotationalMatrix::AT_CRotationalMatrix()
{
	// No initialization required.
}

/********************************************************************
 * Description:	Compute 3D rotational matrix for all EOs in the 
 *				passed-in vector.
 ********************************************************************/
void AT_CRotationalMatrix::Compute3DRotMatrix(AT_VEC_EO const &vec_EO_m, int nimages)
{
	AT_EO		eo;
	Matrix	R(3,3);			// 3D Rotational Matrix
	Matrix	Rx(3,3), Ry(3,3), Rz(3,3);
	int		i;	
	double	x;				// Angle		
	double	cos_x, sin_x;	// cos(Angle), sin(Angle)

	// Assign the number of images into the member variable.
	m_nImages = nimages;

	// Since the vector index will reflect the real image ID, put the 'dummy' 
	// matrix for the image index 0.
	R.Set(0, 0, 0.0);
	R.Set(0, 1, 0.0);
	R.Set(0, 2, 0.0);
	R.Set(1, 0, 0.0);
	R.Set(1, 1, 0.0);
	R.Set(1, 2, 0.0);
	R.Set(2, 0, 0.0);
	R.Set(2, 1, 0.0);
	R.Set(2, 2, 0.0);
	m_vec_RotMatrix.push_back(R);

	// Compute the rotational matrix for each image.
	for (i = 0; i < nimages; i++)
	{
		eo = vec_EO_m[i];

		//		|	1		0		0		|
		// Rx =	|	0	  cos(Om) sin(Om)	|
		//		|	0	  -sin(Om) cos(Om)	|		
		
		x		=	eo.dOmega;	// Omega
		cos_x	=	cos(x);
		sin_x	=	sin(x);

		Rx.Set(0, 0, 1.0);	Rx.Set(0, 1, 0.0);		Rx.Set(0, 2, 0.0);
		Rx.Set(1, 0, 0.0);	Rx.Set(1, 1, cos_x);	Rx.Set(1, 2, sin_x);
		Rx.Set(2, 0, 0.0);	Rx.Set(2, 1, -sin_x);	Rx.Set(2, 2, cos_x);
		
		//		|	cos(Ph)		0		-sin(Ph)	|
		// Ry =	|	  0			1		   0	|
		//		|	sin(Ph)	0		cos(Ph)	|
		
		x		=	eo.dPhi;	// Phi
		cos_x	=	cos(x);
		sin_x	=	sin(x);

		Ry.Set(0, 0, cos_x);	Ry.Set(0, 1, 0.0);	Ry.Set(0, 2, -sin_x);
		Ry.Set(1, 0, 0.0);		Ry.Set(1, 1, 1.0);	Ry.Set(1, 2, 0.0);
		Ry.Set(2, 0, sin_x);	Ry.Set(2, 1, 0.0);	Ry.Set(2, 2, cos_x);

		//		|	cos(Ka)		sin(Ka)	0	|
		// Rz =	|	-sin(Ka)		cos(Ka)		0	|
		//		|	  0			   0		1	|
		
		x		=	eo.dKappa;	// Kappa
		cos_x	=	cos(x);
		sin_x	=	sin(x);

		Rz.Set(0, 0, cos_x);	Rz.Set(0, 1, sin_x);	Rz.Set(0, 2, 0.0);
		Rz.Set(1, 0, -sin_x);	Rz.Set(1, 1, cos_x);	Rz.Set(1, 2, 0.0);
		Rz.Set(2, 0, 0.0);		Rz.Set(2, 1, 0.0);		Rz.Set(2, 2, 1.0);

		// Rotational matrix
		//R = Rx * Ry * Rz;
		R = Rz * Ry * Rx;

		// Store the rotational matrix into the vector.
		m_vec_RotMatrix.push_back(R);
	}

	//ofstream R_res("Results\\AT_Results\\R_res.txt");

	//for(int j =0 ; j<nimages-1 ; j++) {
	//	R_res << j << "\n";
	//	R_res << m_vec_RotMatrix[j].at(0) << "\t" << m_vec_RotMatrix[j].at(1) << "\t"<< m_vec_RotMatrix[j].at(2) << "\n";
	//	R_res << m_vec_RotMatrix[j].at(3) << "\t" << m_vec_RotMatrix[j].at(4) << "\t"<< m_vec_RotMatrix[j].at(5) << "\n";
	//	R_res << m_vec_RotMatrix[j].at(6) << "\t" << m_vec_RotMatrix[j].at(7) << "\t"<< m_vec_RotMatrix[j].at(8) << endl;
	//}
}

/********************************************************************
 * Description:	Compute the derivative of 3D rotational matrix for  
 *				all EOs in the passed-in vector.
 ********************************************************************/
void AT_CRotationalMatrix::ComputeDerivative3DRotMatrix(AT_VEC_EO const &vec_EO_m, int nimages)
{
	AT_EO		eo;
	Matrix	R(3,3);			// 3D Rotational Matrix
	Matrix	Rx(3,3), Ry(3,3), Rz(3,3);
	Matrix	dRx(3,3), dRy(3,3), dRz(3,3);
	int		i;	
	double	om, ph, kp;		// Angle		
	double	cos_om, sin_om, cos_ph, sin_ph, cos_kp, sin_kp;	// cos(Angle), sin(Angle)

	// Since the vector index will reflect the real image ID, put the 'dummy' 
	// matrix for the image index 0.
	R.Set(0, 0, 0.0);	R.Set(0, 1, 0.0);	R.Set(0, 2, 0.0);
	R.Set(1, 0, 0.0);	R.Set(1, 1, 0.0);	R.Set(1, 2, 0.0);
	R.Set(2, 0, 0.0);	R.Set(2, 1, 0.0);	R.Set(2, 2, 0.0);

	m_vec_om_dRotMatrix.push_back(R);
	m_vec_ph_dRotMatrix.push_back(R);
	m_vec_kp_dRotMatrix.push_back(R);

	// Compute the derivative rotational matrix for each image.
	for (i = 0; i < nimages; i++)
	{
		eo = vec_EO_m[i];

		// Omega
		om		=	eo.dOmega;
		cos_om	=	cos(om);
		sin_om	=	sin(om);

		// Phi
		ph		=	eo.dPhi;
		cos_ph	=	cos(ph);
		sin_ph	=	sin(ph);

		// Kappa
		kp		=	eo.dKappa;
		cos_kp	=	cos(kp);
		sin_kp	=	sin(kp);		
		
		Rx.Set(0, 0, 1.0);	Rx.Set(0, 1, 0.0);		Rx.Set(0, 2, 0.0);
		Rx.Set(1, 0, 0.0);	Rx.Set(1, 1, cos_om);	Rx.Set(1, 2, sin_om);
		Rx.Set(2, 0, 0.0);	Rx.Set(2, 1, -sin_om);	Rx.Set(2, 2, cos_om);		

		Ry.Set(0, 0, cos_ph);	Ry.Set(0, 1, 0.0);	Ry.Set(0, 2, -sin_ph);
		Ry.Set(1, 0, 0.0);		Ry.Set(1, 1, 1.0);	Ry.Set(1, 2, 0.0);
		Ry.Set(2, 0, sin_ph);	Ry.Set(2, 1, 0.0);	Ry.Set(2, 2, cos_ph);

		Rz.Set(0, 0, cos_kp);	Rz.Set(0, 1, sin_kp);	Rz.Set(0, 2, 0.0);
		Rz.Set(1, 0, -sin_kp);	Rz.Set(1, 1, cos_kp);	Rz.Set(1, 2, 0.0);
		Rz.Set(2, 0, 0.0);		Rz.Set(2, 1, 0.0);		Rz.Set(2, 2, 1.0);

		dRx.Set(0, 0, 0.0);	dRx.Set(0, 1, 0.0);		dRx.Set(0, 2, 0.0);
		dRx.Set(1, 0, 0.0);	dRx.Set(1, 1, -sin_om);	dRx.Set(1, 2, cos_om);
		dRx.Set(2, 0, 0.0);	dRx.Set(2, 1, -cos_om);	dRx.Set(2, 2, -sin_om);

		dRy.Set(0, 0, -sin_ph); dRy.Set(0, 1, 0.0); dRy.Set(0, 2, -cos_ph);
		dRy.Set(1, 0, 0.0);	    dRy.Set(1, 1, 0.0); dRy.Set(1, 2, 0.0);
		dRy.Set(2, 0, cos_ph);  dRy.Set(2, 1, 0.0); dRy.Set(2, 2, -sin_ph);

		dRz.Set(0, 0, -sin_kp); dRz.Set(0, 1, cos_kp);	dRz.Set(0, 2, 0.0);
		dRz.Set(1, 0, -cos_kp); dRz.Set(1, 1, -sin_kp);	dRz.Set(1, 2, 0.0);
		dRz.Set(2, 0, 0.0);     dRz.Set(2, 1, 0.0);		dRz.Set(2, 2, 0.0);

		R = Rz * Ry * dRx ;
		m_vec_om_dRotMatrix.push_back(R);

		R = Rz * dRy * Rx;
		m_vec_ph_dRotMatrix.push_back(R);

		R = dRz * Ry* Rx;
		m_vec_kp_dRotMatrix.push_back(R);
	}

	//ofstream R_resDom("Results\\AT_Results\\R_resDom.txt");
	//ofstream R_resDph("Results\\AT_Results\\R_resDph.txt");
	//ofstream R_resDkp("Results\\AT_Results\\R_resDkp.txt");

	//for(int j =0 ; j<nimages ; j++) {
	//	R_resDom << j << "\n";
	//	R_resDom << m_vec_om_dRotMatrix[j].at(0) << "\t" << m_vec_om_dRotMatrix[j].at(1) << "\t"<< m_vec_om_dRotMatrix[j].at(2) << "\n";
	//	R_resDom << m_vec_om_dRotMatrix[j].at(3) << "\t" << m_vec_om_dRotMatrix[j].at(4) << "\t"<< m_vec_om_dRotMatrix[j].at(5) << "\n";
	//	R_resDom << m_vec_om_dRotMatrix[j].at(6) << "\t" << m_vec_om_dRotMatrix[j].at(7) << "\t"<< m_vec_om_dRotMatrix[j].at(8) << endl;
	//}

	//for(int j =0 ; j<nimages ; j++) {
	//	R_resDph << j << "\n";
	//	R_resDph << m_vec_ph_dRotMatrix[j].at(0) << "\t" << m_vec_ph_dRotMatrix[j].at(1) << "\t"<< m_vec_ph_dRotMatrix[j].at(2) << "\n";
	//	R_resDph << m_vec_ph_dRotMatrix[j].at(3) << "\t" << m_vec_ph_dRotMatrix[j].at(4) << "\t"<< m_vec_ph_dRotMatrix[j].at(5) << "\n";
	//	R_resDph << m_vec_ph_dRotMatrix[j].at(6) << "\t" << m_vec_ph_dRotMatrix[j].at(7) << "\t"<< m_vec_ph_dRotMatrix[j].at(8) << endl;
	//}


	//for(int j =0 ; j<nimages ; j++) {
	//	R_resDkp << j << "\n";
	//	R_resDkp << m_vec_kp_dRotMatrix[j].at(0) << "\t" << m_vec_kp_dRotMatrix[j].at(1) << "\t"<< m_vec_kp_dRotMatrix[j].at(2) << "\n";
	//	R_resDkp << m_vec_kp_dRotMatrix[j].at(3) << "\t" << m_vec_kp_dRotMatrix[j].at(4) << "\t"<< m_vec_kp_dRotMatrix[j].at(5) << "\n";
	//	R_resDkp << m_vec_kp_dRotMatrix[j].at(6) << "\t" << m_vec_kp_dRotMatrix[j].at(7) << "\t"<< m_vec_kp_dRotMatrix[j].at(8) << endl;
	//}
}

/********************************************************************
 * Description:	Debug 3D rotational matrix for the specified image ID.
 ********************************************************************/
void AT_CRotationalMatrix::Debug3DRotMatrix(int image_id)
{
	int r,c;

	// Return if the specified image index is invalid.
	if (image_id > m_nImages)
	{
		printf("CRotationalMatrix::Debug3DRotMatrix()\n");
		printf("\tThe image index %d is invalid since the largest index is %d.\n", image_id, m_nImages);
		return;
	}

	// Get the 3D rotational matrix of the specified index.
	Matrix R = m_vec_RotMatrix[image_id];

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
void AT_CRotationalMatrix::DebugDerivativeRotMatrix(int image_id)
{
	int r,c;
	Matrix R;

	// Return if the specified image index is invalid.
	if (image_id > m_nImages)
	{
		printf("CRotationalMatrix::DebugDerivativeRotMatrix()\n");
		printf("\tThe image index %d is invalid since the largest index is %d.\n", image_id, m_nImages);
		return;
	}

	printf("Derivative Rotational Matrix of The Image Index %d:\n", image_id);

	// Get the derivative matrix for Omega.
	R = m_vec_om_dRotMatrix[image_id];

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
	R = m_vec_ph_dRotMatrix[image_id];

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
	R = m_vec_kp_dRotMatrix[image_id];

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
