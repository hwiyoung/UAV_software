/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _AT_ROTATIONAL_MATRIX_H_
#define	_AT_ROTATIONAL_MATRIX_H_


#pragma once

#include "AT_Definition.h"
#include "AT_Matrix.h"

using namespace std;

typedef vector< Matrix > AT_VEC_Matrix;	// Vector of the Matrix class

class AT_CRotationalMatrix
{
public:
	// Constructor
	AT_CRotationalMatrix();

	// Compute 3D rotational matrix for all EOs in the passed-in vector.
	void Compute3DRotMatrix(AT_VEC_EO const &vec_EO_m, int nimages);

	// Compute the derivative of 3D rotational matrix for all EOs.
	void ComputeDerivative3DRotMatrix(AT_VEC_EO const &vec_EO_m, int nimages);

	// Debug 3D rotational matrix for the specified image.
	void Debug3DRotMatrix(int image_id);

	// Debug derivative rotational matrix for the specified image.
	void DebugDerivativeRotMatrix(int image_id);

	// Variables declaration.

	AT_VEC_Matrix	m_vec_RotMatrix;		// Vector of the rotational matrix.
	AT_VEC_Matrix	m_vec_om_dRotMatrix;	// Vector of derivative rot. matrix for Omega
	AT_VEC_Matrix	m_vec_ph_dRotMatrix;	// Vector of derivative rot. matrix for Phi
	AT_VEC_Matrix	m_vec_kp_dRotMatrix;	// Vector of derivative rot. matrix for Kappa
	int			m_nImages;				// Number of images in the sequence.

};


#endif // _AT_ROTATIONAL_MATRIX_H_
