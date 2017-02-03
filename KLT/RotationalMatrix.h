/********************************************************
 * Project:			Automated AT
 * Last updated:	22 August 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _ROTATIONAL_MATRIX_H_
#define	_ROTATIONAL_MATRIX_H_


#pragma once

#include "Definition.h"
#include "Matrix.h"

using namespace std;

typedef vector< Matrix > VEC_Matrix;	// Vector of the Matrix class

class CRotationalMatrix
{
public:
	// Constructor
	CRotationalMatrix();

	// Compute 3D rotational matrix for all EOs in the passed-in vector.
	void Compute3DRotMatrix(VEC_EO const &vec_EO_m, int nimages, int ROT_ORDER, int ROT_TYPE);

	// Compute the derivative of 3D rotational matrix for all EOs.
	void ComputeDerivative3DRotMatrix(VEC_EO const &vec_EO_m, int nimages, int ROT_ORDER, int ROT_TYPE);

	// Construct the indexing vector for m_vec_* wrt. Image ID.
	void ConstructIndexVector(VEC_EO const &vec_EO_m);

	// Return the index of the input Image ID.
	int GetIndexFromImageID(int imageID) const;

	// Return the rotational matrix of the input Image ID.
	Matrix GetRotMatrix(int imageID) const;

	// Return the derivative rot. matrix for Omega of the input Image ID.
	Matrix GetDOmegaMatrix(int imageID);

	// Return the derivative rot. matrix for Phi of the input Image ID.
	Matrix GetDPhiMatrix(int imageID);

	// Return the derivative rot. matrix for Kappa of the input Image ID.
	Matrix GetDKappaMatrix(int imageID);

	// Debug 3D rotational matrix for the specified image.
	void Debug3DRotMatrix(int image_id) const;

	// Debug derivative rotational matrix for the specified image.
	void DebugDerivativeRotMatrix(int image_id) const;

	// Variables declaration.

	VEC_Matrix	m_vec_RotMatrix;		// Vector of the rotational matrix.
	VEC_Matrix	m_vec_om_dRotMatrix;	// Vector of derivative rot. matrix for Omega
	VEC_Matrix	m_vec_ph_dRotMatrix;	// Vector of derivative rot. matrix for Phi
	VEC_Matrix	m_vec_kp_dRotMatrix;	// Vector of derivative rot. matrix for Kappa
	VEC_INT		m_vec_idx;				// Vector containing index in m_vec wrt Image ID.
	int			m_nImages;				// Number of images in the sequence.
};


#endif // _ROTATIONAL_MATRIX_H_
