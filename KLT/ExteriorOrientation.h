/********************************************************
 * Project:			Automated AT
 * Last updated:	17 Jaunary 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _EXTERIOR_ORIENTATION_H_
#define	_EXTERIOR_ORIENTATION_H_


#pragma once

#include "Definition.h"
#include "RotationalMatrix.h"

using namespace std;

class CExteriorOrientation
{
public:
	// Constructor
	CExteriorOrientation();

	// Constructor
	CExteriorOrientation(VEC_EO const &vec_eo, int &n_images, bool bCopy, int base_index=0);

	// 'Copy' the passed-in parameters to the member variables.
	void SetMember(VEC_EO const &vec_eo, int &n_images);

	// Set the member parameters.
	void SetMember(VEC_EO const &vec_eo, int &n_images, bool bCopy, int base_index=0);

	// Reset the member parameters.
	void ResetMember(VEC_EO const &vec_eo, int &n_images, bool bCopy, int base_index=0);

	// Construct the indexing vector for m_vec_eo wrt. Image ID.
	void ConstructIndexVector();

	// Return the EO of the specific image.
	EO GetEO(int imageID);	
	
	// Return the index of the EO of the specific image from the passed-in vector.
	int GetEOIndex(VEC_EO const &vec_eo, int imageID);

	// Return the index of the EO of the specific image.
	int GetEOIndex(int imageID);

	// Return the EO of the specific image from the passed-in vector.
	EO GetEO(VEC_EO const &vec_eo, int imageID);

	// Variables declaration.
	VEC_EO	m_vec_eo;		// Vector of EO parameters for each image.
	VEC_INT	m_vec_idx;		// Vector containing index in m_vec wrt Image ID.
	int		m_nImages;		// Number of images in the sequence that reflects "size_eo + 1"
	int		m_base_index;	// The based Image ID for the m_vec_eo vector. 

};

typedef vector<CExteriorOrientation>	VEC_CExteriorOrientation;	// Vector of the CExteriorOrientation object

#endif // _EXTERIOR_ORIENTATION_H_
