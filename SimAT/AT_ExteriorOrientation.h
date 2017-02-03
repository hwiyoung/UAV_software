/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _AT_EXTERIOR_ORIENTATION_H_
#define	_AT_EXTERIOR_ORIENTATION_H_


#pragma once

#include "AT_Definition.h"
#include "AT_RotationalMatrix.h"

using namespace std;

class AT_CExteriorOrientation
{
public:
	// Constructor
	AT_CExteriorOrientation();

	// Constructor
	AT_CExteriorOrientation(AT_VEC_EO const &vec_eo, int &n_images, bool bCopy);

	// 'Copy' the passed-in parameters to the member variables.
	void SetMember(AT_VEC_EO const &vec_eo, int &n_images);

	// Set the member parameters.
	void SetMember(AT_VEC_EO const &vec_eo, int &n_images, bool bCopy);

	// Return the EO of the specific image.
	AT_EO GetEO(int imageID);			

	// Return the index of the EO of the specific image.
	int GetEOIndex(int imageID);

	// Variables declaration.
	AT_VEC_EO	m_vec_eo;		// Vector of EO parameters for each image.
	int		m_nImages;		// Number of images in the sequence that reflects "size_eo + 1"

};

typedef vector<AT_CExteriorOrientation>	AT_VEC_CExteriorOrientation;	// Vector of the CExteriorOrientation object

#endif // _EXTERIOR_ORIENTATION_H_
