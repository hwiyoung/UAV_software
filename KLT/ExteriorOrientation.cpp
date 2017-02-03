/********************************************************
 * Project:			Automated AT
 * Last updated:	26 August 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "GP.h"
#include "ExteriorOrientation.h"
#include "Matrix.h"
#include "stdio.h"
#include <iostream>
#include <CRTDBG.H>	// To use _ASSERT

/********************************************************************
 * Description:	Construction
 ********************************************************************/
CExteriorOrientation::CExteriorOrientation()
{
	m_nImages = 0;
}

CExteriorOrientation::CExteriorOrientation(VEC_EO const &vec_eo, int &n_images, bool bCopy, int base_index)
{
	m_nImages = 0;

	SetMember(vec_eo, n_images, bCopy, base_index);

	ConstructIndexVector();
}

/********************************************************************
 * Description:	Set the member parameters.
 ********************************************************************/
void CExteriorOrientation::SetMember(VEC_EO const &vec_eo, int &n_images)
{
	int i;
	EO	eo;

	// Store the EO parameters into its member variable.
	int vec_eoSize=0;
	vec_eoSize=vec_eo.size();

	for (i = 0; i < vec_eoSize; i++)
	{
		eo = vec_eo[i];
		m_vec_eo.push_back(eo);
	}

	// Store the number of images in the sequence to the member variable.
	m_nImages = n_images;
}

void CExteriorOrientation::SetMember(VEC_EO const &vec_eo, int &n_images, bool bCopy, int base_index)
{
	// Copy the whole EO into the member variable.
	if (bCopy == true)
	{
		SetMember(vec_eo, n_images);
	}
	// Store the images for the specified number of images n_images starting from base_index.
	else
	{
		int i;
		EO	eo;
		
		// Get the index of the specified based image ID.
		int idx_base_vec = GetEOIndex(vec_eo, base_index);

		// Get the last index used in the computation.
		int vec_eoSize=0;
		vec_eoSize=vec_eo.size();

		int idx_last_vec = ( idx_base_vec+n_images > vec_eoSize )? vec_eo.size():idx_base_vec+n_images;

		// Store the EO parameters into its member variable.
		for (i = idx_base_vec; i < idx_last_vec; i++)
		{
			eo = vec_eo[i];
			m_vec_eo.push_back(eo);
			m_nImages++;
		}

		// Assert to validate the accuracy of program.
		_ASSERT( m_nImages == n_images );

	}

	// Store the based Image ID.
	if (base_index != 0)
	{
		m_base_index = base_index;
	}
	else
	{
		EO eo = m_vec_eo[0];
		m_base_index = eo.iImg;
	}
}

/********************************************************************
 * Description:	Reset the member parameters.
 ********************************************************************/
void CExteriorOrientation::ResetMember(VEC_EO const &vec_eo, int &n_images, bool bCopy, int base_index)
{
	// Clear the member variable.
	m_vec_eo.clear();		// Clear the EO parameters of all images.
	m_vec_idx.clear();		// Clear the index vector.
	m_nImages = 0;			// Reset the number of images.
	m_base_index = 0;		// Reset the base index of image.

	// Store the passed-in EO to its member variable.
	SetMember(vec_eo, n_images, bCopy, base_index);

	ConstructIndexVector();
}

/********************************************************************
 * Description:	Construct the indexing vector for m_vec_eo wrt. Image ID.
 ********************************************************************/
void CExteriorOrientation::ConstructIndexVector()
{
	int		i, j;
	int		imageID;
	int		iCurrID;
	EO		eo;

	if (m_vec_idx.size() > 0)
		return;	// The indexing vector has been constructed.

	// Get the number of images in the vector.
	int n_images = m_vec_eo.size();

	// Set current image ID to 0.
	iCurrID = 0;	

	// Iterate through vec_EO_m to fill up the m_vec_idx.
	for (i = 0; i < n_images; i++)
	{
		eo = m_vec_eo[i];

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
 * Description:	Return the index of the EO of the specific image from the passed-in vector.
 ********************************************************************/
int CExteriorOrientation::GetEOIndex(VEC_EO const &vec_eo, int imageID)
{
	int i;
	int vec_eoSize=0;
	vec_eoSize=vec_eo.size();

	for (i = 0; i < vec_eoSize; i++)
	{
		if (imageID == vec_eo[i].iImg)
			break;
	}

	if (i == vec_eo.size())
	{
		cout << "Error! The EO of the specified image does not exist." << endl;
		exit(0);
	}

	return i;
}

/********************************************************************
 * Description:	Return the index of the EO of the specific image.
 ********************************************************************/
int CExteriorOrientation::GetEOIndex(int imageID)
{
	return m_vec_idx[imageID];
}

/********************************************************************
 * Description:	Return the EO of the specific image from the passed-in vector.
 ********************************************************************/
EO CExteriorOrientation::GetEO(VEC_EO const &vec_eo, int imageID)
{
	int idx = GetEOIndex(vec_eo, imageID);
	
	return m_vec_eo[idx];
}

/********************************************************************
 * Description:	Return the EO of the specific image.
 ********************************************************************/
EO CExteriorOrientation::GetEO(int imageID)
{
	int idx = GetEOIndex(imageID);
	
	return m_vec_eo[idx];
}
