/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "AT_GP.h"
#include "AT_ExteriorOrientation.h"
#include "AT_Matrix.h"
#include "stdio.h"
#include <iostream>


/********************************************************************
 * Description:	Construction
 ********************************************************************/
AT_CExteriorOrientation::AT_CExteriorOrientation()
{
	// No operations.
}

AT_CExteriorOrientation::AT_CExteriorOrientation(AT_VEC_EO const &vec_eo, int &n_images, bool bCopy)
{
	SetMember(vec_eo, n_images, bCopy);
}

/********************************************************************
 * Description:	Set the member parameters.
 ********************************************************************/
void AT_CExteriorOrientation::SetMember(AT_VEC_EO const &vec_eo, int &n_images)
{
	int i;
	AT_EO	eo;

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

void AT_CExteriorOrientation::SetMember(AT_VEC_EO const &vec_eo, int &n_images, bool bCopy)
{
	// Copy the whole EO into the member variable.
	if (bCopy == true)
	{
		SetMember(vec_eo, n_images);
	}
	// Store only the images in which its ID is less or equal the maximum n_images.
	else
	{
		int i;
		AT_EO	eo;
		
		m_nImages = 0;

		// Store the EO parameters into its member variable.
		int vec_eoSize=0;
		vec_eoSize=vec_eo.size();
		
		for (i = 0; i < vec_eoSize; i++)
		{
			eo = vec_eo[i];

			if (eo.iImg <= n_images)
			{
				m_vec_eo.push_back(eo);
				m_nImages++;
			}
		}
	}
}

/********************************************************************
 * Description:	Return the index of the EO of the specific image.
 ********************************************************************/
int AT_CExteriorOrientation::GetEOIndex(int imageID)
{
	int i;
	int vec_eoSize=0;
	vec_eoSize=m_vec_eo.size();
	
	for (i = 0; i < vec_eoSize; i++)
	{
		if (imageID == m_vec_eo[i].iImg)
			break;
	}

	if (i == vec_eoSize)
	{
		cout << "Error! The EO of the specified image does not exist." << endl;
		exit(0);
	}

	return i;
}

/********************************************************************
 * Description:	Return the EO of the specific image.
 ********************************************************************/
AT_EO AT_CExteriorOrientation::GetEO(int imageID)
{
	int idx = GetEOIndex(imageID);
	
	return m_vec_eo[idx];
}
