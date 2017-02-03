/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "AT_PhotoCoordinate.h"
#include <iostream>
/********************************************************************
 * Description:	Construction
 ********************************************************************/
AT_CPhotoCoordinate::AT_CPhotoCoordinate()
{
	// No initialization required.
	// Each struct member has implemented its own constructor.
}

AT_CPhotoCoordinate::AT_CPhotoCoordinate(AT_IO const io, AT_IS const is, AT_PS const ps)
{
	m_IO	=	io;		// Principle point
	m_IS	=	is;		// Image dimension
	m_PS	=	ps;		// Pixel size

	// Obtain the image center
	m_ICenter.dX = ((double)is.iWidth - 1.0)/2;
	m_ICenter.dY = ((double)is.iHeight - 1.0)/2;
}

/********************************************************************
 * Description:	Initialize with the value assignment
 ********************************************************************/
void AT_CPhotoCoordinate::Initialization(AT_IO const io, AT_IS const is, AT_PS const ps)
{
	m_IO	=	io;		// Principle point
	m_IS	=	is;		// Image dimension
	m_PS	=	ps;		// Pixel size

	// Obtain the image center
	m_ICenter.dX = ((double)is.iWidth - 1.0)/2;
	m_ICenter.dY = ((double)is.iHeight - 1.0)/2;
	cout << m_ICenter.dX << " " << m_ICenter.dY << endl;
}

/********************************************************************
 * Description:	Convert the pixel coordinate to photo coordinate
 ********************************************************************/
AT_VEC_IP AT_CPhotoCoordinate::ConvertPixel2PhotoCoords(AT_VEC_IP const &vec_ip_i)
{
	int	i, nPoints;
	AT_IP	point;
	AT_VEC_IP	vec_ip_f;		// The photo coordinate vector after conversion.

	// Obtain the number of image points in the vector.
	nPoints = vec_ip_i.size();

	// Calculculate the photo-coordinates for each points (pixel based)
	for (i = 0; i < nPoints; i++)
	{
		point = vec_ip_i[i];

		point.dX = (point.dX - m_ICenter.dX) * m_PS.dX - m_IO.dPPX;
		point.dY = (m_ICenter.dY - point.dY) * m_PS.dY - m_IO.dPPY;

		vec_ip_f.push_back(point);
	}

	return vec_ip_f;
}

/********************************************************************
 * Description:	Convert the photo coordinate to pixel coordinate
 ********************************************************************/
AT_VEC_IP AT_CPhotoCoordinate::ConvertPhoto2PixelCoords(AT_VEC_IP const &vec_ip_f)
{
	int	i, nPoints;
	AT_IP	point;
	AT_VEC_IP	vec_ip_i;		// The pixel coordinate vector after conversion.

	// Obtain the number of image points in the vector.
	nPoints = vec_ip_f.size();

	// Calculculate the pixel-coordinates for each points.
	for (i = 0; i < nPoints; i++)
	{
		point = vec_ip_f[i];

		point.dX = (point.dX + m_IO.dPPX) / m_PS.dX + m_ICenter.dX;
		point.dY = m_ICenter.dY - (point.dY + m_IO.dPPY) / m_PS.dY;

		vec_ip_i.push_back(point);
	}

	return vec_ip_i;
}

/********************************************************************
 * Description:	Convert the number of pixels to metric unit (m)
 ********************************************************************/
AT_PS AT_CPhotoCoordinate::ConvertNumPixel2Metre(double nPixels)
{
	AT_PS	ps;

	// Convert the number of pixels in the x-axis.
	ps.dX = nPixels * m_PS.dX;

	// Convert the number of pixels in the y-axis.
	ps.dY = nPixels * m_PS.dY;
	
	return ps;
}