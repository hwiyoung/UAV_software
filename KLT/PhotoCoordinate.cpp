/********************************************************
 * Project:			Automated AT
 * Last updated:	17 Jaunary 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "PhotoCoordinate.h"

/********************************************************************
 * Description:	Construction
 ********************************************************************/
CPhotoCoordinate::CPhotoCoordinate()
{
	// No initialization required.
	// Each struct member has implemented its own constructor.
}

CPhotoCoordinate::CPhotoCoordinate(IO const io, IS const is, PS const ps)
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
void CPhotoCoordinate::Initialization(IO const io, IS const is, PS const ps)
{
	m_IO	=	io;		// Principle point
	m_IS	=	is;		// Image dimension
	m_PS	=	ps;		// Pixel size

	// Obtain the image center
	m_ICenter.dX = ((double)is.iWidth - 1.0)/2;
	m_ICenter.dY = ((double)is.iHeight - 1.0)/2;
}

/********************************************************************
 * Description:	Convert the pixel coordinate to photo coordinate
 ********************************************************************/
VEC_IP CPhotoCoordinate::ConvertPixel2PhotoCoords(VEC_IP const &vec_ip_i)
{
	int	i, nPoints;
	IP	point;
	VEC_IP	vec_ip_f;		// The photo coordinate vector after conversion.

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
VEC_IP CPhotoCoordinate::ConvertPhoto2PixelCoords(VEC_IP const &vec_ip_f)
{
	int	i, nPoints;
	IP	point;
	VEC_IP	vec_ip_i;		// The pixel coordinate vector after conversion.

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
PS CPhotoCoordinate::ConvertNumPixel2Metre(double nPixels)
{
	PS	ps;

	// Convert the number of pixels in the x-axis.
	ps.dX = nPixels * m_PS.dX;

	// Convert the number of pixels in the y-axis.
	ps.dY = nPixels * m_PS.dY;
	
	return ps;
}