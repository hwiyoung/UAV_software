/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _AT_PHOTOCOORDINATE_H_
#define	_AT_PHOTOCOORDINATE_H_


#pragma once

#include "AT_Definition.h"

using namespace std;

class AT_CPhotoCoordinate
{
public:
	// Constructor
	AT_CPhotoCoordinate();
	AT_CPhotoCoordinate(AT_IO const io, AT_IS const is, AT_PS const ps);

	// Initialize with the value assignment
	void Initialization(AT_IO const io, AT_IS const is, AT_PS const ps);

	// To convert the Pixel Coordinate to Photo Coordinate
	AT_VEC_IP ConvertPixel2PhotoCoords(AT_VEC_IP const &vec_ip_i);

	// To convert the Photo Coordinate to Pixel Coordinate
	AT_VEC_IP ConvertPhoto2PixelCoords(AT_VEC_IP const &vec_ip_f);

	// To convert the number of pixels to metre (m).
	AT_PS ConvertNumPixel2Metre(double nPixels);

	// Variables declaration.

	AT_IO		m_IO;			// IO / Principle point
	AT_IS		m_IS;			// IO / Image dimension
	AT_PS		m_PS;			// IO / Pixel size

	AT_PS		m_ICenter;		// Image center

};


#endif // _AT_PHOTOCOORDINATE_H_
