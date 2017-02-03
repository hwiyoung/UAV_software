/********************************************************
 * Project:			Automated AT
 * Last updated:	17 Jaunary 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _PHOTOCOORDINATE_H_
#define	_PHOTOCOORDINATE_H_


#pragma once

#include "Definition.h"

using namespace std;

class CPhotoCoordinate
{
public:
	// Constructor
	CPhotoCoordinate();
	CPhotoCoordinate(IO const io, IS const is, PS const ps);

	// Initialize with the value assignment
	void Initialization(IO const io, IS const is, PS const ps);

	// To convert the Pixel Coordinate to Photo Coordinate
	VEC_IP ConvertPixel2PhotoCoords(VEC_IP const &vec_ip_i);

	// To convert the Photo Coordinate to Pixel Coordinate
	VEC_IP ConvertPhoto2PixelCoords(VEC_IP const &vec_ip_f);

	// To convert the number of pixels to metre (m).
	PS ConvertNumPixel2Metre(double nPixels);

	// Variables declaration.

	IO		m_IO;			// IO / Principle point
	IS		m_IS;			// IO / Image dimension
	PS		m_PS;			// IO / Pixel size

	PS		m_ICenter;		// Image center

};


#endif // _PHOTOCOORDINATE_H_
