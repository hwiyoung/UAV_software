/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _AT_SIMAT_H_
#define	_AT_SIMAT_H_


#pragma once

#include "AT_Definition.h"
#include "AT_GP.h"
#include "AT_PhotoCoordinate.h"
#include "AT_RotationalMatrix.h"
#include "AT_SeqSimAT.h"
#include <string>

using namespace std;

class AT_CSimAT
{
public:
	// Constructor with configuration file name.
	AT_CSimAT(const string &_CFG_File);

	// Perform the Simultaneous AT.
	void RunSimAT();

	// Read the configuration file according to the passed-in name.
	void ReadConfig();

	// Read the IO file.
	void ReadIOFile();

	// Read the EO file.
	void ReadEOFile();

	// Read the IP/Tiepoints file.
	void ReadIPFile();

	// Call necessary functions prior to perform sim. AT. 
	void Initialization();

	// Initialize EO parameters e.g. 3D rotational matrix.
	void InitializeEO();

	// Initialize ground points.
	void InitializeGP();

	// Initialize member variables prior to processing.
	void InitializeMembers();


	// Perform Simultaneous AT estimation.
 	void SimATEstimation();

	// A wrapper function to get each line from the input file.
	void GetLine(istream &in, char *buf, const int bufSize) const;

	// Variables declaration.

	string	m_CFG_File;		// Configuration file name.
	string	m_IO_File;		// IO file name.
	string	m_EO_File;		// EO file name.
	string	m_IP_File;		// IP file name.

	int		m_max_iter;		// Maximum number of iterations for performing AT
	int		m_ninit_sim_img;// Initial number of images for starting seq. sim AT.
	double	m_dDelta;		// Delta: stop criteria
	double	m_std_IP;		// std_IP: std. of the image point errors (unit: pixels) 
	double	m_std_GPS;		// std_GPS: std. of the GPS errors (unit: cm) 
	double	m_std_INS;		// std_INS: std. of the INS errors (unit: degrees) 

	AT_IO		m_IO;			// IO / Principle point
	AT_IS		m_IS;			// IO / Image dimension
	AT_PS		m_PS;			// IO / Pixel size
	AT_RDC		m_RDC;			// IO / Radial distortion

	AT_VEC_EO	m_EO_m;			// EO parameters measured from GPS/INS
	AT_VEC_IP	m_IP_i;			// List of tie points read from the input image (Pixel Coord)
	AT_VEC_IP	m_IP_f;			// List of tie points converted to 'Photo Coordinate'

	int		m_nImages;		// Number of images in the sequence

	AT_VEC_INT	m_vec_n_IP;		// Number of points for the specified image ID
							// "vec[1] = 10" means the image ID 1 has 10 image points.
							// The index of vector reflects the real image ID.
							// Therefore, the size of vector is one bigger than the total images.

	int		m_nIP_raw;		// Number of image points for whole sequence including
							// those multiple exist in other images.

	AT_CGP		m_GP;			// CGP object to manage ground points.
	AT_CPhotoCoordinate	m_PhotoCoord;	// CPhotoCoordinate object to manage photo coord.
	AT_CRotationalMatrix	m_RotMatrix;	// CRotationalMatrix object to compute 3D rot. matrix.
	AT_CSeqSimAT			m_SeqSimAT;		// CSeqSimAT object to perform sequential AT simultaneous.

};


#endif // _AT_SIMAT_H_
