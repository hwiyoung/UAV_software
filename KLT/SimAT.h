/********************************************************
 * Project:			Automated AT
 * Last updated:	17 Jaunary 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _SIMAT_H_
#define	_SIMAT_H_


#pragma once

#include "Definition.h"
#include "GP.h"
#include "PhotoCoordinate.h"
#include "RotationalMatrix.h"
#include "SeqSimAT.h"
#include "KLTTracker.h"
#include <string>

using namespace std;

class CSimAT
{
public:
	// Constructor with configuration file name.
	CSimAT(const string &_CFG_File);

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

	// Read the IP/Tiepoints file generated from Simulator.
	void ReadIPFile2();

	// Call necessary functions prior to perform sim. AT. 
	void Initialization();

	// Initialize the output tie point file.
	void InitializeTPFile();

	// Initialize EO parameters e.g. 3D rotational matrix.
	void InitializeEO();

	// Initialize ground points.
	void InitializeGP();

	// Initialize member variables prior to processing.
	void InitializeMembers();

	// Convert the units of input parameters to appropriate scale.
	void UnitConversion();

	// Convert the units of std. deviation of error to appropriate scale.
	void StdUnitConversion();

	// Perform Simultaneous AT estimation.
 	void SimATEstimation();

	// This function will be called whenever a new image is captured. 
	void GetNewImage(string image_name, EO eo);

	// Store the pointer to the singleton CKLTTracker object.
	void RegisterKLT(CKLTTracker * pKLT);

	// Get the set of tie-points according to epoch.
	void GetTiePoints();

	// Print out the result of the Simultaneous AT estimation.
	void PrintOutResults();

	// A wrapper function to get each line from the input file.
	void GetLine(istream &in, char *buf, const int bufSize) const;

	// Variables declaration.

	bool	m_bRun;			// Whether or not to run this process.
	bool	m_bOutputTP;	// Whether or not to output the tie points.

	string	m_CFG_File;		// Configuration file name.
	string	m_IO_File;		// IO file name.
	string	m_EO_File;		// EO file name.
	string	m_IP_File;		// IP file name.

	int		m_max_iter;		// Maximum number of iterations for performing AT
	int		m_ninit_sim_img;// Initial number of images for starting seq. sim AT.
	int		m_nmax_sim_img;	// Maximum number of images for seq. AT computation.

	double	m_dDelta;		// Delta: stop criteria
	double	m_std_IP;		// std_IP: std. of the image point errors (unit: pixels) 
	double	m_std_GPS;		// std_GPS: std. of the GPS errors (unit: cm) 
	double	m_std_INS;		// std_INS: std. of the INS errors (unit: degrees) 

	IO		m_IO;			// IO / Principle point
	IS		m_IS;			// IO / Image dimension
	PS		m_PS;			// IO / Pixel size
	RDC		m_RDC;			// IO / Radial distortion

	int		ROT_ORDER;		// Rotational matrix order
							// 1: Rx*Ry*Rz or 2: Rz*Ry*Rx
	int		ROT_TYPE;		// Rotational matrix type: 1 or 2
							// Type 1 
							//			Rx = [1 0 0; 0 cos(om) -sin(om); 0 sin(om) cos(om)]
							//			Ry = [cos(ph) 0 sin(ph); 0 1 0; -sin(ph) 0 cos(ph)]
							//			Rz = [cos(kp) -sin(kp) 0; sin(kp) cos(kp) 0; 0 0 1]
							// Type 2
							//			Rx = [1 0 0; 0 cos(om) sin(om); 0 -sin(om) cos(om)]
							//			Ry = [cos(ph) 0 -sin(ph); 0 1 0; sin(ph) 0 cos(ph)]
							//			Rz = [cos(kp) sin(kp) 0; -sin(kp) cos(kp) 0; 0 0 1]

	VEC_EO	m_EO_m;			// EO parameters measured from GPS/INS
	VEC_IP	m_IP_i;			// List of tie points read from the input image (Pixel Coord)
	VEC_IP	m_IP_f;			// List of tie points converted to 'Photo Coordinate'
	VEC_IP	r_IP_i;
	VEC_IP	r_IP_f;

	int		m_nImages;		// Number of images in the sequence

	VEC_INT	m_vec_n_IP;		// Number of points for the specified image ID
							// "vec[1] = 10" means the image ID 1 has 10 image points.
							// The index of vector reflects the real image ID.
							// Therefore, the size of vector is one bigger than the total images.

	VEC_INT	r_vec_n_IP;

	VEC_INT	m_vec_n_EPOCH;	// Number of points for each epoch ID.
							// The epoch ID starts from 0. For each epoch, it is not necessary
							// to contains tie points from a single image. 

	int		m_nIP_raw;		// Number of image points for whole sequence including
							// those multiple exist in other images.

	CGP		m_GP;			// CGP object to manage ground points.
	CPhotoCoordinate	m_PhotoCoord;	// CPhotoCoordinate object to manage photo coord.
	CRotationalMatrix	m_RotMatrix;	// CRotationalMatrix object to compute 3D rot. matrix.
	CSeqSimAT			m_SeqSimAT;		// CSeqSimAT object to perform sequential AT simultaneous.

	// Provider object
	CKLTTracker *	m_pKLT;	// Pointer to the singleton KLT Tracker object.

	int		m_EPOCH_ID;		// The running epoch ID starting from 0.

};


#endif // _SIMAT_H_
