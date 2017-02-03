/********************************************************
 * Project:			Automated AT
 * Last updated:	24 Jaunary 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _AERIAL_SIMULATOR_H_
#define	_AERIAL_SIMULATOR_H_


#pragma once

#include "Definition.h"
#include "Global.h"
#include "Array.h"
#include "Point.h"
#include "Matrix.h"
#include <CRTDBG.H>
#include "KLTTracker.h"
#include "SimAT.h"

using namespace std;

class CAerialSimulator
{
public:
	// Constructor with configuration file name.
	CAerialSimulator(const string &_CFG_File);

	// Main function to start the Aerial section.
	void StartAerialSimulator();

	// Call necessary functions for initialization.
	void Initialization();

	// Set member variables.
	void SetMembers();

	// Read the configuration file
	void ReadConfig();

	// Read the EO file
	void ReadEOFile();

	// Return the period of time to capture a new image (second)
	double GetExposurePeriod();

	// Store the pointer to the singleton CKLTTracker object.
	void RegisterKLTClient(CKLTTracker * pKLT);

	// Store the pointer to the singleton CSimAT object.
	void RegisterSimAT(CSimAT * pSimAT);

	// Simulate that a new image is captured.
	void RunImageExposure();

	// Handle the process once the exposure event fired.
	void OnExposureEventFired(string image_name, EO eo);
		// Starting from v8.0
	void OnExposureEventFired(string image_name, EO eo, double za);

	// A wrapper function to get each line from the input file.
	void GetLine(istream &in, char *buf, const int bufSize) const;

	// Variables declaration.

	string		m_CFG_File;		// Configuration file name.
	string		m_EO_File;		// EO file name.
	string		m_Za_File;		// Za file name (if ZA_FILE_EXIST=1)
	string		m_img_directory;// Image directory.
	int			ZA_FILE_EXIST;	// If "Za file" exist.
	double		CAPTURE_PERIOD;	// Time period to capture a new image (sec)
	int			TOTAL_IMAGE;	// Total number of images in the sequence.
	DEQ_STRING	m_deq_IMG_File;	// Deque storing the image file name.
	DEQ_EO		m_deq_EO;		// Deque storing the EO parameters.
	DEQ_DOUBLE	m_deq_Za;		// Deque storing Za for each image.
	int			RUNNING_INDEX;	// The index ID of running image.

	// Client objects
	CKLTTracker *	m_pKLT;		// Pointer to the singleton KLT Tracker object.
	CSimAT *		m_pSimAT;	// Pointer to the singleton SimAT object.

};

#endif // _AERIAL_SIMULATOR_H_
