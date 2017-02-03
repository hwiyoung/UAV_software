/********************************************************
 * Project:			Automated AT
 * Last updated:	29 August 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _KLTTRACKER_H_
#define	_KLTTRACKER_H_


#pragma once

#include <cv.h>
#include <cxcore.h>
#include <highgui.h>
#include <stdio.h>
#include <vector>
#include <io.h>
#include <time.h>
#include <windows.h>
#include <iostream>
#include <fstream>
#include <strstream>
#include <string>
#include "Definition.h"
#include "Global.h"
#include "Array.h"
#include "Point.h"
#include "Matrix.h"
#include "cvfast.h"
#include <CRTDBG.H>
//#include "GP.h"
//#include "PhotoCoordinate.h"
//#include "RotationalMatrix.h"

using namespace std;

class CKLTTracker
{
public:
	// Constructor with configuration file name.
	CKLTTracker(const string &_CFG_File);

	// Destructor
	~CKLTTracker();

	// Set member variables
	void SetMembers(int n_init_num_images);

	// Read the KLT configuration file
	void ReadConfig();

	// Read the IO file
	void ReadIOFile();

	// A wrapper function to get each line from the input file.
	void GetLine(istream &in, char *buf, const int bufSize) const;

	// It performs all preliminary tasks prior to run the KLT tracking.
	void Startup();

	// Initialize the constant values required for KLT tracker.
	void InitializeConstants();

	// Initialize the member variables.
	void InitializeMembers();

	// Find subregion ID from the passed-in coordinate.
	int FindSubRegionFromCoord(float y_coord, float x_coord);
	int FindSubRegionFromCoord(float y_coord, float x_coord, /*out*/int &iH, /*out*/int &iW);

	// Find subregion ID from the passed-in row/column ID.
	int FindSubRegionFromID(int iH, int iW);

	// Find the row/column ID of the passed-in region ID.
	void FindRowColFromRegionID(int regionID, /*out*/int &iH, /*out*/int &iW); 

	// Define the ROI of the specified sub-pattern.
	void DefineROI(	int iH, int iW, 
					/*out*/ int &beginX, /*out*/ int &beginY,
					/*out*/ int &lengthX, /*out*/ int &lengthY);

	// Check if at least one block does not locate any tie point.
	bool IsBlankRegionExisted();
	bool IsBlankRegionExisted(int imageID, int * existing_features);
	bool IsBlankRegionExisted(int imageID, int * existing_features, int &n_existing_features);

	// Store existing features into each block of sub-pattern.
	void AssignExistingFeaturesToBlock(	int imageID, 
										vector<vector<VEC_INT> > &vblock_features, 
										int * existing_features, 
										int &n_existing_features);

	// To replace out-of-scope points with valid point from previous image.
	int ReplaceOutOfScopePoints(CvPoint2D32f * cornersA, 
								CvPoint2D32f * cornersB, 
								int * cornersA_GPID,
								int &nPoint);	

	// To remove out-of-scope points.
	int RemoveOutOfScopePoints(	CvPoint2D32f * cornersA, 
								CvPoint2D32f * cornersB, 
								int * cornersA_GPID,
								int &nPoint);	

	// Check if the passed-in point already existed.
	bool IsMatchingPointsExisted(int imageID, CvPoint2D32f point);

	// To add the signficant features to the list of features to be tracked.
	void AddMatchedPointsForTracking(int imageID, CvPoint2D32f * cornersA, int &nPoint);
	void AddMatchedPointsForTracking(	int imageID, 
										CvPoint2D32f * cornersA, 
										int &nPoint, 
										int * cornersA_GPID);

	// To create rotation matrix for both first and second image.
	void Create3DRotMatrix();

	// To calculate the initial guessed positions of the tie points based on collinearity eq.
	void CalInitialGuesses(int corner_count, CvPoint2D32f * cornersA, /*out*/CvPoint2D32f * cornersB);
	void CalInitialGuesses(	int corner_count, CvPoint2D32f * cornersA, /*out*/CvPoint2D32f * cornersB,
							double rot_a, double rot_b, double rot_tx, double rot_ty);

	// To convert the pixel coordinates to the photo coordinates.
	void ConvertPix2PhotoCoords(int corner_count,	
								CvPoint2D32f * cornersA, 
								/*out*/CvPoint2D32f * arrCoords);

	// To convert the pixel coord. to photo coord. for the passed-in pixel.
	CvPoint2D32f ConvertPix2PhotoCoord(CvPoint2D32f inCoord); 

	// To convert the photo coordinates to the pixel coordinates.
	void ConvertPhoto2PixCoords(int corner_count,	
								CvPoint2D32f * arrCoords, 
								/*out*/CvPoint2D32f * cornersB);

	// To convert the photo coordinates to the pixel coordinates.
	CvPoint2D32f ConvertPhoto2PixCoord(CvPoint2D32f inCoord);

	// To project points in the first image to the second image using the col.eq.
	void CollinearityEquation(CvPoint2D32f * arrCoords, int corner_count);
	void CollinearityEquation(CvPoint2D32f * arrCoords, int corner_count, double Z);
	void CollinearityEquation(	/*in*/CvPoint2D32f * inCoords, 
								/*out*/CvPoint2D32f * outCoords, 
								int corner_count, 
								double Z1, double Z2);

	// To compute the distance from ptB to the line (pt1,pt2)
	double Distance2EpiLine(CvPoint2D32f &pt1, CvPoint2D32f &pt2, CvPoint2D32f &ptB);

	// This function will be called whenever a new image is captured. 
	void GetNewImage(string image_name, EO eo);
	void GetNewImage(string image_name, EO eo, double za);

	// This function will perform the tracking. 
	void Track();

	// This function is to perform a 2D similarity transform.
	SIM2D_TRANS_PARAM Sim2DTransform(); // Pre-programmed coord. positions
	SIM2D_TRANS_PARAM Sim2DTransform(	int npoints,
										CvPoint2D32f * coords1, 
										CvPoint2D32f * coords2);
										// Generalize version

	// This function is to remove outliers from characteristics of optical flow.
	void OptFlowOutlierRemoval(/*in,out*/int &corner_count, 
								/*in,out*/CvPoint2D32f * cornersA, 
								/*in,out*/CvPoint2D32f * cornersB,
								/*in,out*/int * cornersA_GPID);

	// This function is to remove outliers by epipolar line.
	void AffineOutlierRemoval(	/*in,out*/int &corner_count, 
							/*in,out*/CvPoint2D32f * cornersA, 
							/*in,out*/CvPoint2D32f * cornersB,
							/*in,out*/int * cornersA_GPID);

	// This function is to remove outliers by epipolar line.
	void EpiOutlierRemoval(	/*in,out*/int &corner_count, 
							/*in,out*/CvPoint2D32f * cornersA, 
							/*in,out*/CvPoint2D32f * cornersB,
							/*in,out*/int * cornersA_GPID);

	// Return the set of tie-points for the latest two images.
	VEC_TP GetTiePoints(int IMG_ID);
	int GetTiePoints(int IMG_ID, /*out*/VEC_IP * vec_tie_point, int last_stored_GP_ID = -1); 

	// Remove the obsoleted tie points from the m_vecTP vector.
	int RemoveObsoletedTiepoints(int iObsoletedIMG);

	// Variables declaration.

	string	m_CFG_File;			// Configuration file name.
	string	m_IO_File;			// IO file name.
	string	m_img_directory;	// Directory that stores the images.
	
	DEQ_EO	m_deqEO;			// Deque of EO parameters for images in each tracking round.
	DEQ_STRING m_deqIMG;		// Deque of image file names in each tracking round.
	DEQ_DOUBLE m_deqZA;			// Deque of ZA values.
								// These deques contains maximum of two EO/image names at a time.

	VEC_TP	m_vecTP;			// Vector of tie points that are the result of tracking.
	int *	m_npoints_pattern_block;	// To store the number of tie points in each pattern block.

	IO		m_IO;				// IO / Principle point
	IS		m_IS;				// IO / Image dimension
	IS		m_scaledIS;			// IO / Image dimension after scaling.
	PS		m_PS;				// IO / Pixel size
	RDC		m_RDC;				// IO / Radial distortion
	PS		m_ICenter;			// Image center

	CvMat	*R1;				// 3D rotation matrix for the first image
	CvMat	*R2;				// 3D rotation matrix for the second image

	int		SCALE_FACTOR;		// Image dimension scale factor.
								// 1: Original, 2: Half size, 4: Quater size, ...
	int		ROT_ORDER;			// Rotational matrix order
								// 1: Rx*Ry*Rz or 2: Rz*Ry*Rx
	int		ROT_TYPE;			// Rotational matrix type: 1 or 2
								// Type 1 
								//			Rx = [1 0 0; 0 cos(om) -sin(om); 0 sin(om) cos(om)]
								//			Ry = [cos(ph) 0 sin(ph); 0 1 0; -sin(ph) 0 cos(ph)]
								//			Rz = [cos(kp) -sin(kp) 0; sin(kp) cos(kp) 0; 0 0 1]
								// Type 2
								//			Rx = [1 0 0; 0 cos(om) sin(om); 0 -sin(om) cos(om)]
								//			Ry = [cos(ph) 0 -sin(ph); 0 1 0; sin(ph) 0 cos(ph)]
								//			Rz = [cos(kp) sin(kp) 0; -sin(kp) cos(kp) 0; 0 0 1]
	bool	m_bOutputTP;		// Whether or not to output the tie points.

	bool	m_bReloadIMG;		// To re-read both image without reusing pre-load image.

	int		H_SPACE;			// The number of pixels for each pattern block in vertical-basis.
	int		W_SPACE;			// The number of pixels for each pattern block in horizontal-basis.
	int		RUNNING_POINT_ID;	// The running number for Ground Point ID.
	int		FAST_MIN_DISTANCE;	// Min distance required for each feature

	// Configuration parameters used by the KLT tracker.

	int		MAX_CORNERS;			// Limit the number of points to track per block
	int		MAX_BLOCK_TP;			// Limit the max. number of TP per block, 0 for ALL
	double	MINIMUM_EIGENVALUE;		// [KLT]Minimum eigen value
	double	MIN_DISTANCE;			// [KLT/FAST]Minimum distance between adjacent good features point
	DETECTOR_TYPE EXTRACTOR;		// [KLT]Feature extracker type: 1 for KLT, 2 for FAST
	TRACKING_MODEL MOTION_MODEL;	// [KLT]Tracking model
	int		BLOCK_AUTOCORR;			// [KLT]Window size of autocorrelation
	int		WIN_SIZE_SUBPIXEL;		// [KLT]Window size for subpixel
	int		WIN_SIZE_TRACK;			// [KLT]Window size for tracking
	int		WIN_SIZE_CCC;			// [KLT]Window size for CCC
	int		SUBPATTERN;				// Pattern of the overlapping P x P e.g. 5x5, 4x4, 3x3
	int		CRI_NUM_ITERATION;		// [KLT]Stop Criteria - iteration
	double	CRI_EPSILON_LIMIT;		// [KLT]Stop Criteria - distance
	int		DEPTH_LEVELS;			// [KLT]Pyramid depth level
	int		FEATURE_ERROR;			// [KLT]Tracking error
	int		FAST_THRESHOLD;			// [FAST]Threshold for FAST feature detector
	bool	FAST_B_NON_MAX_SUPRESS;	// [FAST]1 for Non-maximum suppression, 0 otherwise
	int		FAST_NUM_FAST_PIXEL;	// [FAST]Number of FAST pixel between 9-12
	double	ZA;						// Average terrain elevation
	int		CCC_CHECK;				// Whether or not CCC double checking is performed.
	int		CAL_INITIAL_GUESS;		// Whether or not calculating initial guessed TP is performed.
	int		EPIPOLAR_CHECK;			// Whether or not performing outlier removal by epipolar line.
	double	EPIPOLAR_RATIO;			// Upperbound to reject point with distance > mean+RATIO*std
	double	CCC_THRESHOLD;			// CCC accepting threshold
	double	MIN_ROTATED_ANGLE;		// Minimum angle to rotate image (deg)
	int		KLT_TRACKING_FLAG;		// KLT Tracking flag
	int		COUNT_OUT_SCOPE_TP;		// Whether or not considering out-of-scope tie-point.
	int		AFFINE_REMOVAL;			// Whether or not performing outlier removal by affine model..
	double	AFFINE_RATIO;			// Upperbound to reject point with distance > mean+RATIO*std
	int		OPT_REMOVAL;			// Whether or not performing outlier removal by optical flow.
	double	OPT_RATIO;				// Upperbound to reject point with distance > mean+RATIO*std


	// KLT Tracker variables.
	// The runtime used variables are declared under the class scope to reduce the initialization overhead.
	IplImage* imgA;
	IplImage* imgB;

};


#endif // _KLTTRACKER_H_
