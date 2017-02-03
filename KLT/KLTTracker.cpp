/********************************************************
 * Project:			Automated AT
 * Last updated:	11 November 2012
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "KLTTracker.h"
#include "Global.h"
#include <iostream>
#include <fstream>		// To utilize ifstream for dealing with text file input.
#include <strstream>
#include "stdio.h"
#include <CRTDBG.H>		// To use _ASSERT
#include "Definition.h"

/********************************************************************
 * Description:	Construction
 ********************************************************************/
CKLTTracker::CKLTTracker(const string &_CFG_File)
{
	// Store the configuration file.
	m_CFG_File = _CFG_File; 

	// Initialize the GP ID.
	RUNNING_POINT_ID = 1;

	// Initialize the variables.
	imgA = NULL;	// Set image ptr to NULL
	imgB = NULL;	// Set image ptr to NULL
	R1 = NULL;		// NULL ptr to rotation matrix
	R2 = NULL;		// NULL ptr to rotation matrix
}

/********************************************************************
 * Description:	Destructor
 ********************************************************************/
CKLTTracker::~CKLTTracker()
{
	if (m_npoints_pattern_block)
		delete [] m_npoints_pattern_block;

	if (imgA != NULL)
		cvReleaseImage(&imgA);

	if (imgB != NULL)
		cvReleaseImage(&imgB);

}

/********************************************************************
 * Description:	Set member variables
 ********************************************************************/
void CKLTTracker::SetMembers(int n_init_num_images)
{

}

/********************************************************************
 * Description:	Wrapper function to get each line from the input file.
 ********************************************************************/
void CKLTTracker::GetLine(istream &in, char *buf, const int bufSize) const
{
	while(in.good())
	{
		in.getline(buf, bufSize);
		if(buf[0]!='%') break;
	}
}

/*******************************************************************
 * Description:	Read the KLT configuration file
 ********************************************************************/
void CKLTTracker::ReadConfig()
{
	// Open the configuration file to read.
	ifstream fileInput(m_CFG_File.c_str());

	if (!fileInput)
	{
		cout << "Error! Failed to open the configuration file." << endl;
		exit(0);
	}
	
	const int BUFFER_SIZE = 200;
	char buffer[BUFFER_SIZE+1] = {0};

	// Read the IO file name.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	istrstream(buffer) >> m_IO_File;

	// Read the image directory.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	istrstream(buffer) >> m_img_directory;

	// Read the scaling factor.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &SCALE_FACTOR);

	// Read the rotational matrix order.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &ROT_ORDER);

	// Read the rotational matrix type.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &ROT_TYPE);

	// Read MAX_CORNERS: the number of points to track.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &MAX_CORNERS);

	// Read MAX_BLOCK_TP: the maximum number of tie points for each block.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &MAX_BLOCK_TP);

	// Read TRACKER: type of feature detector either 1) KLT or 2) FAST.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &EXTRACTOR);

	// Read MOTION_MODEL: type of KLT tracking model {1,2,3,4}.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &MOTION_MODEL);

	// Read MINIMUM_EIGENVALUE: Minimum eigen value.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%lf", &MINIMUM_EIGENVALUE);

	// Read MIN_DISTANCE: Minimum distance between adjacent good features point.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%lf", &MIN_DISTANCE);
	
	// Read BLOCK_AUTOCORR: Window size of autocorrelation.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &BLOCK_AUTOCORR);

	// Read WIN_SIZE_SUBPIXEL: Window size for subpixel.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &WIN_SIZE_SUBPIXEL);

	// Read WIN_SIZE_TRACK: Window size for tracking.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &WIN_SIZE_TRACK);

	// Read WIN_SIZE_CCC: Window size for CCC.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &WIN_SIZE_CCC);

	// Read SUBPATTERN: Pattern of the overlapping P x P e.g. 5x5, 4x4, 3x3.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &SUBPATTERN);

	// Read CRI_NUM_ITERATION: Stop Criteria - iteration.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &CRI_NUM_ITERATION);

	// Read CRI_EPSILON_LIMIT: Stop Criteria - distance.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%lf", &CRI_EPSILON_LIMIT);

	// Read DEPTH_LEVELS: Pyramid depth level.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &DEPTH_LEVELS);

	// Read FEATURE_ERROR: Tracking error.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &FEATURE_ERROR);

	// Read FAST_THRESHOLD: Threshold for FAST feature extractor.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &FAST_THRESHOLD);

	// Read NON_MAX_SUPPRESSION: 1 for non-maximum suppression, 0 otherwise.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &FAST_B_NON_MAX_SUPRESS);

	// Read NUM_FAST_PIXEL: between 9-12.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &FAST_NUM_FAST_PIXEL);

	// Read minimum angle to rotate image (degree).
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%lf", &MIN_ROTATED_ANGLE);

	// Read if speeding up KLT by initial guessed TP will be performed.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &CAL_INITIAL_GUESS);

	// Read if considering out-of-scope tie-point.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &COUNT_OUT_SCOPE_TP);

	// Read if double check by CCC will be performed.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &CCC_CHECK);

	// Read the CCC threshold for double check if CCC_CHECK flag is set.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%lf", &CCC_THRESHOLD);

	// Read if outlier removal by epipolar line will be performed.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &EPIPOLAR_CHECK);
	
	// Read ratio of std of distance to epipolar line 
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%lf", &EPIPOLAR_RATIO);

	// Read if outlier removal by affine model will be performed.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &AFFINE_REMOVAL);
	
	// Read ratio of std of distance to affine model 
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%lf", &AFFINE_RATIO);

	// Read ZA: Average terrain elevation.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%lf", &ZA);

	// Read the flag whether to output the tie points.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d", &m_bOutputTP);

	// Close the file after finish reading.
	fileInput.close();
}

/********************************************************************
 * Description:	Read the IO file
 * NOTE:		NOT USED
 ********************************************************************/
void CKLTTracker::ReadIOFile()
{
	// Open the IO file to read.
	ifstream fIO(m_IO_File.c_str());

	if (!fIO)
	{
		cout << "Error! Failed to open the IO file." << endl;
		exit(0);
	}
	
	const int BUFFER_SIZE = 200;
	char buffer[BUFFER_SIZE+1] = {0};

	// Read the principle point.
	GetLine(fIO, buffer, BUFFER_SIZE);
	sscanf(buffer, "%lf,%lf,%lf", &m_IO.dPPX, &m_IO.dPPY, &m_IO.dF);

	// Read the radial distortion.
	GetLine(fIO, buffer, BUFFER_SIZE);
	sscanf(buffer, "%lg,%lg,%lg", &m_RDC.dK1, &m_RDC.dK2, &m_RDC.dK3);
	
	// Read the resolution / image dimension
	GetLine(fIO, buffer, BUFFER_SIZE);
	sscanf(buffer, "%d,%d", &m_IS.iWidth, &m_IS.iHeight);

	// Read the pixel size (mm)
	GetLine(fIO, buffer, BUFFER_SIZE);
	sscanf(buffer, "%lf,%lf", &m_PS.dX, &m_PS.dY);

	// Close the file after complete reading.
	fIO.close();
}

/********************************************************************
 * Description:	It performs all preliminary tasks prior to run the
 *				KLT tracking.
 ********************************************************************/
void CKLTTracker::Startup()
{
	// Read the KLT configuration file.
	ReadConfig();

	// Read the IO file.
	ReadIOFile();

	// Initialize the constant values required for KLT tracker.
	InitializeConstants();

	// Initialize the member variables.
	InitializeMembers();
}

/********************************************************************
 * Description:	Initialize the constant values required for KLT tracker.
 ********************************************************************/
void CKLTTracker::InitializeConstants()
{
	//START----- PERFORMANCE IMPROVEMENT SOLUTION -----

	// Adjust the parameters / image dimension
	// CAUTION: Need to handle the case when the image size is ODD number.
	m_scaledIS.iWidth = m_IS.iWidth/SCALE_FACTOR;	
	m_scaledIS.iHeight = m_IS.iHeight/SCALE_FACTOR;

	//END----- PERFORMANCE IMPROVEMENT SOLUTION -----

	// The size of 'vertical' subregion in the pattern block. 
	H_SPACE = m_scaledIS.iHeight/SUBPATTERN;

	// The size of 'horizontal' subregion in the pattern block. 
	W_SPACE = m_scaledIS.iWidth/SUBPATTERN;

	// Obtain the image center
	m_ICenter.dX = ((double)m_IS.iWidth - 1.0)/2;
	m_ICenter.dY = ((double)m_IS.iHeight - 1.0)/2;

	// Assign flag to re-read both image
	m_bReloadIMG = true;		

	// Minimum distance between FAST features
	FAST_MIN_DISTANCE = (int)(MIN_DISTANCE*MIN_DISTANCE);

}

/********************************************************************
 * Description:	Initialize the member variables.
 ********************************************************************/
void CKLTTracker::InitializeMembers()
{
	// Initialize the pattern block to store the number of tie points in each block.
	m_npoints_pattern_block = new int[SUBPATTERN * SUBPATTERN];

	if (NULL == m_npoints_pattern_block)
	{
		printf("Error! Failed to allocate memory.");
		exit(0);
	}

	if (true == m_bOutputTP)
	{
		// Intialialize the TP file.
		//FILE * fTP = fopen("TP.txt", "w");
		//fclose(fTP);
	}
}

/********************************************************************
 * Description:	For an image, the whole area will be divided into sub-region
 *				following the defined pattern e.g. 3x3 or 4x4. 
 *				This function check if the passed-in coordinate falls
 *				into which region ID. 
 * Note:		The region ID is calculated as
 *				RegionID = iH * SUBPATTERN + iW
 ********************************************************************/
int CKLTTracker::FindSubRegionFromCoord(float y_coord, float x_coord)
{
	int iH, iW;

	// Calculate the row/column of the pattern block.
	iH = (int)(y_coord/H_SPACE);
	iW = (int)(x_coord/W_SPACE);

	iH = (iH >= SUBPATTERN)?SUBPATTERN-1:iH;
	iW = (iW >= SUBPATTERN)?SUBPATTERN-1:iW;

	// Calculate and return the subregion ID.
	return iH*SUBPATTERN + iW;
}

int CKLTTracker::FindSubRegionFromCoord(float y_coord, float x_coord, /*out*/int &iH, /*out*/int &iW)
{
	// Calculate the row/column of the pattern block.
	iH = (int)(y_coord/H_SPACE);
	iW = (int)(x_coord/W_SPACE);

	iH = (iH >= SUBPATTERN)?SUBPATTERN-1:iH;
	iW = (iW >= SUBPATTERN)?SUBPATTERN-1:iW;

	// Calculate and return the subregion ID.
	return iH*SUBPATTERN + iW;
}

/********************************************************************
 * Description:	The function returns the subregion ID in which the
 *				caller function passes in the row/column ID of  the
 *				pattern block.
 * Note:		The region ID is calculated as
 *				RegionID = iH * SUBPATTERN + iW
 ********************************************************************/
int CKLTTracker::FindSubRegionFromID(int iH, int iW)
{
	// Calculate and return the subregion ID.
	return iH*SUBPATTERN + iW;
}

/********************************************************************
 * Description:	The function returns the row/colmn ID of the passed-in
 *				region ID.
 * Note:		The region ID is calculated as
 *				RegionID = iH * SUBPATTERN + iW
 ********************************************************************/
void CKLTTracker::FindRowColFromRegionID(int regionID, /*out*/int &iH, /*out*/int &iW)
{
	iH = regionID/SUBPATTERN;
	iW = regionID%SUBPATTERN;
}

/********************************************************************
 * Description:	For an image, the whole area will be divided into sub-region
 *				following the defined pattern e.g. 3x3 or 4x4. 
 *				For the specified region ID, this function returns the 
 *				beginning x, beginning y, length x, and length y, in order. 
 ********************************************************************/
void CKLTTracker::DefineROI(int iH, int iW, 
							/*out*/ int &beginX, /*out*/ int &beginY,
							/*out*/ int &lengthX, /*out*/ int &lengthY)
{
	// Define the beginning coordinates.
	beginX = iW * W_SPACE;
	beginY = iH * H_SPACE;

	// Define the length of region in x-axis.
	if (iW < SUBPATTERN) 
		lengthX = W_SPACE;
	else
		lengthX = m_scaledIS.iWidth - beginX;

	// Define the length of region in y-axis.
	if (iH < SUBPATTERN)
		lengthY = H_SPACE;
	else
		lengthY = m_scaledIS.iHeight - beginY;

}

/********************************************************************
 * Description:	The function is to check if any block in the image
 *				pattern does not locate any tie point.
 *				If at least one block does not have tie points located,
 *				the function return false.
 ********************************************************************/
bool CKLTTracker::IsBlankRegionExisted()
{
	int i;

	for (i = 0; i < SUBPATTERN*SUBPATTERN; i++)
	{
		if (0 == m_npoints_pattern_block[i])
		{
			return true;
		}
	}
	
	return false;
}

bool CKLTTracker::IsBlankRegionExisted(int imageID, int * existing_features, int &n_existing_features)
{
	int i;
	VEC_TP::iterator iter;

	// Reset the number of tie points in the subregion (existing_features)
	memset(existing_features, 0, SUBPATTERN * SUBPATTERN * sizeof(int));

	n_existing_features = 0;

	// Iterate through m_vecTP to count the number of tie points in each subregion.
	for (iter = m_vecTP.begin(); iter != m_vecTP.end(); iter++)
	{
		// Increasing the number of tie points in the specified 'imageID'.
		if (imageID == (*iter).imageID)
		{
			existing_features[(*iter).regionID]++;
			n_existing_features++;
		}
	}
	
	// If the pattern is 1x1, always force the program to extract new features.
	if (1 == SUBPATTERN)
		return true;

	// If at least one subregion contains no tie point, return true.
	for (i = 0; i < SUBPATTERN * SUBPATTERN; i++)
	{
		if (0 == existing_features[i])
		{
			return true;
		}
	}
	
	return false;
}

/********************************************************************
 * Description:	The function is to store existing features into each block
 *				of sub-pattern according to its region ID.
 ********************************************************************/
void CKLTTracker::AssignExistingFeaturesToBlock(int imageID, 
												vector<vector<VEC_INT> > &vblock_features, 
												int * existing_features, 
												int &n_existing_features)
{
	//int i;
	int iH, iW;
	VEC_TP::iterator iter;

	// Reset the number of tie points in the subregion (existing_features)
	memset(existing_features, 0, SUBPATTERN * SUBPATTERN * sizeof(int));

	n_existing_features = 0;

	// Iterate through m_vecTP to store existing features into each subregion.
	for (iter = m_vecTP.begin(); iter != m_vecTP.end(); iter++)
	{
		// Increasing the number of tie points in the specified 'imageID'.
		if (imageID == (*iter).imageID)
		{
			// Obtain the row/col ID of the feature's region ID.
			FindRowColFromRegionID((*iter).regionID, /*out*/iH, /*out*/iW);

			// Put the point ID into the sub-pattern.
			(vblock_features[iH][iW]).push_back((*iter).pointID);
			
			// Count the number of features in each region.
			existing_features[(*iter).regionID]++;
			
			// Increment the number of features.
			n_existing_features++;
		}
	}

}

/********************************************************************
 * Description:	To check if the passed-in point already existed.
 ********************************************************************/
bool CKLTTracker::IsMatchingPointsExisted(int imageID, CvPoint2D32f point)
{
	VEC_TP::iterator iter;
	double dist;

	// Iterate through the m_vecTP to validate if the point existed.
	for (iter = m_vecTP.begin(); iter != m_vecTP.end(); iter++) 
	{
		if (imageID == (*iter).imageID)
		{
			CvPoint2D32f mpcoord = (*iter).coord;

			// Check if the point already existed by
			//	calculating the Euclidian distance between point.
			dist =	(point.y - mpcoord.y) * (point.y - mpcoord.y) + 
					(point.x - mpcoord.x) * (point.x - mpcoord.x);
			
			if DBL_LT(dist, MIN_DISTANCE)
				return true;

		}
	}

	return false;
}

/********************************************************************
 * Description:	To add the signficant features to the list of features
 *				to be tracked.
 ********************************************************************/
void CKLTTracker::AddMatchedPointsForTracking(	int imageID, 
												CvPoint2D32f * cornersA, 
												int &nPoint, 
												int * cornersA_GPID)
{
	VEC_TP::iterator iter;

	// Reset the number of points to zero.
	nPoint = 0;

	// Continue searching tie points in the specified 'imageID'.
	for (iter = m_vecTP.begin(); iter != m_vecTP.end(); iter++)
	{
		if (imageID == (*iter).imageID)
		{
			CvPoint2D32f mpcoord = (*iter).coord;
			
			// Assign the tie point to the vector.
			cornersA[nPoint] = mpcoord;

			// Store the 'GPID' of the tie points in cornersA into cornersA_GPID
			cornersA_GPID[nPoint] = (*iter).pointID;

			nPoint++;			
		}
	}
}

/********************************************************************
 * Description:	This function is to calculate the initial guessed 
 *				positions of the tie points based on the collinearity
 *				equation.
 ********************************************************************/
void CKLTTracker::CalInitialGuesses(int corner_count, 
									CvPoint2D32f * cornersA, 
									/*out*/CvPoint2D32f * cornersB)
{
	// Initialize a temporaly buffer to store calculated coordinates.
	CvPoint2D32f * arrCoords;		
	arrCoords = new CvPoint2D32f[corner_count];

	// Convert the pixel coordinates to the photo coordinates.
	ConvertPix2PhotoCoords(corner_count, cornersA, arrCoords);

	// Project points in the first image to the second image using collinearity equation.
	CollinearityEquation(arrCoords, corner_count);

	// Convert the photo coordinates to the pixel coordinates.
	ConvertPhoto2PixCoords(corner_count, arrCoords, cornersB);

	// Free the memory.
	delete [] arrCoords;
}

void CKLTTracker::CalInitialGuesses(int corner_count, 
									CvPoint2D32f * cornersA, 
									/*out*/CvPoint2D32f * cornersB,
									double rot_a, double rot_b, 
									double rot_tx, double rot_ty)
{
	// Initialize a temporaly buffer to store calculated coordinates.
	CvPoint2D32f * arrCoords;		
	arrCoords = new CvPoint2D32f[corner_count];

	// Convert the pixel coordinates to the photo coordinates.
	ConvertPix2PhotoCoords(corner_count, cornersA, arrCoords);

	// Project points in the first image to the second image using collinearity equation.
	CollinearityEquation(arrCoords, corner_count);

	// Convert the photo coordinates to the pixel coordinates.
	ConvertPhoto2PixCoords(corner_count, arrCoords, cornersB);

	// Transform to rotated image coordinate
	for (int i = 0; i < corner_count; i++)
	{
		double x = cornersB[i].x;
		double y = cornersB[i].y;
		cornersB[i].x = (float)(rot_a*x + rot_b*y + rot_tx);
		cornersB[i].y = (float)(-rot_b*x + rot_a*y + rot_ty);
	}

	// Free the memory.
	delete [] arrCoords;
}

/********************************************************************
 * Description:	This function is to compute the epipolar line and  
 *				classify points which is far from the line as outliers.
 ********************************************************************/
void CKLTTracker::EpiOutlierRemoval(/*in,out*/int &corner_count, 
									/*in,out*/CvPoint2D32f * cornersA, 
									/*in,out*/CvPoint2D32f * cornersB,
									/*in,out*/int * cornersA_GPID)
{
	double deltaZ = 100.0;		// Pre-defined the terrain delta.
	double Z1 = m_deqZA[IMG_A] - deltaZ;	// The first terrain elevation.
	double Z2 = m_deqZA[IMG_A] + deltaZ;	// The second terrain elevation.

	double dDist;
	double dSum, dSumSQ, dMean, dStd;

	CvPoint2D32f * cornersApxl;	// Pixel coords of cornersA.
	CvPoint2D32f * epipoint1;	// Corresponding points computed with Z1.
	CvPoint2D32f * epipoint2;	// Corresponding points computed with Z2.
	double * arrDistance;		// Array storing distances from point to epi. line.

	cornersApxl = new CvPoint2D32f[corner_count];
	epipoint1 = new CvPoint2D32f[corner_count];
	epipoint2 = new CvPoint2D32f[corner_count];
	arrDistance = new double[corner_count];

	// Compute a set of corresponding points for both elevations.
	{
		// Convert the pixel coord. to the photo coord.
		ConvertPix2PhotoCoords(corner_count, cornersA, cornersApxl);

		// Compute corresponding points in the second image with Z1.
		CollinearityEquation(	/*in*/cornersApxl, 
								/*out*/epipoint1, 
								corner_count, 
								Z1, Z1);
		ConvertPhoto2PixCoords(corner_count, epipoint1, epipoint1);

		// Compute corresponding points in the second image with Z2.
		CollinearityEquation(	/*in*/cornersApxl, 
								/*out*/epipoint2, 
								corner_count, 
								Z2, Z2);
		ConvertPhoto2PixCoords(corner_count, epipoint2, epipoint2);

	}

	dSum	= 0.0;
	dSumSQ	= 0.0;

	// Compute the distance from point to epipolar line.
	for (int i = 0; i < corner_count; i++)
	{
		dDist = Distance2EpiLine(epipoint1[i], epipoint2[i], cornersB[i]);
		arrDistance[i] = dDist;

		// Store the data for later computation
		dSum += dDist;
		dSumSQ += dDist*dDist;
	}

	// Compute mean and standard deviation.
	dStd = 0.0;
	dMean = dSum/corner_count;
	if (corner_count > 1)
		dStd = sqrt((corner_count*dSumSQ - (dSum*dSum))/corner_count/(corner_count-1));

	// === HOW TO DISCARD MISMATCH POINTS ====//

	// Set the upperbound for the reject threshold
	double dUpperThreshold = dMean + EPIPOLAR_RATIO*dStd;

			printf("\t\tAverage distance to epipolar line %lf\n", dMean);
			printf("\t\tStandard deviation of distance to epipolar line %lf\n", dStd);
			printf("\t\tEpipolar reject threshold %lf\n", dUpperThreshold);

	// Discard points likely to be mismatched
	int iNumValidPoint = -1;
	for (int i = 0; i < corner_count; i++)
	{
			// DEBUG
			printf("\t\tCornerB[%d,GP%d] (%lf,%lf) is far from epipolar line by %lf ", 
						i, cornersA_GPID[i], cornersB[i].x, cornersB[i].y, arrDistance[i]);
			// DEBUG


		// Reject the point if it is greater than the upperbound.
		if DBL_GT(arrDistance[i],dUpperThreshold)
		{
			continue;
		}

		// Increment the index of valid point.
		iNumValidPoint = iNumValidPoint + 1; 

		// If the current index is not the last valid point
		// (which means there is a hole before this current point),
		// re-arrange the cornersA and cornersB array by move this
		// current point up.
		if (i != iNumValidPoint)
		{
			// Move the current point the lastest valid position.
			cornersA[iNumValidPoint].x = cornersA[i].x;
			cornersA[iNumValidPoint].y = cornersA[i].y;
			cornersB[iNumValidPoint].x = cornersB[i].x;
			cornersB[iNumValidPoint].y = cornersB[i].y;
			cornersA_GPID[iNumValidPoint] = cornersA_GPID[i];
		}
	}

	// Update the number of matched points.
	corner_count = iNumValidPoint+1;
}

/********************************************************************
 * Description:	This function is to compute the affine transformation and  
 *				classify points which is far from the line as outliers.
 ********************************************************************/
void CKLTTracker::AffineOutlierRemoval(/*in,out*/int &corner_count, 
									/*in,out*/CvPoint2D32f * cornersA, 
									/*in,out*/CvPoint2D32f * cornersB,
									/*in,out*/int * cornersA_GPID)
{
	double dDist;
	double dSum, dSumSQ, dMean, dStd;
	int N_PARAMS;	// Number of unknown parameters
	int i;
	double a0,a1,a2,a3,a4,a5,a6,a7;

	CvMat * Y;	// The Y matrix
	CvMat * A;	// Design matrix
	CvMat * X;	// Solution

	// To compute the solution
	CvMat * At;		// Transpose of A
	CvMat * tmp1;
	CvMat * tmp2;

	N_PARAMS = 6; // For affine transformation


	double * arrDistance;		// Array storing distances from point to epi. line.
	arrDistance = new double[corner_count];

	// Construct the matrices.
	Y = cvCreateMat(corner_count*2, 1, CV_32FC1);
	A = cvCreateMat(corner_count*2, N_PARAMS, CV_32FC1);
	X = cvCreateMat(N_PARAMS, 1, CV_32FC1);

	// Set up the matrix
	for(i=0; i<corner_count; i++ ) 
	{
		// Add 1 feature from CornerB to the Y matrix 
		cvmSet(Y, i*2, 0, cornersB[i].x);	
		cvmSet(Y, i*2+1, 0, cornersB[i].y);	

		// AFFINE_TRANS
		{
			// Add 1 feature from CornerA to the A matrix 
			cvmSet(A, i*2, 0, cornersA[i].x);	
			cvmSet(A, i*2, 1, cornersA[i].y);	
			cvmSet(A, i*2, 2, 1);	
			cvmSet(A, i*2, 3, 0);	
			cvmSet(A, i*2, 4, 0);	
			cvmSet(A, i*2, 5, 0);	
			cvmSet(A, i*2+1, 0, 0);	
			cvmSet(A, i*2+1, 1, 0);	
			cvmSet(A, i*2+1, 2, 0);	
			cvmSet(A, i*2+1, 3, cornersA[i].x);	
			cvmSet(A, i*2+1, 4, cornersA[i].y);	
			cvmSet(A, i*2+1, 5, 1);	
		}
	}

	// Construct the matrices.
	At = cvCreateMat(N_PARAMS, corner_count*2, CV_32FC1);	// A transpose
	tmp1 = cvCreateMat(N_PARAMS, N_PARAMS, CV_32FC1);		// At*A, inv(At*A)
	tmp2 = cvCreateMat(N_PARAMS, 1, CV_32FC1);				// At*y

	// Compute the solution
	// X = inv(At * A) * (At * Y)
	cvTranspose(A, At);
	cvMatMul(At,A,tmp1);
	cvInvert(tmp1, tmp1);
	cvMatMul(At,Y,tmp2);
	cvMatMul(tmp1,tmp2,X);

	// Motion transformation parameters
	a0 = cvmGet(X, 0, 0);		a1 = cvmGet(X, 1, 0);
	a2 = cvmGet(X, 2, 0);		a3 = cvmGet(X, 3, 0);
	a4 = cvmGet(X, 4, 0);		a5 = cvmGet(X, 5, 0);

		printf("\t\t Affine: %lf %lf %lf\n", a0, a1, a2);
		printf("\t\t         %lf %lf %lf\n", a3, a4, a5);

	//****************************************
	// Partition the feature into Fin or Fout
	//****************************************

	double TransX1, TransY1;

	dSum	= 0.0;
	dSumSQ	= 0.0;

	for(i=0; i<corner_count; i++ ) 
	{
		// AFFINE_TRANS
		{
			TransX1 =	a0 * cornersA[i].x +
						a1 * cornersA[i].y +				
						a2;

			TransY1 =	a3 * cornersA[i].x +
						a4 * cornersA[i].y +				
						a5;
		}
		
		// Compute the distance to the transformation matrix
		double distX = cornersB[i].x - TransX1;
		double distY = cornersB[i].y - TransY1;
		dDist = sqrt(distX*distX + distY*distY);
		// Store the data for later computation
		arrDistance[i] = dDist;
		dSum += dDist;
		dSumSQ += dDist*dDist;
	}

	// Compute mean and standard deviation.
	dStd = 0.0;
	dMean = dSum/corner_count;
	if (corner_count > 1)
		dStd = sqrt((corner_count*dSumSQ - (dSum*dSum))/corner_count/(corner_count-1));

	// Set the upperbound for the reject threshold
	double dUpperThreshold = dMean + AFFINE_RATIO*dStd;

			printf("\t\tAverage distance to affine model %lf\n", dMean);
			printf("\t\tStandard deviation of distance to affine model %lf\n", dStd);
			printf("\t\tAffine model reject threshold %lf\n", dUpperThreshold);

	// Discard points likely to be mismatched
	int iNumValidPoint = -1; //Index of the last valid pair of points.
	for (int i = 0; i < corner_count; i++)
	{
		// Reject the point if it is greater than the upperbound.
		if DBL_GT(arrDistance[i],dUpperThreshold)
		{
				// DEBUG
			printf("REJECTED A[ID%d,%lf,%lf] and B[%lf,%lf] by %lf\n",
				cornersA_GPID[i],cornersA[i].x, cornersA[i].y,
				cornersB[i].x, cornersB[i].y,arrDistance[i]);
				// DEBUG

			continue;
		}

		// Increment the index of valid point.
		iNumValidPoint = iNumValidPoint + 1; 

		// If the current index is not the last valid point
		// (which means there is a hole before this current point),
		// re-arrange the cornersA and cornersB array by move this
		// current point up.
		if (i != iNumValidPoint)
		{
			// Move the current point the lastest valid position.
			cornersA[iNumValidPoint].x = cornersA[i].x;
			cornersA[iNumValidPoint].y = cornersA[i].y;
			cornersB[iNumValidPoint].x = cornersB[i].x;
			cornersB[iNumValidPoint].y = cornersB[i].y;
			cornersA_GPID[iNumValidPoint] = cornersA_GPID[i];
		}
	}

	// Update the number of matched points.
	corner_count = iNumValidPoint+1;

	delete [] arrDistance;
}

/********************************************************************
 * Description:	This function is to compute the magnitude and  
 *				direction of optical flows and classify points 
 *				that are different as outliers.
 ********************************************************************/
void CKLTTracker::OptFlowOutlierRemoval(/*in,out*/int &corner_count, 
									/*in,out*/CvPoint2D32f * cornersA, 
									/*in,out*/CvPoint2D32f * cornersB,
									/*in,out*/int * cornersA_GPID)
{
	int i;
	// To compute optical flow magnitude and direction
	double dMag, dSumMag, dSumSQMag, dMeanMag, dStdMag;
	double dDir, dSumDir, dSumSQDir, dMeanDir, dStdDir;
	double dX, dY;
	double dDummy;

	double * arrMag;		// Array storing magnitude.
	arrMag = new double[corner_count];
	double * arrDir;		// Array storing direction.
	arrDir = new double[corner_count];

	dSumMag	= 0.0;
	dSumSQMag	= 0.0;
	dSumDir	= 0.0;
	dSumSQDir	= 0.0;

	for(i=0; i<corner_count; i++ ) 
	{
		dX = cornersA[i].x - cornersB[i].x;
		dY = cornersA[i].y - cornersB[i].y;
		dDummy = dX*dX + dY*dY;
		dSumSQMag += dDummy; // Store data here to avoid double operation
		dMag = sqrt(dDummy);	// Magnitude of the optical flow vector
		dDir = rad2deg(atan2(dY,dX)); // Direction of the optical flow vector

		// Store the data for later computation
		arrMag[i] = dMag;
		dSumMag += dMag;
		arrDir[i] = dDir;
		dSumDir += dDir;
		dSumSQDir += dDir*dDir;

	}

	// Compute mean and standard deviation for magnitude
	dStdMag = 0.0;
	dMeanMag = dSumMag/corner_count;
	if (corner_count > 1)
		dStdMag = sqrt((corner_count*dSumSQMag - (dSumMag*dSumMag))/corner_count/(corner_count-1));

	// Compute mean and standard deviation for direction
	dStdDir = 0.0;
	dMeanDir = dSumDir/corner_count;
	if (corner_count > 1)
		dStdDir = sqrt((corner_count*dSumSQDir - (dSumDir*dSumDir))/corner_count/(corner_count-1));


	// Set the upperbound for the reject threshold
	double dUpperThresholdMag = dMeanMag + OPT_RATIO*dStdMag;
	double dLowerThresholdMag = dMeanMag - OPT_RATIO*dStdMag;

			printf("\t\tAverage magnitude optical flow %lf\n", dMeanMag);
			printf("\t\tStandard deviation of magnitude %lf\n", dStdMag);
			printf("\t\tMagnitude reject threshold %lf, %lf\n", 
						dUpperThresholdMag, dLowerThresholdMag);

	double dUpperThresholdDir = dMeanDir + OPT_RATIO*dStdDir;
	double dLowerThresholdDir = dMeanDir - OPT_RATIO*dStdDir;

			printf("\t\tAverage direction optical flow %lf\n", dMeanDir);
			printf("\t\tStandard deviation of direction %lf\n", dStdDir);
			printf("\t\tDirection reject threshold %lf, %lf\n", 
						dUpperThresholdDir, dLowerThresholdDir);

	// Discard points likely to be mismatched
	int iNumValidPoint = -1; //Index of the last valid pair of points.
	for (int i = 0; i < corner_count; i++)
	{
		// Reject the point if it is greater than the upperbound.
		if ( DBL_GT(arrMag[i],dUpperThresholdMag) || // Mag > upper threshold
			 DBL_LT(arrMag[i],dLowerThresholdMag) || // Mag < lower threshold
			 DBL_GT(arrDir[i],dUpperThresholdDir) || // Dir > upper threshold
			 DBL_LT(arrDir[i],dLowerThresholdDir) )  // Dir < lower threshold
		{
				// DEBUG
			printf("REJECTED A[ID%d,%lf,%lf] and B[%lf,%lf] by %lf,%lf deg\n",
				cornersA_GPID[i],cornersA[i].x, cornersA[i].y,
				cornersB[i].x, cornersB[i].y, arrMag[i], arrDir[i]);
				// DEBUG

			continue;
		}

		// Increment the index of valid point.
		iNumValidPoint = iNumValidPoint + 1; 

		// If the current index is not the last valid point
		// (which means there is a hole before this current point),
		// re-arrange the cornersA and cornersB array by move this
		// current point up.
		if (i != iNumValidPoint)
		{
			// Move the current point the lastest valid position.
			cornersA[iNumValidPoint].x = cornersA[i].x;
			cornersA[iNumValidPoint].y = cornersA[i].y;
			cornersB[iNumValidPoint].x = cornersB[i].x;
			cornersB[iNumValidPoint].y = cornersB[i].y;
			cornersA_GPID[iNumValidPoint] = cornersA_GPID[i];
		}
	}

	// Update the number of matched points.
	corner_count = iNumValidPoint+1;

	delete [] arrMag;
	delete [] arrDir;
}

/********************************************************************
 * Description:	This function is to compute the distance from the ptB
 *				in the 2nd image to its corresponding epipolar line
 *				computed from pt1 and pt2.
 * 
 * Line equation: y = mx+b --> mx-y+b = 0
 * where	m = (y2-y1)/(x2-x1)
 *			b = y1 - m*x1
 *
 * Distance to line: 
 *
 *			Dist =	|m(xB)-(yB)+ b|
 *					---------------
 *					sqrt(m*m + 1)
 *
 ********************************************************************/
double CKLTTracker::Distance2EpiLine(CvPoint2D32f &pt1, 
									 CvPoint2D32f &pt2, 
									 CvPoint2D32f &ptB)
{
	double m, b, dist;

	if ((int)pt2.x == (int)pt1.x )
	{
		// Vertical line: x = pt1.x
		m = 1.0;
		b = pt1.x;
	}
	else
	{
		m = (pt2.y - pt1.y)/(pt2.x - pt1.x);
		b = pt1.y - m*pt1.x;
	}
	
	dist = fabs( m*ptB.x - ptB.y + b)/sqrt(m*m + 1);

	return dist;
}

/********************************************************************
 * Description:	This function is to compute the 2D similarity transform
 *				parameters from a given set of image points.
 *
 * 2D similarity transform
 * -----------------------
 *
 *	|x21|	|x11	-y11	1	0|
 *	|y21|	|y11	x11		0	1|		|a11|
 *	|x22|	|x12	-y12	1	0|		|a12|
 *	|y22| =	|y12	x12		0	1|	X	|tx |
 *	 ...			  ...				|ty |
 *	|x2n|	|x1n	-y1n	1	0|
 *	|y2n|	|y1n	x1n		0	1|
 *
 *    y   =   Ax
 *	  x	  = inv(AtA)(Aty)
 *		
 * At least 2 pairs of corresponding points are required.
 *
 *		----------
 *		*		 |
 *		---------*
 * Coordinates: 1: (m_IS.iHeight/2,0) 
 *				2: (m_IS.iHeight,m_IS.iWidth)
 *
 ********************************************************************/
SIM2D_TRANS_PARAM CKLTTracker::Sim2DTransform()
{
	SIM2D_TRANS_PARAM Sim2DParam;	// 2D similarity transform parameters
	CvMat * A;						// Design matrix
	CvMat * y;						// The y matrix
	CvMat * x;						// Solution matrix
	CvMat * At;						// Transpose of A
	CvMat * tmp1;					// At*A, inv(At*A)
	CvMat * tmp2;					// At*y
	CvPoint2D32f* coords1;			// Pixel coordinates in Image#1
	CvPoint2D32f* coords2;			// Pixel coordinates in Image#2

	int npoints = 2;				// The number of points to be 

	// 1. Input the coordinate of points to compute 2D transform.
	coords1 = new CvPoint2D32f[npoints];
	coords2 = new CvPoint2D32f[npoints];
	// For coordinate: (m_IS.iHeight/2,0) 
	coords1[0].x = 0.0;	
	coords1[0].y = (float)(m_scaledIS.iHeight/2.0);
	// For coordinate: (m_IS.iHeight,m_IS.iWidth) 
	coords1[1].x = (float)m_scaledIS.iWidth;	
	coords1[1].y = (float)m_scaledIS.iHeight;

	// 2. Obtain the corresponding points in the
	//	  second image using the collinearity 
	//	  equation.

	// Convert the pixel coordinates to the photo coordinates.
	ConvertPix2PhotoCoords(npoints, coords1, coords2);

	// Project points in the first image to the second image using collinearity equation.
	CollinearityEquation(coords2, npoints);

	// Convert the photo coordinates to the pixel coordinates.
	ConvertPhoto2PixCoords(npoints, coords2, /*pixel coord*/coords2);

	// 3. Perform the 2D similarity transform
	Sim2DParam = Sim2DTransform(npoints, coords1, coords2);

	return Sim2DParam;
}

SIM2D_TRANS_PARAM CKLTTracker::Sim2DTransform(	int npoints,
												CvPoint2D32f * coords1, 
												CvPoint2D32f * coords2)
{
	SIM2D_TRANS_PARAM Sim2DParam;	// 2D similarity transform parameters
	CvMat * A;						// Design matrix
	CvMat * y;						// The y matrix
	CvMat * x;						// Solution matrix

	// Construct the matrices.
	A = cvCreateMat(npoints*2, 4, CV_32FC1);
	y = cvCreateMat(npoints*2, 1, CV_32FC1);
	
	// Input values to the matrix
	for (int i = 0; i < npoints; i++)
	{
		// A matrix - 2rows x 4cols for each point		
		// 	|x11	y11	1	0|
		cvmSet(A, i*2, 0, coords2[i].x); 
		cvmSet(A, i*2, 1, coords2[i].y); 
		cvmSet(A, i*2, 2, 1.0); 
		cvmSet(A, i*2, 3, 0.0); 
		//	|y11	-x11	0	1|		
		cvmSet(A, i*2+1, 0, coords2[i].y); 
		cvmSet(A, i*2+1, 1, -coords2[i].x); 
		cvmSet(A, i*2+1, 2, 0.0); 
		cvmSet(A, i*2+1, 3, 1.0); 
		
		// y matrix - 2rows x 1col for each point
		//	|x21; y21|
		cvmSet(y, i*2, 0, coords1[i].x); 
		cvmSet(y, i*2+1, 0, coords1[i].y);  
	}

	CvMat * At;		// Transpose of A
	CvMat * tmp1;
	CvMat * tmp2;

	// Construct the matrices.
	At = cvCreateMat(4, npoints*2, CV_32FC1);	// A transpose
	tmp1 = cvCreateMat(4, 4, CV_32FC1);			// At*A, inv(At*A)
	tmp2 = cvCreateMat(4, 1, CV_32FC1);			// At*y
	x = cvCreateMat(4, 1, CV_32FC1);			// inv(At*A)*(At*y)

	// Compute the solution
	// x = inv(At * A) * (At * y)
	cvTranspose(A, At);
	cvMatMul(At,A,tmp1);
	cvInvert(tmp1, tmp1);
	cvMatMul(At,y,tmp2);
	cvMatMul(tmp1,tmp2,x);

	float a = (float)cvmGet(x, 0, 0);
	float b = (float)cvmGet(x, 1, 0);

	// Compute scale
	Sim2DParam.scale = sqrt(a*a + b*b);

	// Compute rotation angle
	Sim2DParam.angle = (float)rad2deg((atan2(b,a)));

	// Translation vector
	Sim2DParam.tx = (float)cvmGet(x, 2, 0);
	Sim2DParam.ty = (float)cvmGet(x, 3, 0);

		// DEBUG
		printf("scale = %lf\n", Sim2DParam.scale);
		printf("angle = %lf\n", Sim2DParam.angle);
		printf("tx = %lf\n", Sim2DParam.tx);
		printf("ty = %lf\n", Sim2DParam.ty);
		// DEBUG

	// Release the memory
	cvReleaseMat(&A);
	cvReleaseMat(&y);
	cvReleaseMat(&x);
	cvReleaseMat(&At);
	cvReleaseMat(&tmp1);
	cvReleaseMat(&tmp2);

	return Sim2DParam;
}

/********************************************************************
 * Description:	To remove out-of-scope points computed after
 *				calculating initial guesses and rearrange
 *				relating arrays.
 ********************************************************************/
int CKLTTracker::RemoveOutOfScopePoints(	CvPoint2D32f * cornersA, 
											CvPoint2D32f * cornersB, 
											int * cornersA_GPID,
											int &nPoint)
{
	int i;
	int idxCurrent;	// Current index position to store the next benign point.
	int nOutPoints;	// Number of points that are out of scope.

	idxCurrent = 0;

	for (i = 0; i < nPoint; i++)
	{
		// Validate initial guessed tie-points.
		if ( DBL_LT(cornersB[i].y, 0.0) ||
			 DBL_LT(cornersB[i].x, 0.0) ||
			 DBL_GT(cornersB[i].y, double(m_scaledIS.iHeight+1)) ||
			 DBL_GT(cornersB[i].x, double(m_scaledIS.iWidth+1)) )
		{
			// Skip invalid point.
		}
		else
		{
			if (i != idxCurrent)
			{
				// Manage cornersA
				cornersA[idxCurrent].x = cornersA[i].x;
				cornersA[idxCurrent].y = cornersA[i].y;
			
				// Manage cornersB
				cornersB[idxCurrent].x = cornersB[i].x;
				cornersB[idxCurrent].y = cornersB[i].y;

				// Manage cornersA GPID
				cornersA_GPID[idxCurrent] = cornersA_GPID[i];
			}

			// Increment the current index.
			idxCurrent++;
		}

	}

	nOutPoints = nPoint - idxCurrent;
	
	nPoint = idxCurrent;

	return nOutPoints;
}

int CKLTTracker::ReplaceOutOfScopePoints(	CvPoint2D32f * cornersA, 
											CvPoint2D32f * cornersB, 
											int * cornersA_GPID,
											int &nPoint)
{
	int i;
	int idxCurrent;	// Current index position to store the next benign point.
	int nOutPoints;	// Number of points that are out of scope.

	idxCurrent = 0;

	for (i = 0; i < nPoint; i++)
	{
		// Validate initial guessed tie-points.
		if ( DBL_LT(cornersB[i].y, 0.0) ||
			 DBL_LT(cornersB[i].x, 0.0) ||
			 DBL_GT(cornersB[i].y, double(m_scaledIS.iHeight+1)) ||
			 DBL_GT(cornersB[i].x, double(m_scaledIS.iWidth+1)) )
		{
				cornersB[i].x = cornersA[i].x;
				cornersB[i].y = cornersA[i].y;
		}

	}

	nOutPoints = 0;

	return nOutPoints;
}

/********************************************************************
 * Description:	This function is to convert the pixel coordinates 
 *				stored in the cornerA to the photo coordinates
 *				and store them in the arrCoords variable.
 ********************************************************************/
void CKLTTracker::ConvertPix2PhotoCoords(	int corner_count,	
											CvPoint2D32f * cornersA, 
											/*out*/CvPoint2D32f * arrCoords)
{
	CvPoint2D32f coord;

	// Calculculate the photo-coordinates for each points (pixel based)
	for (int i = 0; i < corner_count; i++)
	{
		// Obtain the points in scaled coordinate.
		coord = cornersA[i];

		// Scale to original coordinate.
		coord.x = coord.x * SCALE_FACTOR;
		coord.y = coord.y * SCALE_FACTOR;		

		// Convert to photo-coordinate and store into the output array.
		arrCoords[i].x = (float)((coord.x - m_ICenter.dX) * m_PS.dX - m_IO.dPPX);
		arrCoords[i].y = (float)((m_ICenter.dY - coord.y) * m_PS.dY - m_IO.dPPY);
	}
}

CvPoint2D32f CKLTTracker::ConvertPix2PhotoCoord(CvPoint2D32f inCoord) 
{
	CvPoint2D32f outCoord;

	// Scale to original coordinate.
	inCoord.x = inCoord.x * SCALE_FACTOR;
	inCoord.y = inCoord.y * SCALE_FACTOR;

	// Calculculate the photo-coordinates for passed-in pixel coord.
	outCoord.x = (float)((inCoord.x - m_ICenter.dX) * m_PS.dX - m_IO.dPPX);
	outCoord.y = (float)((m_ICenter.dY - inCoord.y) * m_PS.dY - m_IO.dPPY);

	return outCoord;
}

/********************************************************************
 * Description:	This function is to convert the photo coordinates 
 *				stored in the arrCoords to the pixel coordinates
 *				and store them in the cornersB.
 ********************************************************************/
void CKLTTracker::ConvertPhoto2PixCoords(	int corner_count,	
											CvPoint2D32f * arrCoords, 
											/*out*/CvPoint2D32f * cornersB)
{
	CvPoint2D32f coord;

	// Calculculate the pixel coordinates from the photo coordinates.
	for (int i = 0; i < corner_count; i++)
	{
		// Obtain the photo coordinate.
		coord = arrCoords[i];

		// Convert to pixel coordinate.
		coord.x = (float)(m_ICenter.dX + (coord.x + m_IO.dPPX) / m_PS.dX);
		coord.y = (float)(m_ICenter.dY - (coord.y + m_IO.dPPY) / m_PS.dY);

		// Scale the coordinate and store into the output array.
		cornersB[i].x = coord.x / SCALE_FACTOR;
		cornersB[i].y = coord.y / SCALE_FACTOR;
	}
}

CvPoint2D32f CKLTTracker::ConvertPhoto2PixCoord(CvPoint2D32f inCoord) 
{
	CvPoint2D32f outCoord;

	// Calculculate the photo-coordinates for passed-in pixel coord.
	outCoord.x = (float)(m_ICenter.dX + (inCoord.x + m_IO.dPPX) / m_PS.dX);
	outCoord.y = (float)(m_ICenter.dY - (inCoord.y + m_IO.dPPY) / m_PS.dY);

	// Scale the pixel coordinate.
	outCoord.x = outCoord.x / SCALE_FACTOR;
	outCoord.y = outCoord.y / SCALE_FACTOR;

	return outCoord;
}

/********************************************************************
 * Description:	This function is to create 3D rotation matrix for
 *				both pair of images.
 ********************************************************************/
void CKLTTracker::Create3DRotMatrix()
{
	EO		eo;
	double	x;				// Angle		
	double	cos_x, sin_x;	// cos(Angle), sin(Angle)
	CvMat	*Rx, *Ry, *Rz;	// Rotational matrix


	// 1. Initialize the rotational matrix.
	Rx		= cvCreateMat(3, 3, CV_32FC1);
	Ry		= cvCreateMat(3 ,3, CV_32FC1);
	Rz		= cvCreateMat(3, 3, CV_32FC1);

	// 2. Manage R1 rotaion matrix

	// Not re-create R1 but point to the value of previous R2
	if ((m_bReloadIMG == true) && (R2 != NULL))
	{
		if (R1 == NULL)
		{ // This will never be entered.
			R1 = cvCreateMat(3 ,3, CV_32FC1);
		}

		// Copy the rotation matrix R2 to R1
		cvCopy(R2, R1);			
		
	}
	else // Re-create the rotation matrix R1
	{
		if (R1 == NULL)
			R1 = cvCreateMat(3 ,3, CV_32FC1);

		// Create 3D rotation matrix R1
		eo = m_deqEO[IMG_A];

		// Construct rotational matrix for Image#1
		if (ROT_TYPE == 1)
		{
			//		|	1		0		0		|
			// Rx =	|	0	  cos(Om) -sin(Om)	|
			//		|	0	  sin(Om) cos(Om)	|		
			
			x		=	eo.dOmega;
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			cvmSet(Rx, 0, 0, 1.0);
			cvmSet(Rx, 0, 1, 0.0);
			cvmSet(Rx, 0, 2, 0.0);
			cvmSet(Rx, 1, 0, 0.0);
			cvmSet(Rx, 1, 1, cos_x);
			cvmSet(Rx, 1, 2, -sin_x);
			cvmSet(Rx, 2, 0, 0.0);
			cvmSet(Rx, 2, 1, sin_x);
			cvmSet(Rx, 2, 2, cos_x);
			
			//		|	cos(Ph)		0		sin(Ph)	|
			// Ry =	|	  0			1		   0	|
			//		|	-sin(Ph)	0		cos(Ph)	|
			
			x		=	eo.dPhi;
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			cvmSet(Ry, 0, 0, cos_x);
			cvmSet(Ry, 0, 1, 0.0);
			cvmSet(Ry, 0, 2, sin_x);
			cvmSet(Ry, 1, 0, 0.0);
			cvmSet(Ry, 1, 1, 1.0);
			cvmSet(Ry, 1, 2, 0.0);
			cvmSet(Ry, 2, 0, -sin_x);
			cvmSet(Ry, 2, 1, 0.0);
			cvmSet(Ry, 2, 2, cos_x);

			//		|	cos(Ka)		-sin(Ka)	0	|
			// Rz =	|	sin(Ka)		cos(Ka)		0	|
			//		|	  0			   0		1	|
			
			x		=	eo.dKappa;	// Kappa
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			cvmSet(Rz, 0, 0, cos_x);
			cvmSet(Rz, 0, 1, -sin_x);
			cvmSet(Rz, 0, 2, 0.0);
			cvmSet(Rz, 1, 0, sin_x);
			cvmSet(Rz, 1, 1, cos_x);
			cvmSet(Rz, 1, 2, 0.0);
			cvmSet(Rz, 2, 0, 0.0);
			cvmSet(Rz, 2, 1, 0.0);
			cvmSet(Rz, 2, 2, 1.0);
		}
		else // ROT_TYPE = 2
		{
			//		|	1			0		0		|
			// Rx =	|	0		 cos(Om)	sin(Om)	|
			//		|	0		 -sin(Om)	cos(Om)	|		
			
			x		=	eo.dOmega;
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			cvmSet(Rx, 0, 0, 1.0);
			cvmSet(Rx, 0, 1, 0.0);
			cvmSet(Rx, 0, 2, 0.0);
			cvmSet(Rx, 1, 0, 0.0);
			cvmSet(Rx, 1, 1, cos_x);
			cvmSet(Rx, 1, 2, sin_x);
			cvmSet(Rx, 2, 0, 0.0);
			cvmSet(Rx, 2, 1, -sin_x);
			cvmSet(Rx, 2, 2, cos_x);
			
			//		|	cos(Ph)		0		-sin(Ph)|
			// Ry =	|	  0			1		   0	|
			//		|	sin(Ph)		0		cos(Ph)	|
			
			x		=	eo.dPhi;
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			cvmSet(Ry, 0, 0, cos_x);
			cvmSet(Ry, 0, 1, 0.0);
			cvmSet(Ry, 0, 2, -sin_x);
			cvmSet(Ry, 1, 0, 0.0);
			cvmSet(Ry, 1, 1, 1.0);
			cvmSet(Ry, 1, 2, 0.0);
			cvmSet(Ry, 2, 0, sin_x);
			cvmSet(Ry, 2, 1, 0.0);
			cvmSet(Ry, 2, 2, cos_x);

			//		|	cos(Ka)		sin(Ka)		0	|
			// Rz =	|	-sin(Ka)	cos(Ka)		0	|
			//		|	  0			   0		1	|
			
			x		=	eo.dKappa;	// Kappa
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			cvmSet(Rz, 0, 0, cos_x);
			cvmSet(Rz, 0, 1, sin_x);
			cvmSet(Rz, 0, 2, 0.0);
			cvmSet(Rz, 1, 0, -sin_x);
			cvmSet(Rz, 1, 1, cos_x);
			cvmSet(Rz, 1, 2, 0.0);
			cvmSet(Rz, 2, 0, 0.0);
			cvmSet(Rz, 2, 1, 0.0);
			cvmSet(Rz, 2, 2, 1.0);
		}

		if (ROT_ORDER == 1) // R = Rx * Ry * Rz
		{
			cvMatMul(Rx, Ry, R1);	// R = Rx * Ry
			cvMatMul(R1, Rz, R1);		// R = (Rx * Ry) * Rz
		}
		else // R = Rz * Ry * Rx
		{
			cvMatMul(Rz, Ry, R1);	// R = Rz * Ry
			cvMatMul(R1, Rx, R1);		// R = (Rz * Ry) * Rx
		}

	} // end if 

	// 2. Create the rotation matrix R2

	eo = m_deqEO[IMG_B];

	if (R2 == NULL)
		R2 = cvCreateMat(3 ,3, CV_32FC1);

	// Construct rotational matrix for Image#2
	{
		if (ROT_TYPE == 1)
		{
			//		|	1		0		0		|
			// Rx =	|	0	  cos(Om) -sin(Om)	|
			//		|	0	  sin(Om) cos(Om)	|		
			
			x		=	eo.dOmega;
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			cvmSet(Rx, 0, 0, 1.0);
			cvmSet(Rx, 0, 1, 0.0);
			cvmSet(Rx, 0, 2, 0.0);
			cvmSet(Rx, 1, 0, 0.0);
			cvmSet(Rx, 1, 1, cos_x);
			cvmSet(Rx, 1, 2, -sin_x);
			cvmSet(Rx, 2, 0, 0.0);
			cvmSet(Rx, 2, 1, sin_x);
			cvmSet(Rx, 2, 2, cos_x);
			
			//		|	cos(Ph)		0		sin(Ph)	|
			// Ry =	|	  0			1		   0	|
			//		|	-sin(Ph)	0		cos(Ph)	|
			
			x		=	eo.dPhi;
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			cvmSet(Ry, 0, 0, cos_x);
			cvmSet(Ry, 0, 1, 0.0);
			cvmSet(Ry, 0, 2, sin_x);
			cvmSet(Ry, 1, 0, 0.0);
			cvmSet(Ry, 1, 1, 1.0);
			cvmSet(Ry, 1, 2, 0.0);
			cvmSet(Ry, 2, 0, -sin_x);
			cvmSet(Ry, 2, 1, 0.0);
			cvmSet(Ry, 2, 2, cos_x);

			//		|	cos(Ka)		-sin(Ka)	0	|
			// Rz =	|	sin(Ka)		cos(Ka)		0	|
			//		|	  0			   0		1	|
			
			x		=	eo.dKappa;	// Kappa
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			cvmSet(Rz, 0, 0, cos_x);
			cvmSet(Rz, 0, 1, -sin_x);
			cvmSet(Rz, 0, 2, 0.0);
			cvmSet(Rz, 1, 0, sin_x);
			cvmSet(Rz, 1, 1, cos_x);
			cvmSet(Rz, 1, 2, 0.0);
			cvmSet(Rz, 2, 0, 0.0);
			cvmSet(Rz, 2, 1, 0.0);
			cvmSet(Rz, 2, 2, 1.0);
		}
		else // ROT_TYPE = 2
		{
			//		|	1			0		0		|
			// Rx =	|	0		 cos(Om)	sin(Om)	|
			//		|	0		 -sin(Om)	cos(Om)	|		
			
			x		=	eo.dOmega;
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			cvmSet(Rx, 0, 0, 1.0);
			cvmSet(Rx, 0, 1, 0.0);
			cvmSet(Rx, 0, 2, 0.0);
			cvmSet(Rx, 1, 0, 0.0);
			cvmSet(Rx, 1, 1, cos_x);
			cvmSet(Rx, 1, 2, sin_x);
			cvmSet(Rx, 2, 0, 0.0);
			cvmSet(Rx, 2, 1, -sin_x);
			cvmSet(Rx, 2, 2, cos_x);
			
			//		|	cos(Ph)		0		-sin(Ph)|
			// Ry =	|	  0			1		   0	|
			//		|	sin(Ph)		0		cos(Ph)	|
			
			x		=	eo.dPhi;
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			cvmSet(Ry, 0, 0, cos_x);
			cvmSet(Ry, 0, 1, 0.0);
			cvmSet(Ry, 0, 2, -sin_x);
			cvmSet(Ry, 1, 0, 0.0);
			cvmSet(Ry, 1, 1, 1.0);
			cvmSet(Ry, 1, 2, 0.0);
			cvmSet(Ry, 2, 0, sin_x);
			cvmSet(Ry, 2, 1, 0.0);
			cvmSet(Ry, 2, 2, cos_x);

			//		|	cos(Ka)		sin(Ka)		0	|
			// Rz =	|	-sin(Ka)	cos(Ka)		0	|
			//		|	  0			   0		1	|
			
			x		=	eo.dKappa;	// Kappa
			cos_x	=	cos(x);
			sin_x	=	sin(x);

			cvmSet(Rz, 0, 0, cos_x);
			cvmSet(Rz, 0, 1, sin_x);
			cvmSet(Rz, 0, 2, 0.0);
			cvmSet(Rz, 1, 0, -sin_x);
			cvmSet(Rz, 1, 1, cos_x);
			cvmSet(Rz, 1, 2, 0.0);
			cvmSet(Rz, 2, 0, 0.0);
			cvmSet(Rz, 2, 1, 0.0);
			cvmSet(Rz, 2, 2, 1.0);
		}

		if (ROT_ORDER == 1) // R = Rx * Ry * Rz
		{
			cvMatMul(Rx, Ry, R2);	// R = Rx * Ry
			cvMatMul(R2, Rz, R2);	// R = (Rx * Ry) * Rz
		}
		else // R = Rz * Ry * Rx
		{
			cvMatMul(Rz, Ry, R2);	// R = Rz * Ry
			cvMatMul(R2, Rx, R2);	// R = (Rz * Ry) * Rx
		}
	}

	cvReleaseMat(&Rx);
	cvReleaseMat(&Ry);
	cvReleaseMat(&Rz);

}

/********************************************************************
 * Description:	This function is to map one point in the first image
 *				to the corresponding point in the second image based
 *				on the collinearity equation.
 ********************************************************************/
void CKLTTracker::CollinearityEquation(/*in,out*/CvPoint2D32f * arrCoords, int corner_count)
{
	CollinearityEquation(arrCoords, arrCoords, corner_count, m_deqZA[IMG_A], m_deqZA[IMG_B]);
}

void CKLTTracker::CollinearityEquation(/*in,out*/CvPoint2D32f * arrCoords, int corner_count, double Z)
{
	CollinearityEquation(arrCoords, arrCoords, corner_count, Z, Z);
}

void CKLTTracker::CollinearityEquation(/*in*/CvPoint2D32f * inCoords, 
									   /*out*/CvPoint2D32f * outCoords, 
									   int corner_count, 
									   double ZA1,
									   double ZA2)
{
	double	r11, r12, r13;
	double	r21, r22, r23;
	double	r31, r32, r33;
	double	T1, T2, T3;
	double	dgroundX, dgroundY;
	CvPoint2D32f coord;
	int		i;
	EO		eo;

	// Create R1 and R2 matrices if they are not existent
	if ((R1 == NULL) || (R2 == NULL))
		Create3DRotMatrix();

	// STEP 1. Convert the pixel coordinates in Image#1 to the ground coordinates.
	//		|	r11	r12	r13	|
	// R =	|	r21	r22	r23	|
	//		|	r31	r32	r33	|

	r11 = cvmGet(R1, 0, 0);	r12 = cvmGet(R1, 0, 1);	r13 = cvmGet(R1, 0, 2);
	r21 = cvmGet(R1, 1, 0);	r22 = cvmGet(R1, 1, 1);	r23 = cvmGet(R1, 1, 2);
	r31 = cvmGet(R1, 2, 0);	r32 = cvmGet(R1, 2, 1);	r33 = cvmGet(R1, 2, 2);

	// Get the EO of Image#A
	eo = m_deqEO[IMG_A];

	// Compute the ground coordinates.

	for (i = 0; i < corner_count; i++)
	{
		coord = inCoords[i];

		T1 = r11*coord.x + r21*coord.y - r31*m_IO.dF;

		T2 = r12*coord.x + r22*coord.y - r32*m_IO.dF;

		T3 = r13*coord.x + r23*coord.y - r33*m_IO.dF;

		dgroundX = (ZA1 - eo.dZc)*T1/T3 + eo.dXc;
		dgroundY = (ZA1 - eo.dZc)*T2/T3 + eo.dYc;

		outCoords[i].x = (float)dgroundX;
		outCoords[i].y = (float)dgroundY;
	}

	// STEP 2. Convert the ground coordinates to the pixel coordinates in Image#2.

	r11 = cvmGet(R2, 0, 0);	r12 = cvmGet(R2, 0, 1);	r13 = cvmGet(R2, 0, 2);
	r21 = cvmGet(R2, 1, 0);	r22 = cvmGet(R2, 1, 1);	r23 = cvmGet(R2, 1, 2);
	r31 = cvmGet(R2, 2, 0);	r32 = cvmGet(R2, 2, 1);	r33 = cvmGet(R2, 2, 2);

	// Get the EO of Image#B
	eo = m_deqEO[IMG_B];

	// Compute the image coordinate.

	for (i = 0; i < corner_count; i++)
	{
		dgroundX = outCoords[i].x;
		dgroundY = outCoords[i].y;

		T1 = (dgroundX - eo.dXc)*r11 + (dgroundY - eo.dYc)*r12 + (ZA2 - eo.dZc)*r13;
		T2 = (dgroundX - eo.dXc)*r21 + (dgroundY - eo.dYc)*r22 + (ZA2 - eo.dZc)*r23;
		T3 = (dgroundX - eo.dXc)*r31 + (dgroundY - eo.dYc)*r32 + (ZA2 - eo.dZc)*r33;

		coord.x = (float)(- m_IO.dF*T1/T3);
		coord.y = (float)(- m_IO.dF*T2/T3);

		outCoords[i] = coord;
	}
}

void CKLTTracker::GetNewImage(string image_name, EO eo)
{
	// Pass the average ZA to the function
	GetNewImage(image_name, eo, ZA);
}

void CKLTTracker::GetNewImage(string image_name, EO eo, double za)
{
	// Store the new capture data to the deque which allows a maximum
	// of two sets of data at a time.	
	
	// Store the new capture data to the deque.
	m_deqEO.push_back(eo);
	m_deqIMG.push_back(image_name);
	m_deqZA.push_back(za);

	// Maintain the deque to store a maximum of two sets of data at a time.
	if (m_deqEO.size() > MAX_TRACKING_IMG)
	{
		// Obtain the obsoleted image ID.
		int iObsoletedIMG = m_deqEO[0].iImg;
	
		// Remove the obsoleted images.
		m_deqEO.pop_front();
		m_deqIMG.pop_front();
		m_deqZA.pop_front();

		// Remove the obsoleted tie points from the m_vecTP vector.
		(void) RemoveObsoletedTiepoints(iObsoletedIMG);

		// Reuse previous read image pointer.
		m_bReloadIMG = false;		
	}

	// Assert to ensure the deque stores a maximum of two images at a time.
	_ASSERT( m_deqEO.size() <= MAX_TRACKING_IMG );

	// If there exist two images, perform tracking.
	if (m_deqEO.size() == MAX_TRACKING_IMG)
	{
		Track();
	}

}

/********************************************************************
 * Description:	This function will perform the tracking. 
 ********************************************************************/
void CKLTTracker::Track()
{
	IplImage* pyrA = NULL;
	IplImage* pyrB = NULL;
	IplImage* eig_image = NULL;
	IplImage* tmp_image = NULL;
	IplImage* roiA = NULL;
	CvPoint2D32f* cornersA;
	CvPoint2D32f* cornersB;
	int * cornersA_GPID;			// To store the GPID of the corresponding tie points in cornersA.
	int corner_count;				// The number of corners (features)
	int NUM_UPPERBOUND;				// FAST::The upperbound number of features in subpattern
	int iH, iW, i, p;
	int r, c;						// Row and Column index of sub-pattern block.
	int xROI, yROI;					// Beginning coordinates of ROI (top-left coordinate)
	int widthROI, heightROI;		// Dimension of ROI
	int * existing_features;		// To store the number of tie points in each pattern block.
	int n_existing_features;		// To store the number of existing features.
	string imgA_file, imgB_file;	// String file name of Image A and Image B
	EO	eoA, eoB;					// EO parameters for Image A and Image B
	bool bRequireFindingGoodFeatures;	// A flag to indicate if the feature selection process will be run.
//	VEC_INT vec_num_matched_A(9, 0);	// To store the number of matched points for each block in imgA.
//	VEC_INT vec_num_matched_B(9, 0);	// To store the number of matched points for each block in imgB.
	VEC_INT vec_num_matched_A/*(16, 0)*/;	// To store the number of matched points for each block in imgA.
	VEC_INT vec_num_matched_B/*(16, 0)*/;	// To store the number of matched points for each block in imgB.
		//
		vec_num_matched_A.resize(SUBPATTERN*SUBPATTERN, 0);
		vec_num_matched_B.resize(SUBPATTERN*SUBPATTERN, 0);
		//
	IplImage *rotated_img = NULL;
	CvMat* rot_matrix = NULL;

	SIM2D_TRANS_PARAM sim2Dparam;	// 2D similarity transform parameters between two images.
	
	// Measure the computational time
	LARGE_INTEGER N_FREQUENCY;
	LARGE_INTEGER N_BEGIN_TIME;
	LARGE_INTEGER N_END_TIME;
	QueryPerformanceFrequency(&N_FREQUENCY);
	QueryPerformanceCounter(&N_BEGIN_TIME);

	/************** START THE KLT TRACKING **************/
	
	// Get the tracking image file names.
	imgA_file = m_deqIMG[IMG_A];
	imgB_file = m_deqIMG[IMG_B];

	// Get the corresponding EO of the tracking images.
	eoA = m_deqEO[IMG_A];
	eoB = m_deqEO[IMG_B];

		// DEBUG
		printf("****************************************************\n");
		printf("*****          KLT Image Matching Start        *****\n");
		printf("****************************************************\n");
		printf("Tracking image: [%d] %s\n", eoA.iImg, imgA_file.c_str());
		printf("Tracking image: [%d] %s\n", eoB.iImg, imgB_file.c_str());
		// DEBUG

	// Construct the actual path and filename of the tracking images.
	imgA_file = m_img_directory + imgA_file;
	imgB_file = m_img_directory + imgB_file;

	//START----- PERFORMANCE IMPROVEMENT SOLUTION -----

	// Load the 1st image and 2nd image (A and B) into the CV image structures.
	if (1 == SCALE_FACTOR) // Load the original image size.
	{

		if (m_bReloadIMG == true)
		{
			imgA = cvLoadImage(imgA_file.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
			imgB = cvLoadImage(imgB_file.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
		}
		else // Reuse the initialized image by swapping pointers.
		{
			// Release the existing A image.
			if (imgA != NULL)
			{
				cvReleaseImage(&imgA);
				imgA = NULL;
			}

			// Set the image A ptr to the previous allocation of image B.
			imgA = imgB;
		
			// Load the new image B.
			imgB = cvLoadImage(imgB_file.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
		}

	}
	else // Load the scaling image size.
	{
		if (m_bReloadIMG == true) // Re-read both images.
		{
			IplImage *	srcA;		// Temporary storage for Image#A
			IplImage *	srcB;		// Temporary storage for Image#B
			CvSize		img_sz;		// Image dimension

			// Load the images into the CV image structure.
			srcA = cvLoadImage(imgA_file.c_str(), CV_LOAD_IMAGE_GRAYSCALE);
			srcB = cvLoadImage(imgB_file.c_str(), CV_LOAD_IMAGE_GRAYSCALE);

			// Construct the scaled images.
			imgA = cvCreateImage(cvSize(m_scaledIS.iWidth, m_scaledIS.iHeight), srcA->depth, srcA->nChannels);
			cvResize(srcA, imgA, CV_INTER_LINEAR);

			imgB = cvCreateImage(cvSize(m_scaledIS.iWidth, m_scaledIS.iHeight), srcB->depth, srcB->nChannels);
			cvResize(srcB, imgB, CV_INTER_LINEAR);

			// Release the temporary storages.
			if (srcA != NULL)
				cvReleaseImage(&srcA);
			if (srcB != NULL)
				cvReleaseImage(&srcB);

		}
		else // Reuse the initialized image by swapping pointers.
		{
			// Release the existing A image.
			if (imgA != NULL)
			{
				cvReleaseImage(&imgA);
				imgA = NULL;
			}

			// Set the image A ptr to the previous allocation of image B.
			imgA = imgB;

			// --- Load and scale the new image B.

			IplImage *	srcB;		// Temporary storage for Image#B
			CvSize		img_sz;		// Image dimension

			// Load the images into the CV image structure.
			srcB = cvLoadImage(imgB_file.c_str(), CV_LOAD_IMAGE_GRAYSCALE);

			// Construct the scaled images.
			imgB = cvCreateImage(cvSize(m_scaledIS.iWidth, m_scaledIS.iHeight), srcB->depth, srcB->nChannels);
			cvResize(srcB, imgB, CV_INTER_LINEAR);

			// Release the temporary storages.
			if (srcB != NULL)
				cvReleaseImage(&srcB);
		
		}
	}

	//END----- PERFORMANCE IMPROVEMENT SOLUTION -----

	//START----- EXTRACTING GOOD FEATURES -----

		// DEBUG
		int iNewPoints = 0;
		// DEBUG

	// Initialize the pattern block to store features (Point ID) locating in each block.
	vector<vector<VEC_INT> > vblock_features(SUBPATTERN, vector<VEC_INT>(SUBPATTERN));

	// Initialize the pattern block to store the number of tie points in each block.
	existing_features = new int[SUBPATTERN * SUBPATTERN];

	// Store existing features into each sub-pattern of vblock_features.
	AssignExistingFeaturesToBlock(eoA.iImg, vblock_features, existing_features, n_existing_features);

			// DEBUG
			printf("\t Number of existing points in first image is %d\n", n_existing_features);
			// DEBUG

	for (iH = 0; iH < SUBPATTERN; iH++)
	{
		for (iW = 0; iW < SUBPATTERN; iW++)
		{
			// If the number of features already reaches the MAX_CORNERS, just skip.
			int numberOfFeatures = 0;
			numberOfFeatures = (vblock_features[iH][iW]).size();

			if (numberOfFeatures >= MAX_CORNERS)
			{
				continue;
			}

			// ----- Create ROI for sub-region (iH,iW)

			// Obtain the ROI for subpattern [iH][iW].
			DefineROI(iH, iW, xROI, yROI, widthROI, heightROI);

			// Sub-divide image A to obtain ROI for the region [iH][iW].
			cvSetImageROI(imgA, cvRect(xROI, yROI, widthROI, heightROI));

			// Create a canvas for the ROI on the image A.
			roiA = cvCreateImage(cvGetSize(imgA), imgA->depth, imgA->nChannels);

			// Copy the ROI to the recent created canvas.
			cvCopy(imgA, roiA, NULL);

			// Reset the ROI to recover the original imgA.
			cvResetImageROI(imgA);

			// ----- End of creating ROI for sub-region (iH,iW)

			// ----- Perform the good features to track for the ROI (iH,iW)

			// Obtain the size of the ROI
			CvSize img_sz = cvGetSize( roiA );

			// Perform the good feature to track in the ROI (Image A).
			eig_image = cvCreateImage( img_sz, IPL_DEPTH_32F, 1 );
			tmp_image = cvCreateImage( img_sz, IPL_DEPTH_32F, 1 );
			corner_count = MAX_CORNERS - (vblock_features[iH][iW]).size();
			cornersA = new CvPoint2D32f[ corner_count ];

			_ASSERT(cornersA != NULL);

			if (EXTRACTOR == KLT)
			{
				cvGoodFeaturesToTrack(
						roiA,				// Original image
						eig_image,			// Contains minimal eigenvalues
						tmp_image,
						cornersA,			// Contains results after run
						&corner_count,		// The number of corner
						MINIMUM_EIGENVALUE,	// Minimum eigenvalue
						MIN_DISTANCE,		// Minimum distance: multiple points within small region
						0,
						BLOCK_AUTOCORR,		// Block size for autocorrelation
						0,					// Use Shi and Tomasi
						0.04				// Used by Harris = NOT USED HERE
					);
			}
			else // EXTRACTOR == FAST
			{
				CvPoint* cornersA_FAST;		// Immediated features extracted by FAST
				CvPoint* cornersA_IN;		// FAST features that are in [MIN_DISTANCE] distance
				int		ii, jj, kk;			// Iteration index
				bool	bNonAccept;			// Flag to discard feature

				// Maximum number of points allowed in this subregion 
				NUM_UPPERBOUND = corner_count;	

				// Initialize cornersA_IN
				cornersA_IN = new CvPoint [corner_count];

				// Extract FAST features
				cvCornerFast(imgA, 
						FAST_THRESHOLD, 
						FAST_NUM_FAST_PIXEL, 
						FAST_B_NON_MAX_SUPRESS, 
						&corner_count, 
						&cornersA_FAST
					);

				// Since FAST extracts a huge number of features,
				//   discard features that are closer than FAST_MIN_DISTANCE
				if (corner_count > 0)
				{
					cornersA_IN[0] = cornersA_FAST[0];	// Always accept the first feature
					jj = 1;								// Total number of accepted features
				}

				// Iterate through the list of extracted FAST features
				for (ii = 1; ii < corner_count; ii++)
				{	
					if (jj >= NUM_UPPERBOUND)
						break;

					bNonAccept = false;

					// Iterate through the list of accepted FAST features
					for (kk = 0; kk < jj; kk++)
					{
						int dx = cornersA_FAST[ii].x - cornersA_IN[kk].x;
						int dy = cornersA_FAST[ii].y - cornersA_IN[kk].y;
						int dist = dx * dx + dy * dy;

						// If dist is less than "minimum distance", discard the feature.
						if (dist < FAST_MIN_DISTANCE)
						{
							bNonAccept = true;
							break;
						}
					}

					if (!bNonAccept)
					{
						cornersA_IN[jj] = cornersA_FAST[ii];		
						jj++;
					}
				} // End for ii

				// The final number of FAST features
				corner_count = jj;

				// Convert CvPoint to CvPoint2D32f to be compliant with KLT tracker
				for (ii = 0; ii < corner_count; ii++)
				{
					cornersA[ii] = cvPointTo32f(cornersA_IN[ii]);				
				}

				delete [] cornersA_IN;
				free(cornersA_FAST);
			}

			cvFindCornerSubPix(
					roiA,
					cornersA,
					corner_count,
					cvSize(WIN_SIZE_SUBPIXEL,WIN_SIZE_SUBPIXEL),
					cvSize(-1,-1),
					cvTermCriteria(CV_TERMCRIT_ITER|CV_TERMCRIT_EPS,CRI_NUM_ITERATION,CRI_EPSILON_LIMIT)
				);

			// Get the region ID from the coordinate.
			int regionID = 	FindSubRegionFromID(iH, iW);

			// Validate and store the extracted features.
			for (p=0; p<corner_count; p++)
			{
				// Apply the initial coordinate according to ROI.
				cornersA[p].x = xROI + cornersA[p].x;
				cornersA[p].y = yROI + cornersA[p].y;

				// Continue to next point if either of its values is negative.
				if ( DBL_LT(cornersA[p].x, 1e-10) || DBL_LT(cornersA[p].y, 1e-10) )
					continue;

				// If the point already existed, continue with other points.
				if (IsMatchingPointsExisted(eoA.iImg, cornersA[p]))
					continue;

				// Otherwise, store the points in the matching results.
				TIE_POINT mp;
				mp.pointID = RUNNING_POINT_ID;
				mp.imageID = eoA.iImg;
				mp.regionID = regionID;
				mp.coord = cornersA[p];
				m_vecTP.push_back(mp);
				
				// Increment the world point ID
				RUNNING_POINT_ID++;

					// DEBUG
					iNewPoints++;
					// DEBUG	

				// Store the point in the subregion array as well.
				(vblock_features[iH][iW]).push_back(p);
			}

			// Clear the memory
			delete [] cornersA;
			cornersA = NULL;

			if (eig_image != NULL)
				cvReleaseImage(&eig_image);
			if (tmp_image != NULL)
				cvReleaseImage(&tmp_image);
			if (roiA != NULL)
				cvReleaseImage(&roiA);
		}
	}
			// DEBUG
			printf("\t Number of new points added in first image is %d\n", iNewPoints);

			// DEBUG
			printf("\t Good features to track in the first image\n");
			for (iH = 0; iH < SUBPATTERN; iH++)
			{
				printf("\t\t");
				for (iW = 0; iW < SUBPATTERN; iW++)
				{
					printf("%d\t", (vblock_features[iH][iW]).size());
				}
				printf("\n");
			}
			// DEBUG

	//END----- EXTRACTING GOOD FEATURES -----

	// Create 3D rotation matrix iff
	// 1. Initial guessed tie-point FLAG is set	--- or ---
	// 2. The second image undergoes huge rotation - by *purely* evaluating EO diff
	if (CAL_INITIAL_GUESS == 1)
	{
		Create3DRotMatrix();

		// Compute 2D Similarity transformation between Image#1 and #2
		sim2Dparam = Sim2DTransform();
	}


	// WE HAVE TO KEEP imgB FOR THE NEXT ROUND PROCESSING!!
	// THE REMAINING SOURCECODE OF THIS FUNCTION
	// WILL REFER TO rotated_img.

	/********* The code below is invalid! **********/
	/*
	// Assign the rotated_img back to imgB.
	if (imgB != NULL)
	{
		cvReleaseImage(&imgB);
		imgB = NULL;
	}

	// Set the image B ptr to the rotated image.
	imgB = rotated_img;

	cvReleaseImage(&rotated_img);
	*/
	// ----- Start performing the improved KLT image matching.
	{
		// Double size of the array (*2) to allow the case of MAX_BLOCK_TP

		// Reconstruct the "cornerA" to store the points to be tracked.
		cornersA = new CvPoint2D32f[ (MAX_CORNERS*2) * SUBPATTERN * SUBPATTERN];

		// Reconstruct the array storing the tie points' GPID.
		cornersA_GPID = new int[ (MAX_CORNERS*2) * SUBPATTERN * SUBPATTERN ];

		_ASSERT(cornersA != NULL);

		// Add the significant points in the 1st image to the list to be tracked.
		AddMatchedPointsForTracking(eoA.iImg, &cornersA[0], corner_count, &cornersA_GPID[0]);

		// DEBUG
			printf("\t Total number of features to track is %d\n", corner_count);
		// DEBUG

#ifdef __LSM_DEBUGGING		

#endif //__LSM_DEBUGGING

		cornersB = new CvPoint2D32f[ corner_count ];

		// Define the tracking flag whether or not initial guessed positions exist.
		KLT_TRACKING_FLAG = (CAL_INITIAL_GUESS==0)? 0 : CV_LKFLOW_INITIAL_GUESSES;

		// Construct the tracking parameters
		char * features_found = new char [ corner_count ];
		float * feature_errors = new float [ corner_count ];
		CvSize pyr_sz = cvSize( imgA->width+8, imgA->height/3 );

		pyrA = cvCreateImage( pyr_sz, IPL_DEPTH_32F, 1 );
		pyrB = cvCreateImage( pyr_sz, IPL_DEPTH_32F, 1 );

		/***************** TRACKING **********************/

		if (MOTION_MODEL == TRANSLATION_ROTATED_IMG)
		{ // Translation model tracking on rotated image

			// Rotate image iff its angle is greater than the minimum angle.
			if DBL_GT(fabs(sim2Dparam.angle), MIN_ROTATED_ANGLE)
			{
					printf("\t%s is rotated.\n", imgB_file.c_str());

				CvPoint2D32f src_center;
				rot_matrix = cvCreateMat(2,3,CV_32FC1);
				src_center.x = imgB->width/2.0F;
				src_center.y = imgB->height/2.0F;

				// Compute rotational matrix and store in rot_matrix
				cv2DRotationMatrix(	src_center, 
									sim2Dparam.angle /*degree*/,
									1.0 /*scale*/,
									rot_matrix);

				rotated_img = cvCloneImage(imgB);
				rotated_img->origin = imgB->origin;
				cvZero(rotated_img);
				
				// Warp the result and store in rotated_img
				cvWarpAffine(imgB, rotated_img, rot_matrix);

			}

			if (rotated_img != NULL)
			{
				printf("...Translation tracking with rotated image\n");

				if (CAL_INITIAL_GUESS == 1)
				{
					CalInitialGuesses(	corner_count, cornersA, cornersB, 
										cvGet2D(rot_matrix, 0, 0).val[0],
										cvGet2D(rot_matrix, 0, 1).val[0],
										cvGet2D(rot_matrix, 0, 2).val[0],
										cvGet2D(rot_matrix, 1, 2).val[0]);

					int nOutPt;
					if (COUNT_OUT_SCOPE_TP == 1)
						nOutPt = ReplaceOutOfScopePoints(cornersA, cornersB, cornersA_GPID, corner_count);
					else
						nOutPt = RemoveOutOfScopePoints(cornersA, cornersB, cornersA_GPID, corner_count);

					// DEBUG
						printf("\t Total number of out-of-scope points is %d\n", nOutPt);
					// DEBUG
				}

				cvCalcOpticalFlowPyrLK(
						imgA,					// Image at time t
						rotated_img,			// Image at time t+1
						pyrA,					// Buffer to store input image A
						pyrB,					// Buffer to store input image B
						cornersA,				// Point in image A to be track
						cornersB,				// Point in new location
						corner_count,
						cvSize( WIN_SIZE_TRACK,WIN_SIZE_TRACK ),	// Window size
						DEPTH_LEVELS,			// Level
						features_found,			// status
						feature_errors,
						cvTermCriteria( CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, CRI_NUM_ITERATION, CRI_EPSILON_LIMIT),
						KLT_TRACKING_FLAG//CV_LKFLOW_INITIAL_GUESSES//0
					);

				// Convert the coordinate in rotated image back to original image
				double a11, a12, a21, a22, tx, ty;
				double orix, oriy;
				a11 = cvGet2D(rot_matrix, 0, 0).val[0];
				a12 = cvGet2D(rot_matrix, 0, 1).val[0];
				a21 = cvGet2D(rot_matrix, 1, 0).val[0];
				a22 = cvGet2D(rot_matrix, 1, 1).val[0];
				tx = cvGet2D(rot_matrix, 0, 2).val[0];
				ty = cvGet2D(rot_matrix, 1, 2).val[0];

				for(i=0; i<corner_count; i++ ) 
				{
					orix = a11*(cornersB[i].x - tx) + a21*(cornersB[i].y - ty);
					oriy = a12*(cornersB[i].x - tx) + a22*(cornersB[i].y - ty);

					cornersB[i].x = (float)orix;
					cornersB[i].y = (float)oriy;
				}

				// Release the rotational matrix
				cvReleaseMat(&rot_matrix);

			}
			else // simple perform tracking on non-translated image
			{
				printf("...Translation tracking\n");

				if (CAL_INITIAL_GUESS == 1)
				{
					CalInitialGuesses(corner_count, cornersA, cornersB);

					int nOutPt;
					if (COUNT_OUT_SCOPE_TP == 1)
						nOutPt = ReplaceOutOfScopePoints(cornersA, cornersB, cornersA_GPID, corner_count);
					else
						nOutPt = RemoveOutOfScopePoints(cornersA, cornersB, cornersA_GPID, corner_count);

					// DEBUG
						printf("\t Total number of out-of-scope points is %d\n", nOutPt);
					// DEBUG
				}

				cvCalcOpticalFlowPyrLK(
						imgA,					// Image at time t
						imgB,					// Image at time t+1
						pyrA,					// Buffer to store input image A
						pyrB,					// Buffer to store input image B
						cornersA,				// Point in image A to be track
						cornersB,				// Point in new location
						corner_count,
						cvSize( WIN_SIZE_TRACK,WIN_SIZE_TRACK ),	// Window size
						DEPTH_LEVELS,			// Level
						features_found,			// status
						feature_errors,
						cvTermCriteria( CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, CRI_NUM_ITERATION, CRI_EPSILON_LIMIT),
						KLT_TRACKING_FLAG//CV_LKFLOW_INITIAL_GUESSES//0
					);			
			}

		}
		else if ( (MOTION_MODEL == TRANSLATION_ROTATED_WIN) && (CAL_INITIAL_GUESS==1) )
		{
			if (CAL_INITIAL_GUESS == 1)
			{
				CalInitialGuesses(corner_count, cornersA, cornersB);

				int nOutPt;
				if (COUNT_OUT_SCOPE_TP == 1)
					nOutPt = ReplaceOutOfScopePoints(cornersA, cornersB, cornersA_GPID, corner_count);
				else
					nOutPt = RemoveOutOfScopePoints(cornersA, cornersB, cornersA_GPID, corner_count);

				// DEBUG
					printf("\t Total number of out-of-scope points is %d\n", nOutPt);
				// DEBUG
			}

			if DBL_GT(fabs(sim2Dparam.angle), MIN_ROTATED_ANGLE)
			{
				printf("...Translation tracking with rotated window\n");

				/* Define rotation matrix as used in OpenCV
				 *
				 *		| cos(angle)	sin(angle) |
				 *		| -sin(angle)	cos(angle) |
				 *
				 */

				float * matrices = (float *) calloc(4,sizeof(float));

				matrices[0] = matrices[3] = (float) cos(deg2rad(sim2Dparam.angle));
				matrices[2] = (float) sin(deg2rad(sim2Dparam.angle));
				matrices[1] = -matrices[2];

					// DEBUG
					printf("matrices\n");
					printf("\t%f\t%f\n", matrices[0], matrices[1]);
					printf("\t%f\t%f\n", matrices[2], matrices[3]);
					// DEBUG


				cvCalcRotatedWinFlowPyrLK(
						imgA,					// Image at time t
						imgB,					// Image at time t+1 t
						pyrA,					// Buffer to store input image A
						pyrB,					// Buffer to store input image B
						cornersA,				// Point in image A to be track
						cornersB,				// Point in new location
						matrices,
						corner_count,
						cvSize( WIN_SIZE_TRACK,WIN_SIZE_TRACK ),	// Window size
						DEPTH_LEVELS,			// Level
						features_found,			// status
						feature_errors,
						cvTermCriteria( CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, CRI_NUM_ITERATION, CRI_EPSILON_LIMIT),
						KLT_TRACKING_FLAG//CV_LKFLOW_INITIAL_GUESSES//0
					);

				if (matrices != NULL)
					free(matrices);
			}
			else
			{ // The rotation angle is less than MIN_ROTATED_ANGLE

				printf("...Translation tracking\n");

				cvCalcOpticalFlowPyrLK(
						imgA,					// Image at time t
						imgB,					// Image at time t+1
						pyrA,					// Buffer to store input image A
						pyrB,					// Buffer to store input image B
						cornersA,				// Point in image A to be track
						cornersB,				// Point in new location
						corner_count,
						cvSize( WIN_SIZE_TRACK,WIN_SIZE_TRACK ),	// Window size
						DEPTH_LEVELS,			// Level
						features_found,			// status
						feature_errors,
						cvTermCriteria( CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, CRI_NUM_ITERATION, CRI_EPSILON_LIMIT),
						KLT_TRACKING_FLAG//CV_LKFLOW_INITIAL_GUESSES//0
					);			
			
			}
		
		} // END TRANSLATION_ROTATED_WIN
		else if (MOTION_MODEL == AFFINE)
		{
			printf("...Affine tracking\n");

			/* Define rotation matrix as used in OpenCV
			 *
			 *		| cos(angle)	sin(angle) |
			 *		| -sin(angle)	cos(angle) |
			 *
			 */

			float * matrices;

			if (CAL_INITIAL_GUESS==1)
			{
				CalInitialGuesses(corner_count, cornersA, cornersB);

				int nOutPt;
				if (COUNT_OUT_SCOPE_TP == 1)
					nOutPt = ReplaceOutOfScopePoints(cornersA, cornersB, cornersA_GPID, corner_count);
				else
					nOutPt = RemoveOutOfScopePoints(cornersA, cornersB, cornersA_GPID, corner_count);

				// DEBUG
					printf("\t Total number of out-of-scope points is %d\n", nOutPt);
				// DEBUG

				matrices = (float *) calloc(corner_count*4,sizeof(float));
/*
				matrices[0] = matrices[3] = (float) cos(deg2rad(sim2Dparam.angle));
				matrices[1] = (float) sin(deg2rad(sim2Dparam.angle));
				matrices[2] = -matrices[1];
*/
				matrices[0] = matrices[3] = (float) cos(deg2rad(sim2Dparam.angle));
				matrices[2] = (float) sin(deg2rad(sim2Dparam.angle));
				matrices[1] = -matrices[2];

				for (i = 1; i < corner_count; i++)
				{
					matrices[i*4] = matrices[0];
					matrices[i*4+1] = matrices[1];
					matrices[i*4+2] = matrices[2];
					matrices[i*4+3] = matrices[3];
				}

					// DEBUG
					printf("matrices\n");
					printf("\t%f\t%f\n", matrices[0], matrices[1]);
					printf("\t%f\t%f\n", matrices[2], matrices[3]);
					// DEBUG

			}
			else
			{
				matrices = (float *) calloc(corner_count*4,sizeof(float));
			}

			cvCalcAffineFlowPyrLK(
					imgA,					// Image at time t
					imgB,					// Image at time t+1 t
					pyrA,					// Buffer to store input image A
					pyrB,					// Buffer to store input image B
					cornersA,				// Point in image A to be track
					cornersB,				// Point in new location
					&matrices[0],
					corner_count,
					cvSize( WIN_SIZE_TRACK,WIN_SIZE_TRACK ),	// Window size
					DEPTH_LEVELS,			// Level
					features_found,			// status
					feature_errors,
					cvTermCriteria( CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, CRI_NUM_ITERATION, CRI_EPSILON_LIMIT),
					KLT_TRACKING_FLAG//CV_LKFLOW_INITIAL_GUESSES//0
				);		

			if (matrices)
				free(matrices);

		}
		else // TRANSLATION
		{
			printf("...Translation tracking\n");

			if (CAL_INITIAL_GUESS == 1)
			{
				CalInitialGuesses(corner_count, cornersA, cornersB);

				int nOutPt;
				if (COUNT_OUT_SCOPE_TP == 1)
					nOutPt = ReplaceOutOfScopePoints(cornersA, cornersB, cornersA_GPID, corner_count);
				else
					nOutPt = RemoveOutOfScopePoints(cornersA, cornersB, cornersA_GPID, corner_count);

				// DEBUG
					printf("\t Total number of out-of-scope points is %d\n", nOutPt);
				// DEBUG
			}

			cvCalcOpticalFlowPyrLK(
					imgA,					// Image at time t
					imgB,					// Image at time t+1
					pyrA,					// Buffer to store input image A
					pyrB,					// Buffer to store input image B
					cornersA,				// Point in image A to be track
					cornersB,				// Point in new location
					corner_count,
					cvSize( WIN_SIZE_TRACK,WIN_SIZE_TRACK ),	// Window size
					DEPTH_LEVELS,			// Level
					features_found,			// status
					feature_errors,
					cvTermCriteria( CV_TERMCRIT_ITER | CV_TERMCRIT_EPS, CRI_NUM_ITERATION, CRI_EPSILON_LIMIT),
					KLT_TRACKING_FLAG//CV_LKFLOW_INITIAL_GUESSES//0
				);		
		
		}
			// DEBUG
			int iNumMatchedPoint = 0;
			// DEBUG

		/******************** START OUTLIER REMOVAL ****************************/

		// Remove invalid coordinates
		{
			// Index of the last valid pair of points.
			int iNumValidPoint = -1;

			for(i=0; i<corner_count; i++ ) 
			{
				if( (features_found[i]==0) || DBL_LT(feature_errors[i], CCC_THRESHOLD) || 
					DBL_LT(cornersB[i].x, 1e-10) || DBL_LT(cornersB[i].y, 1e-10) ||
					DBL_GT(cornersB[i].x, m_scaledIS.iWidth) || DBL_GT(cornersB[i].y, m_scaledIS.iHeight) ) 
				{

					#ifdef __LSM_DEBUGGING

						if (features_found[i]==0)
							printf("\tFeatures[%d,GP%d]:B(%lf,%lf) in 1st image _not found_ with error %f\n", 
									i, cornersA_GPID[i], cornersB[i].x, cornersB[i].y, feature_errors[i]);
						
						if DBL_LT(feature_errors[i], CCC_THRESHOLD)
							printf("\tFeatures[%d,GP%d]:B(%lf,%lf) in 1st image found _with CCC %f_\n", 
									i, cornersA_GPID[i], cornersB[i].x, cornersB[i].y, feature_errors[i]);
					
					#endif //__LSM_DEBUGGING

					continue;
				}

				// Increment the index of valid point.
				iNumValidPoint = iNumValidPoint + 1; 

				// If the current index is not the last valid point
				// (which means there is a hole before this current point),
				// re-arrange the cornersA and cornersB array by move this
				// current point up.
				if (i != iNumValidPoint)
				{
					// Move the current point the lastest valid position.
					cornersA[iNumValidPoint].x = cornersA[i].x;
					cornersA[iNumValidPoint].y = cornersA[i].y;
					cornersB[iNumValidPoint].x = cornersB[i].x;
					cornersB[iNumValidPoint].y = cornersB[i].y;
					cornersA_GPID[iNumValidPoint] = cornersA_GPID[i];
				}
			} // end for each corner

			// Update the number of points.
			if (corner_count > 0) // AutomatedAT_2008_8.1.7.9_FIXED
				corner_count = iNumValidPoint+1;

			// DEBUG
				printf("\t Number of valid features is %d\n", corner_count);
			// DEBUG

		}

		if (CCC_CHECK == 1)
		{
			// Index of the last valid pair of points.
			int iNumValidPoint = -1;

			for(i=0; i<corner_count; i++ ) 
			{
				//***** START : Perform Cross Correlation Check *****//
				{
				
					int		n_pixels = 0;	// The number of pixels counted for CCC computation.
					int		xt, yt;			// The x,y index for window size.
					int		valA, valB;		// The pixel value for image A and image B.
					CvScalar s;
					double	d_std_A, d_std_B, d_std_AB;
					double	d_sum_sqare_A = 0.0, d_sum_A = 0.0, d_sum_AxB = 0.0;
					double	d_sum_sqare_B = 0.0, d_sum_B = 0.0;
					double	d_ccc;

					// Looping for variance calculation
					for (yt = -WIN_SIZE_CCC; yt < WIN_SIZE_CCC; yt++)
					{
						for (xt = -WIN_SIZE_CCC; xt < WIN_SIZE_CCC; xt++)
						{
							// Validate the boundary of window.
							if ( DBL_LT(cornersA[i].y + yt, 0.0) ||
								 DBL_LT(cornersA[i].x + xt, 0.0) ||
								 DBL_GT(cornersA[i].y + yt, double(m_scaledIS.iHeight-1)) ||
								 DBL_GT(cornersA[i].x + xt, double(m_scaledIS.iWidth-1)) ||
								 DBL_LT(cornersB[i].y + yt, 0.0) ||
								 DBL_LT(cornersB[i].x + xt, 0.0) ||
								 DBL_GT(cornersB[i].y + yt, double(m_scaledIS.iHeight-1)) ||
								 DBL_GT(cornersB[i].x + xt, double(m_scaledIS.iWidth-1)) )
							{
								continue;
							}

							// Count the number of pixels for CCC computation.
							n_pixels++;

							// Get value of the image A at the matching position.
							s = cvGet2D(imgA, cornersA[i].y + yt, cornersA[i].x + xt);
							valA = (int)s.val[0];

							if (rotated_img != NULL)
							{
							// Get value of the image B at the matching position.
							s = cvGet2D(/*imgB*/rotated_img, cornersB[i].y + yt, cornersB[i].x + xt);
							valB = (int)s.val[0];
							}
							else
							{
							// Get value of the image B at the matching position.
							s = cvGet2D(imgB, cornersB[i].y + yt, cornersB[i].x + xt);
							valB = (int)s.val[0];
							
							}

							// Variance calculation for image A.
							d_sum_sqare_A += valA * valA;
							d_sum_A += valA;

							// Variance calculation for image B.
							d_sum_sqare_B += valB * valB;
							d_sum_B += valB;

							d_sum_AxB += valA * valB;
						}
					}

					if (n_pixels == 0)
						continue;

					// Standard deviation of Image A
					d_std_A		= sqrt((n_pixels * d_sum_sqare_A - d_sum_A*d_sum_A)/n_pixels/(n_pixels-1));

					// Standard deviation of Image B
					d_std_B		= sqrt((n_pixels * d_sum_sqare_B - d_sum_B*d_sum_B)/n_pixels/(n_pixels-1));

					// Covariance between A and B
					d_std_AB	= (n_pixels*d_sum_AxB - d_sum_A*d_sum_B)/n_pixels/(n_pixels-1);

					d_ccc = d_std_AB/(d_std_A*d_std_B);			
					
					//***** END : Perform Cross Correlation Check *****//				

					// Perform validation in comparison with the pre-defined threshold.
					if DBL_LT(d_ccc, CCC_THRESHOLD)
					{
					
						#ifdef	__LSM_DEBUGGING

							printf("\tFeatures[%d, GP%d]:B(x%lf,y%lf) in 1st image _not match_ with ccc %lf\n", 
									i, cornersA_GPID[i], cornersB[i].x, cornersB[i].y, d_ccc);
						
						#endif //__LSM_DEBUGGING

						continue;
					}

				}	// end performing CCC_CHECK for each point

				// Increment the index of valid point.
				iNumValidPoint = iNumValidPoint + 1; 

				// If the current index is not the last valid point
				// (which means there is a hole before this current point),
				// re-arrange the cornersA and cornersB array by move this
				// current point up.
				if (i != iNumValidPoint)
				{
					// Move the current point the lastest valid position.
					cornersA[iNumValidPoint].x = cornersA[i].x;
					cornersA[iNumValidPoint].y = cornersA[i].y;
					cornersB[iNumValidPoint].x = cornersB[i].x;
					cornersB[iNumValidPoint].y = cornersB[i].y;
					cornersA_GPID[iNumValidPoint] = cornersA_GPID[i];
				}

			} // for loop for each corner

			// Update the number of points.
			if (corner_count > 0) // AutomatedAT_2008_8.1.7.9_FIXED
				corner_count = iNumValidPoint+1;

				// DEBUG
					printf("\t Number of features after CCC is is %d\n", corner_count);
				// DEBUG


		} // end if CCC_CHECK

		// AutomatedAT_2008_8.1.7.9.4
		if (OPT_REMOVAL == 1)
		{
			// Outlier removal by optical flow
			if (corner_count > 0) // AutomatedAT_2008_8.1.7.9_FIXED
				OptFlowOutlierRemoval(corner_count, cornersA, cornersB, cornersA_GPID);

			// DEBUG
				printf("\t Number of features after optical flow removal is %d\n", corner_count);
			// DEBUG

		} // end if OPT_REMOVAL

		if (AFFINE_REMOVAL == 1)
		{
			// Outlier removal by AFFINE MODEL
			if (corner_count > 0) // AutomatedAT_2008_8.1.7.9_FIXED
				AffineOutlierRemoval(corner_count, cornersA, cornersB, cornersA_GPID);

			// DEBUG
				printf("\t Number of features after affine removal is %d\n", corner_count);
			// DEBUG
		} // end if AFFINE_REMOVAL

		if (EPIPOLAR_CHECK == 1)
		{	

			// Outlier removal by epipolar line
			if (corner_count > 0) // AutomatedAT_2008_8.1.7.9_FIXED
				EpiOutlierRemoval(corner_count, cornersA, cornersB, cornersA_GPID);

			// DEBUG
				printf("\t Number of features after epipolar removal is %d\n", corner_count);
			// DEBUG


		} // end if EPIPOLAR_CHECK


		// Add the matching results into the matching table.	
		for(i=0; i<corner_count; i++ ) 
		{
			TIE_POINT mp;
			mp.pointID = cornersA_GPID[i];
			mp.imageID = eoB.iImg;
			mp.regionID = FindSubRegionFromCoord(cornersB[i].y, cornersB[i].x);
			mp.coord = cornersB[i];

			// Check the number of matched prior to store the TP to the matched vector.
			if ((MAX_BLOCK_TP == 0) || (vec_num_matched_B[mp.regionID] < MAX_BLOCK_TP))
			{
				m_vecTP.push_back(mp);

				// Increment the number of matched points.
				vec_num_matched_B[mp.regionID] += 1;

				// DEBUG
				iNumMatchedPoint++;
				// DEBUG

					// DEBUG
					{
						// Record the matched points under the image A's view point.
						CvPoint2D32f coord = cornersA[i];
						int regionID = FindSubRegionFromCoord(cornersA[i].y, cornersA[i].x);
						vec_num_matched_A[regionID] += 1;
					}
					// DEBUG
			}
		} // end for each corner


				printf("Total number of matched point is %d\n", iNumMatchedPoint);

			// DEBUG
	
				printf("\t\tMatched point pattern in Image A\n");
				for (r = 0; r < SUBPATTERN; r++)
				{
					printf("\t\t");
					for (c = 0; c < SUBPATTERN; c++)
					{
						i = FindSubRegionFromID(r, c);
						printf("%d\t", vec_num_matched_A[i]);
					}
					printf("\n");
				}

				//---

				printf("\t\tMatched point pattern in Image B\n");
				for (r = 0; r < SUBPATTERN; r++)
				{
					printf("\t\t");
					for (c = 0; c < SUBPATTERN; c++)
					{
						i = FindSubRegionFromID(r, c);
						printf("%d\t", vec_num_matched_B[i]);
					}
					printf("\n");
				}

			// DEBUG

		try
		{
			if (existing_features != NULL)
				delete [] existing_features;
			if (features_found != NULL)
				delete [] features_found;
			if (feature_errors != NULL)
				delete [] feature_errors;
			if (pyrA != NULL)
				cvReleaseImage(&pyrA);
			if (pyrB != NULL)
				cvReleaseImage(&pyrB);
			if (rotated_img != NULL)
				cvReleaseImage(&rotated_img);
			if ((/*iNumMatchedPoint*/corner_count != 0) && (cornersA != NULL))
			{
				delete [] cornersA;
				cornersA = NULL;
			}
			if ((/*iNumMatchedPoint*/corner_count != 0) && (cornersB != NULL))
			{
				delete [] cornersB;
				cornersB = NULL;
			}
		}
		catch(...)
		{
		
		}
	}

	// ----- End of performing the improved KLT image matching.

	/*************** END THE KLT TRACKING ***************/

	// Get the finishing time
	QueryPerformanceCounter(&N_END_TIME);

		printf("\t THE TOTAL TIME USED IS %lf milliseconds\n", 
			(double)(N_END_TIME.QuadPart-N_BEGIN_TIME.QuadPart)*1000/N_FREQUENCY.QuadPart);
}

/********************************************************************
 * Description:	Return the set of tie-point, for the specified IMG_ID
 *				starting from the last stored GP ID, in the form of
 *				the passed-in (by pointer) parameters. 
 *				The function stores the tie points of the specified 
 *				image ID from the m_vecTP vector (declared as VEC_TP) 
 *				into the passed-in vector which is declared as VEC_IP.
 *
 * Return:		The number of tie points stored into the vector.
 ********************************************************************/
int CKLTTracker::GetTiePoints(int IMG_ID, /*out*/VEC_IP * vec_tie_point, int last_stored_GP_ID)
{
	IP					ip;			// Image point 		
	VEC_TP::iterator	iter;		// Iterator for VEC_TP
	int					n_points;	// Total number of tie points.
	CvPoint2D32f		photoCoord;	// Photo coordinate of a point
	FILE *				f;			// Output file.

	// Initialize the number of tie points.
	n_points = 0;

	// Open the file to debug out the TP.txt
	if (true == m_bOutputTP)
	{
		//f = fopen("..\\Matching_Results\\TP.txt", "a");
	}

	// Iterate through the vector of tie points, m_vecTP.


	for (iter = m_vecTP.begin(); iter != m_vecTP.end(); iter++)
	{
		// Store the tie point with the specified image ID into the passed-in vector.
		if ( ((*iter).imageID == IMG_ID) && ((*iter).pointID > last_stored_GP_ID) )
		{
			// Construct the IP object from TP object.
			ip.iImg		= (*iter).imageID;
			ip.iObj		= (*iter).pointID;
			
			// Convert from pixel to photo coordinate.
			photoCoord = ConvertPix2PhotoCoord((*iter).coord);
			
			ip.dX		= photoCoord.x;
			ip.dY		= photoCoord.y;

			// Store the IP object into the passed-in vector.
			vec_tie_point->push_back(ip);

			// Increment the number of tie points.
			n_points++;


			// Debug out tie points into the TP.txt file.
			if (true == m_bOutputTP)
			{
				//fprintf(f, "%d\t%d\t%f\t%f\n", ip.iImg, ip.iObj, (*iter).coord.x * SCALE_FACTOR, (*iter).coord.y * SCALE_FACTOR);
			}

		}
	}

	// Close the TP.txt file after writing.
	if (true == m_bOutputTP)
	{
		//fclose(f);
	}

	return n_points;
}


/********************************************************************
 * Description:	Remove the obsoleted tie points so that the m_vecTP 
 *				vector contains only the set of tie points extracted
 *				from the latest two images.
 ********************************************************************/
int CKLTTracker::RemoveObsoletedTiepoints(int iObsoletedIMG)
{
	int					nObsoletedTP;		// number of obsoleted TP
	VEC_TP::iterator	iter;				// iterator.

	nObsoletedTP = 0;

	// Iterate through the vector of tie points, m_vecTP.
	for (iter = m_vecTP.begin(); iter != m_vecTP.end(); )
	{
		// Remove the obsoleted tie points by checking its image ID.
		if ( (*iter).imageID <= iObsoletedIMG )
		{
			// After 'erase', it returns a pointer to the next valid iterator.
			iter = m_vecTP.erase(iter);
			
			// Count the number of obsoleted TP
			nObsoletedTP++;
		}
		else
			 iter++;
	}

	return nObsoletedTP;
}
