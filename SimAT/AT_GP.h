/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _AT_GP_H_
#define	_AT_GP_H_


#pragma once

#include "AT_Definition.h"
#include "AT_RotationalMatrix.h"
#include "AT_Bits.h"

using namespace std;

class AT_CGP
{
public:
	// Constructor
	AT_CGP();

	// Clone all member variables.
	AT_CGP Clone(bool bAll);

	// Make the list of ground point ID
	void MakeGPList(AT_VEC_IP const &vec_ip_i, int nimages);

	// Compute the 'rough' initial approximation of GP. 
	void GPInitApproximation(	AT_VEC_IP const &vec_ip_f, 
								AT_IO const &io, 
								AT_VEC_EO const &vec_eo,
								int const nimages,
								AT_CRotationalMatrix const &rot,
								int const n_max_iter);

	// Compute the initial approximation of GP. 
	void GPInitApproximation(	AT_VEC_IP const &vec_ip_f, 
								AT_IO const &io, 
								AT_VEC_EO const &vec_eo,
								int const nimages,
								AT_CRotationalMatrix const &rot,
								int const n_max_iter,
								double const threshold,
								bool b_fine_approx);

	// Compute the Ap and y matrices.
	void ComputeAy(	/*in*/	AT_IP const &ip_f,
					/*in*/	AT_IO const &io, 
					/*in*/	AT_EO const &eo, 
					/*in*/	AT_GP const &gp, 
					/*in*/	Matrix const &R, // Rotational Matrix	
					/*out*/	Matrix &Ap,
					/*out*/	Matrix &y);

	// Find the pair of images that contribute largest distance between PCs.
	bool FindLargestDistancePair(	int iGP, 
									AT_VEC_EO const &eo, 
									int const n_img_sequence, 
									/*out*/ AT_PAIR_INT &pair_img);

	// Return the EO of each image.
	AT_EO GetEO(AT_VEC_EO const &vec_eo, int imageID);

	// Get the image point vector of the specified GP ID & Image ID.
	Matrix GetImagePointVector(		int iGP,
									int imageID,
									AT_VEC_IP const &vec_ip_f, 
									AT_IO const &io);

	// compute GP from a pair of tie points.
	
	AT_GP ComputeGPFromATiePoint(		int iGP,
									AT_PAIR_INT &pair_image,
									AT_VEC_IP const &vec_ip_f, 
									AT_VEC_EO const &vec_eo,
									AT_IO const &io, 
									AT_CRotationalMatrix const &rot);

	// Fill the 'image bits' into the vector of image bits (m_vec_imagebits)
	void FillImageBits(AT_IP const &ip);

	// Return the image point which has the specified Img ID and GP ID.
	AT_IP GetImagePoint(AT_VEC_IP const &vec_ip, int ImageID, int GPID, int iStartHint=0);

	// Return the image point which has the specified Img ID and GP ID.
	int GetImagePointIndex(AT_VEC_IP const &vec_ip, int ImageID, int GPID, int iStartHint=0);

	// Return the index of the specified GP ID in the passed-in vector.
	int GetIndexGP(AT_VEC_INT const &vec_id, int gpID);

	// Get the number of images that this GP exists.
	int GetNumImagesGPExisted(int gpID);

	//	Debug out the distinct object ID in the m_vec_GPID
	void DebugDistinctGP();

	// Debug out the intial approximation of GP.
 	void DebugInitGP(bool bRoughGP = true);

	// Variables declaration.

	int			m_n_distinct_GP;		// Number of distinct ground points
	AT_VEC_INT		m_vec_distinct_GPID;	// Vector of distinct ground point ID
	AT_VEC_INT		m_vec_nimg_GP;			// Vector of number of images that GP exist

	AT_VEC_GP		m_vec_initGP;			// Vector of initial approximation of ground points (X,Y,Z) (rough)
	AT_VEC_GP		m_vec_initGP2;			// Vector of initial approximation of ground points (X,Y,Z) (fine)
	AT_VEC_INT		m_vec_initGPID;			// Vector of index of GP that corresponding to the init approx.

	// Vector with #object point size in which each element contains array of 
	// bits in which '1' means the object exists in that image.
	AT_VEC_DWORD_BITS		m_vec_imagebits;	

	AT_CBits		m_CBits;		// Object of CBits for helper function

};

typedef vector<AT_CGP>	AT_VEC_CGP;		// Vector of the CGP object


#endif // _GP_H_
