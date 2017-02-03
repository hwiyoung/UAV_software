/********************************************************
 * Project:			Automated AT
 * Last updated:	17 Jaunary 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _GP_H_
#define	_GP_H_


#pragma once

#include "Definition.h"
#include "RotationalMatrix.h"
#include "Bits.h"

using namespace std;

class CGP
{
public:
	// Constructor
	CGP();

	// Constructor with member variables
	CGP(VEC_INT const &vec_n_IP);

	// Clone all member variables.
	CGP Clone(bool bAll);

	// Set the member variables.
	void Set(VEC_INT const &vec_n_IP);

	// Reset all member variables
	void ResetMember(VEC_INT const &vec_n_IP);

	// Make all elements in m_vec_initGP1 and m_vec_initGP2 to be equal
	void ReplicateGP(bool bRough);

	// Make the list of ground point ID
	void MakeGPList(VEC_IP const &vec_ip_i,
					int idx_start_ip,
					int idx_last_ip,
					int start_img_id,
					int last_img_id,
					int nimages);

	// Compute the 'rough' initial approximation of GP. 
	void GPInitApproximation(	VEC_IP const &vec_ip_f, 
								int idx_start_ip,
								int idx_last_ip,
								IO const &io, 
								VEC_EO const &vec_eo,
								int start_img_id,
								int last_img_id,
								int const nimages,
								CRotationalMatrix const &rot,
								int const n_max_iter);

	// Compute the initial approximation of GP. 
	void GPInitApproximation(	VEC_IP const &vec_ip_f, 
								int idx_start_ip,
								int idx_last_ip,
								IO const &io, 
								VEC_EO const &vec_eo,
								int start_img_id,
								int last_img_id,
								int const nimages,
								CRotationalMatrix const &rot,
								int const n_max_iter,
								double const threshold,
								bool b_fine_approx);

	// Compute the Ap and y matrices.
	void ComputeAy(	/*in*/	IP const &ip_f,
					/*in*/	IO const &io, 
					/*in*/	EO const &eo, 
					/*in*/	GP const &gp, 
					/*in*/	Matrix const &R, // Rotational Matrix	
					/*out*/	Matrix &Ap,
					/*out*/	Matrix &y);

	// Find the pair of images that contribute largest distance between PCs.
	bool FindLargestDistancePair(	int GPID, 
									VEC_EO const &vec_eo, 
									int start_img_id,
									int last_img_id,
									int const n_img_sequence, 
									/*out*/ PAIR_INT &pair_img);

	// Return the EO of each image.
	EO GetEO(VEC_EO const &vec_eo, int imageID);

	// Get the image point vector of the specified GP ID & Image ID.
	Matrix GetImagePointVector(	int GPID,
								int imageID,
								VEC_IP const &vec_ip_f, 
								IO const &io);

	Matrix GetImagePointVector(	int GPID,
								int imageID,
								VEC_IP const &vec_ip_f, 
								int idx_start_ip,
								int idx_last_ip, 
								IO const &io);

	// compute GP from a pair of tie points.
	GP ComputeGPFromATiePoint(	int GPID,
								PAIR_INT &pair_image,
								VEC_IP const &vec_ip_f, 
								VEC_EO const &vec_eo,
								int idx_start_ip,
								int idx_last_ip, 
								IO const &io, 
								CRotationalMatrix const &rot);

	// Fill the 'image bits' into the vector of image bits (m_vec_imagebits)
	void FillImageBits(IP const &ip);

	// Return the image point which has the specified Img ID and GP ID.
	IP GetImagePoint(VEC_IP const &vec_ip, int ImageID, int GPID, int iStartHint=0);

	// Return the image point which has the specified Img ID and GP ID.
	int GetImagePointIndex(VEC_IP const &vec_ip, int ImageID, int GPID, int iStartHint=0);

	// Return the index of the specified GP ID in the passed-in vector.
	int GetIndexGP(VEC_INT const &vec_id, int gpID);

	// Return the index of the specified GP ID in the distinct GP vector.
	int GetIndexDistinctGP(int gpID);

	// Return the index of the specified GP ID from the list of GP Init.
	int GetIndexInitGP(int gpID);

	// Get the number of images that this GP exists.
	int GetNumImagesGPExisted(int gpID);

	//	Debug out the distinct object ID in the m_vec_GPID
	void DebugDistinctGP();

	// Debug out the intial approximation of GP.
 	void DebugInitGP(bool bRoughGP = true);

	// Variables declaration.

	int			m_n_distinct_GP;		// Number of distinct ground points
	VEC_INT		m_vec_distinct_GPID;	// Vector of distinct ground point ID
	VEC_INT		m_vec_nimg_GP;			// Vector of number of images that GP exist

	VEC_GP		m_vec_initGP;			// Vector of initial approximation of ground points (X,Y,Z) (rough)
	VEC_GP		m_vec_initGP2;			// Vector of initial approximation of ground points (X,Y,Z) (fine)
	VEC_INT		m_vec_initGPID;			// Vector of index of GP that corresponding to the init approx.

	VEC_INT		m_vec_n_validGP;		// Vector of number of valid GPs (TP shared > 2 images).
										// The index of vector reflects the real image ID.
	
	// Vector with #object point size in which each element contains array of 
	// bits in which '1' means the object exists in that image.
	VEC_DWORD_BITS		m_vec_imagebits;	

	int					m_base_img_id;	// The based Image ID for m_vec_imagebits.
	int					m_base_gp_id;	// The based GP ID.

	CBits		m_CBits;				// Object of CBits for helper function

};

typedef vector<CGP>	VEC_CGP;		// Vector of the CGP object


#endif // _GP_H_
