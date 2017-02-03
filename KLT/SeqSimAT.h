/********************************************************
 * Project:			Automated AT
 * Last updated:	31 August 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _SEQ_SIM_AT_H_
#define	_SEQ_SIM_AT_H_


#pragma once

#include "Definition.h"
#include "GP.h"
#include "ExteriorOrientation.h"
#include "Matrix.h"

using namespace std;

typedef vector< Matrix > VEC_Matrix;	// Vector of the Matrix class

class CSeqSimAT
{
public:
	// Constructor
	CSeqSimAT();

	// Set the member parameters
	void Set(	int nImages, 
				int max_iter, 
				int ninit_sim_img, 
				int nmax_sim_img,
				double dDelta,
				double std_IP,
				double std_GPS,
				double std_INS,
				int rot_order,
				int rot_type);

	// Assign the initial approximation of EO.
	void EOInitialization(VEC_EO &vec_eo);

	// Estimate sequential AT simultaneous processing based on 'epoch'
	void SeqSimEstimation(	VEC_INT	const &vec_n_EPOCH,	
							VEC_IP const &vec_ip_f, 
							VEC_INT const &vec_n_IP,
							IO const &io, 
							VEC_EO const &vec_eo,
							int latest_image_id);

	// Perform Simultaneous AT tasks
	bool SimultaneousAT(	/*in*/	int start_img_id,
							/*in*/	int last_img_id,
							/*in*/	int n_run_img,
							/*in*/	int idx_start_ip,
							/*in*/	int idx_last_ip,
							/*in*/	int n_run_ip,
							/*in*/	VEC_IP const &vec_ip_f, 
							/*in*/	IO const &io, 
							/*in*/	VEC_EO const &vec_eo, 
							/*in*/	CRotationalMatrix const &objRot,
							/*out*/	CExteriorOrientation &objEO,
							/*out*/	CGP &objGP,	/* cover GP_id and Kt_p */
							/*out*/	Matrix &Qee,
							/*out*/	Matrix &Qep,
							/*out*/	Matrix &Qpp,
							/*out*/ int &n_iterations,
							/*out*/ double &d_error);

	// Compute Ae, Ap and y matrices.
	void ComputeAy(	/*in*/	IP const &ip_f,
					/*in*/	IO const &io, 
					/*in*/	EO const &eo, 
					/*in*/	GP const &gp, 
					/*in*/	Matrix const &R,
					/*in*/	Matrix const &dROm,
					/*in*/	Matrix const &dRPh,
					/*in*/	Matrix const &dRKp,
					/*out*/	Matrix &Ae,
					/*out*/	Matrix &Ap,
					/*out*/	Matrix &y);

	// Return the approximated EO of the specific image.
	EO GetInitEO(int imageID);			

	// Return the index of the approximated EO of the specific image.
	int GetInitEOIndex(int imageID);

	// Debug the AT results.
	void DebugSimATResults();

	// Variables declaration.

	int		m_max_iter;			// Maximum number of iterations for performing AT
	int		m_ninit_sim_img;	// Initial number of images for starting seq. sim AT.
	int		m_nmax_sim_img;		// Maximum number of images for seq. AT computation.
	double	m_dDelta;			// Delta: stop criteria
	double	m_std_IP;			// std_IP: std. of the image point errors
	double	m_std_GPS;			// std_GPS: std. of the GPS errors
	double	m_std_INS;			// std_INS: std. of the INS errors
	int		m_nImages;			// Total number of images in the sequence
	int		m_nRound;			// Total running round.
	int		m_RoundIdx;			// The running round index.
	int		ROT_ORDER;			// Rotational matrix order
								// 1 : Rx*Ry*Rz or 2 : Rz*Ry*Rx
	int		ROT_TYPE;			// Rotational matrix type: 1 or 2
								// Type 1 
								//			Rx = [1 0 0; 0 cos(om) -sin(om); 0 sin(om) cos(om)]
								//			Ry = [cos(ph) 0 sin(ph); 0 1 0; -sin(ph) 0 cos(ph)]
								//			Rz = [cos(kp) -sin(kp) 0; sin(kp) cos(kp) 0; 0 0 1]
								// Type 2
								//			Rx = [1 0 0; 0 cos(om) sin(om); 0 -sin(om) cos(om)]
								//			Ry = [cos(ph) 0 -sin(ph); 0 1 0; sin(ph) 0 cos(ph)]
								//			Rz = [cos(kp) sin(kp) 0; -sin(kp) cos(kp) 0; 0 0 1]

	VEC_INT	vec_idx_eo_pc;		// Vector of image index for EO: PC constraints.
	VEC_INT	vec_idx_eo_at;		// Vector of image index for EO: AC constraints.
	bool	m_bInitConstraint;	// Flag to initialize constraint.

	// Output of the simultaneous AT
	VEC_DWORD m_vec_elapsed_time;	// The vector of elapsed time for each round.
	VEC_INT	m_vec_iters;			// The vector of no. iterations for each round.
	VEC_CGP	m_vec_objGP;			// Vector of the CGP object containing the adjusted GP.
	VEC_CGP	m_vec_objGPi;			// Vector of the CGP object containing init GP for comparison.
	VEC_CExteriorOrientation m_vec_objEO;	// Vector of the CExteriorOrientation object containing the adjusted EO.
	VEC_BOOL m_vec_result;			// Vector of flags if the Sim AT operation success.
	VEC_DOUBLE m_vec_error;			// Vector of error measured with pre-defined threshold.

};


#endif // _SEQ_SIM_AT_H_
