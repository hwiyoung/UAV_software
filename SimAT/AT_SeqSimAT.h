/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _AT_SEQ_SIM_AT_H_
#define	_AT_SEQ_SIM_AT_H_

#pragma once

#include "AT_Definition.h"
#include "AT_GP.h"
#include "AT_ExteriorOrientation.h"
#include "AT_Matrix.h"

using namespace std;

class AT_CSeqSimAT
{
public:
	// Constructor
	AT_CSeqSimAT();

	// Set the member parameters
	void Set(	int nImages, 
				int max_iter, 
				int ninit_sim_img, 
				double dDelta,
				double std_IP,
				double std_GPS,
				double std_INS);

	// Estimate sequential AT simultaneous
	void SeqSimEstimation(	AT_VEC_IP const &vec_ip_f, 
							AT_VEC_INT const &vec_n_IP,
							AT_IO const &io, 
							AT_VEC_EO const &vec_eo);

	// Perform Simultaneous AT tasks
	bool SimultaneousAT(	/*in*/	int n_run_img,
							/*in*/	int n_run_ip,
							/*in*/	AT_VEC_IP const &vec_ip_f, 
							/*in*/	AT_IO const &io, 
							/*in*/	AT_VEC_EO const &vec_eo, 
							/*in*/	AT_CRotationalMatrix const &objRot,
							/*out*/	AT_CExteriorOrientation &objEO,
							/*out*/	AT_CGP &objGP,	/* cover GP_id and Kt_p */
							/*out*/	Matrix &Qee,
							/*out*/	Matrix &Qep,
							/*out*/	Matrix &Qpp,
							/*out*/ int &n_iterations,
							/*out*/ double &d_error);

	// Compute Ae, Ap and y matrices.
	void ComputeAy(	/*in*/	AT_IP const &ip_f,
					/*in*/	AT_IO const &io, 
					/*in*/	AT_EO const &eo, 
					/*in*/	AT_GP const &gp, 
					/*in*/	Matrix const &R,
					/*in*/	Matrix const &dROm,
					/*in*/	Matrix const &dRPh,
					/*in*/	Matrix const &dRKp,
					/*out*/	Matrix &Ae,
					/*out*/	Matrix &Ap,
					/*out*/	Matrix &y);

	// Debug the AT results.
	void DebugSimATResults();

	void GetLine(istream &in, char *buf, const int bufSize) const;

	// Variables declaration.

	int		m_max_iter;			// Maximum number of iterations for performing AT
	int		m_ninit_sim_img;	// Initial number of images for starting seq. sim AT.
	double	m_dDelta;			// Delta: stop criteria
	double	m_std_IP;			// std_IP: std. of the image point errors
	double	m_std_GPS;			// std_GPS: std. of the GPS errors
	double	m_std_INS;			// std_INS: std. of the INS errors
	int		m_nImages;			// Total number of images in the sequence
	int		m_nRound;			// Total running round.
	AT_DEQ_STRING	m_deq_IMG_File;
	// Output of the simultaneous AT
	AT_VEC_DWORD m_vec_elapsed_time;	// The vector of elapsed time for each round.
	AT_VEC_INT	m_vec_iters;			// The vector of no. iterations for each round.
	AT_VEC_CGP	m_vec_objGP;			// Vector of the CGP object containing the adjusted GP.
	AT_VEC_CGP	m_vec_objGPi;			// Vector of the CGP object containing init GP for comparison.
	AT_VEC_CExteriorOrientation m_vec_objEO;	// Vector of the CExteriorOrientation object containing the adjusted EO.
	AT_VEC_BOOL m_vec_result;			// Vector of flags if the Sim AT operation success.
	AT_VEC_DOUBLE m_vec_error;			// Vector of error measured with pre-defined threshold.


};


#endif // _AT_SEQ_SIM_AT_H_
