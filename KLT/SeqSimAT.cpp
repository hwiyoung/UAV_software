/********************************************************
 * Project:			Automated AT
 * Last updated:	31 August 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "SeqSimAT.h"
#include "Matrix.h"
#include "stdio.h"
#include <iostream>
#include <CRTDBG.H>	// To use _ASSERT

/********************************************************************
 * Description:	Construction
 ********************************************************************/
CSeqSimAT::CSeqSimAT()
{
	// Initialize the running round index.
	m_RoundIdx = 0;

	// Initialize the flag that constraints have not been initialized yet.
	m_bInitConstraint = false;
}

/********************************************************************
 * Description:	Set the member parameters.
 ********************************************************************/
void CSeqSimAT::Set(int nImages, 
					int max_iter, 
					int ninit_sim_img, 
					int nmax_sim_img,
					double dDelta,
					double std_IP,
					double std_GPS,
					double std_INS,
					int rot_order,
					int rot_type)
{
	m_nImages		=	nImages;
	m_max_iter		=	max_iter;
	m_ninit_sim_img	=	ninit_sim_img;
	m_nmax_sim_img	=	nmax_sim_img;
	m_dDelta		=	dDelta;	
	m_std_IP		=	std_IP;
	m_std_GPS		=	std_GPS;
	m_std_INS		=	std_INS;
	ROT_ORDER		=	rot_order;
	ROT_TYPE		=	rot_type;
}

/********************************************************************
 * Description:	Estimate sequential AT simultaneous processing based
 *				on 'epoch'
 ********************************************************************/
void CSeqSimAT::SeqSimEstimation(	VEC_INT	const &vec_n_EPOCH,	
									VEC_IP const &vec_ip_f, 
									VEC_INT const &vec_n_IP,
									IO const &io, 
									VEC_EO const &vec_eo,
									int latest_image_id)
{
	int round;			// Round of iteration.
	int	n_run_img;		// The last number of image in each round.
	int	n_run_ip;		// Number of images points in processing.
	int	n_images;		// The total number of images in computation.
	DWORD start, end;	// Measuring computational time.
	int	id_start_img;	// Starting image ID to run sequential AT.
	int id_start_ip;	// Starting ip index in the 'vec_ip_f' vector.
	int id_last_ip;		// Last ip index in the 'vec_ip_f' vector for computation.
	int i;
	bool bSuccess;

	// The round of processing.
	round = m_RoundIdx;

	// Assign the latest image ID to the variable to maintain the existing code.
	n_run_img = latest_image_id;

	// Obtain the starting Image ID in which only 'max. number of images' are in computation.
	id_start_img	= 1;
	if (n_run_img > m_nmax_sim_img)
		id_start_img = n_run_img - m_nmax_sim_img + 1;

	// Get the index of the starting image point in the 'vec_ip_f' vector.
	id_start_ip = 0;
	for (i = 1; i < id_start_img; i++)			// Image index starts at 1.
		id_start_ip = id_start_ip + vec_n_IP[i];

	// Get the index of the last image point for computation.
	id_last_ip = id_start_ip;
	for (i = id_start_img; i <= n_run_img; i++)
		id_last_ip = id_last_ip + vec_n_IP[i];
	
	id_last_ip = id_last_ip - 1;

	// Get the number of image points in the processing.
	n_run_ip = id_last_ip - id_start_ip;

	// Get the number of images actually used in the processing.
	n_images = n_run_img - id_start_img + 1;

	//********* Specify the index of images for EO constraints *********//

	if (m_bInitConstraint == false)
	{
		// Initialize all constraints.
		for (i = id_start_img; i <= n_run_img; i++)
		{
			// EO constraint : Perspective center
			vec_idx_eo_pc.push_back(i);
	
			// EO constraint : Attitude
			vec_idx_eo_at.push_back(i);
		}	

		m_bInitConstraint = true;
	}
	else
	{
		// Add the current image to be a constraint.
		vec_idx_eo_pc.push_back(n_run_img);
		vec_idx_eo_at.push_back(n_run_img);
	}

	// Remove the obsolete image constraint when it reached m_nmax_sim_img.
	
	// Remove if the first image in the AT constraint list is obsolete.
	i = vec_idx_eo_at[0];		
	if (i < id_start_img)
		vec_idx_eo_at.erase( vec_idx_eo_at.begin()); 

	// Remove if the first image in the PC constraint list is obsolete.
	i = vec_idx_eo_pc[0];		
	if (i < id_start_img)
		vec_idx_eo_pc.erase( vec_idx_eo_pc.begin()); 
	
	//***************** Performing simultaneous AT *****************//

	printf("====================================================\n");
	printf("          Simultaneous AT processing with %d images.\n", round, n_images);
	printf("          Starting from Image#%d to Image#%d.\n", id_start_img, n_run_img);
	printf("====================================================\n");

	Matrix	Qee, Qep, Qpp;
	int		n_iterations = 0;
	double	d_error = -1.0;
	CGP						objGP(vec_n_IP);	
	CExteriorOrientation	objEO(vec_eo, n_images, false, id_start_img);
	CRotationalMatrix		objRot;
	objRot.Compute3DRotMatrix(objEO.m_vec_eo, objEO.m_nImages, ROT_ORDER, ROT_TYPE);
	objRot.ComputeDerivative3DRotMatrix(objEO.m_vec_eo, objEO.m_nImages, ROT_ORDER, ROT_TYPE);

	// Get the tick of starting time.
	start = GetTickCount();

	// Perform Simultaneous AT
	bSuccess = SimultaneousAT(	/*in*/	id_start_img,
								/*in*/	n_run_img,
								/*in*/	n_images,
								/*in*/	id_start_ip,
								/*in*/	id_last_ip,
								/*in*/	n_run_ip,
								/*in*/	vec_ip_f, 
								/*in*/	io, 
								/*in*/	vec_eo, 
								/*in*/	objRot,
								/*out*/	objEO,
								/*out*/	objGP,	/* cover GP_id and Kt_p */
								/*out*/	Qee,
								/*out*/	Qep,
								/*out*/	Qpp,
								/*out*/ n_iterations,
								/*out*/ d_error);

	if (bSuccess == false)
	{
		// Remove the last element from the constraint list.			
		vec_idx_eo_at.erase( vec_idx_eo_at.begin() + (vec_idx_eo_at.size()-1) ); 	

		// Reset the output EO and GP objects.
		objEO.ResetMember(vec_eo, n_images, false, id_start_img);
		objGP.ResetMember(vec_n_IP);

		// Perform Simultaneous AT
		bSuccess = SimultaneousAT(	/*in*/	id_start_img,
									/*in*/	n_run_img,
									/*in*/	n_images,
									/*in*/	id_start_ip,
									/*in*/	id_last_ip,
									/*in*/	n_run_ip,
									/*in*/	vec_ip_f, 
									/*in*/	io, 
									/*in*/	vec_eo, 
									/*in*/	objRot,
									/*out*/	objEO,
									/*out*/	objGP,	/* cover GP_id and Kt_p */
									/*out*/	Qee,
									/*out*/	Qep,
									/*out*/	Qpp,
									/*out*/ n_iterations,
									/*out*/ d_error);

		if (bSuccess == false)
		{
			// Remove the last image from the constraint list.	
			vec_idx_eo_pc.erase( vec_idx_eo_pc.begin() + (vec_idx_eo_pc.size()-1) );  	

			// Reset the output EO and GP objects.
			objEO.ResetMember(vec_eo, n_images, false, id_start_img);
			objGP.ResetMember(vec_n_IP);

			// Perform Simultaneous AT
			bSuccess = SimultaneousAT(	/*in*/	id_start_img,
										/*in*/	n_run_img,
										/*in*/	n_images,
										/*in*/	id_start_ip,
										/*in*/	id_last_ip,
										/*in*/	n_run_ip,
										/*in*/	vec_ip_f, 
										/*in*/	io, 
										/*in*/	vec_eo, 
										/*in*/	objRot,
										/*out*/	objEO,
										/*out*/	objGP,	/* cover GP_id and Kt_p */
										/*out*/	Qee,
										/*out*/	Qep,
										/*out*/	Qpp,
										/*out*/ n_iterations,
										/*out*/ d_error);
		}
	}

	// Get the tick of end time.
	end = GetTickCount();

	// Record the Sequential AT result for further processing.
	// NOTE: This may slow down the process

			// Store the flag if the Sim AT operation success or not.
			m_vec_result.push_back(bSuccess);

			// Store the adjusted EO (CExteriorOrientation) into the vector.
			m_vec_objEO.push_back(objEO);

			// Store the adjusted GP (CGP object) into the vector.
			m_vec_objGP.push_back(objGP);

			// Store the no. iterations for this round to the vector.
			m_vec_iters.push_back(n_iterations);

			// Store the error wrt the pre-defined threshold to the vector.
			m_vec_error.push_back(d_error);

			// Store the elapsed time into the vector.
			DWORD elapsed_time = end - start;
			m_vec_elapsed_time.push_back(elapsed_time);

			// Compute GP init for comparison.
			CGP objGPi(vec_n_IP);
			
			objGPi.MakeGPList(	vec_ip_f, 
								id_start_ip, 
								id_last_ip, 
								id_start_img, 
								n_run_img, 
								n_images);

			objGPi.GPInitApproximation(	vec_ip_f, 
										id_start_ip, 
										id_last_ip, 
										io, 
										vec_eo, 
										id_start_img, 
										n_run_img, 
										n_images, 
										objRot, 
										m_max_iter, 
										m_dDelta, 
										true);

			// Store the GP from initial approximation for comparison.
			m_vec_objGPi.push_back(objGPi);

	// End of recording.

	printf("Round %d consumes %d milliseconds for processing.\n\n", round, elapsed_time);

	// Increment the running round index.
	m_RoundIdx++;

	m_nRound = m_RoundIdx;
}

/********************************************************************
 * Description:	Perform Simultaneous AT tasks.
 ********************************************************************/
bool CSeqSimAT::SimultaneousAT(	/*in*/	int start_img_id,
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
								/*out*/ double &d_error)
{
	// Variable declaration
	Matrix	Nee, Nep, Npp, Ce, Cp, iNpp, Q;		// Normal matrices.
	int r, c, k, n_IP, n_total_GP, n_IM;
	int	imgID, gpID, idxIM, idxGP;
	int idx;
	int n_idx_eo_pc, n_idx_eo_at;
	double thCnd;								// Matrix condition factor.

	// Make the list of distinct ground point ID
	objGP.MakeGPList(vec_ip_f, idx_start_ip, idx_last_ip, start_img_id, last_img_id, n_run_img);

		#ifdef __LSM_DEBUGGING

			// Debug out the list of GP.
			objGP.DebugDistinctGP();

		#endif // __LSM_DEBUGGING

	// Initial approximation for GPs
	bool bFineApprox = true; // false

	objGP.GPInitApproximation(	vec_ip_f, 
								idx_start_ip, 
								idx_last_ip,
								io, 
								vec_eo, 
								start_img_id,
								last_img_id,
								n_run_img, 
								objRot, 
								m_max_iter, 
								m_dDelta, 
								bFineApprox);


		#ifdef __LSM_DEBUGGING

			// Debug out the inital approximated GPs.
			objGP.DebugInitGP(true);
			objGP.DebugInitGP(false);

		#endif //__LSM_DEBUGGING

	// Get the number of GPs that exists more than 2 images.
	n_total_GP = objGP.m_vec_initGP.size();

	if (n_total_GP == 0)
	{
		printf("Operation Failed!\n");
		printf("Not a single pair of tie points existed.\n");
		
		return false;
	}

	// Store the fine approximated GP to enable recovery.
	objGP.ReplicateGP(true);

	// Iterate until the error is smaller than the delta threshold or maximum iterations is reached.
	for (k = 0; k < m_max_iter; k++)
	{
		printf("\tIteration %d\n", k);

		// Initialize the normal matrix and projected observations.
		Nee.SetDim(n_run_img * N_UNKNOWN_EO, n_run_img * N_UNKNOWN_EO);
		Nep.SetDim(n_run_img * N_UNKNOWN_EO, n_total_GP * N_UNKNOWN_GP);
		Npp.SetDim(n_total_GP * N_UNKNOWN_GP, n_total_GP * N_UNKNOWN_GP);
		Ce.SetDim(n_run_img * N_UNKNOWN_EO, 1);
		Cp.SetDim(n_total_GP * N_UNKNOWN_GP, 1);
		iNpp.SetDim(n_total_GP * N_UNKNOWN_GP, n_total_GP * N_UNKNOWN_GP);
		Nee.Clear();
		Nep.Clear();
		Npp.Clear();
		Ce.Clear();
		Cp.Clear();
		iNpp.Clear();

		// Define vector of Ay matrices for each IP 
		VEC_Matrix vec_Ae, vec_Ap, vec_y;

		// Define vector for EO's PC and EO's attitude constraints.
		VEC_Matrix vec_yepi, vec_yeai;

		// Initialize the rotational matrix
		CRotationalMatrix rot;
		rot.Compute3DRotMatrix(objEO.m_vec_eo, objEO.m_nImages, ROT_ORDER, ROT_TYPE);
		rot.ComputeDerivative3DRotMatrix(objEO.m_vec_eo, objEO.m_nImages, ROT_ORDER, ROT_TYPE);
			
		// ***** Compute the normal matrix ***** //

		// Iterate through each image point until the last point in of the current image.
		for (n_IP = idx_start_ip; n_IP <= idx_last_ip; n_IP++)
		{
			Matrix	Ae, Ap, y(2,1);
			Matrix	M;			// Temporary
			double	d;			// Temporary
			EO		eo;			// Temporary
			GP		gp;			// Temporary

			gpID	= vec_ip_f[n_IP].iObj;			// GP ID

			// Continue if this GP ID exists less than two images.
			if (objGP.GetNumImagesGPExisted(gpID) < 2)
				continue;
			
			imgID	= vec_ip_f[n_IP].iImg;				// Image ID
			idxIM	= objEO.GetEOIndex(imgID);			// Index of Image ID
			idxGP	= objGP.GetIndexInitGP(gpID);		// Index of the GPID.
			gp		= objGP.m_vec_initGP2[idxGP];		// GP coordinates.
			eo		= objEO.GetEO(imgID);				// EO of Image ID

			// Compute the Ae, Ap and y matrices.
			ComputeAy(	/*in*/	vec_ip_f[n_IP],
						/*in*/	io, 
						/*in*/	eo, 
						/*in*/	gp, 
						/*in*/	rot.GetRotMatrix(imgID),
						/*in*/	rot.GetDOmegaMatrix(imgID),
						/*in*/	rot.GetDPhiMatrix(imgID),
						/*in*/	rot.GetDKappaMatrix(imgID),
						/*out*/	Ae,
						/*out*/	Ap,
						/*out*/	y);

			// Record the computed results into the defined vectors.
			vec_Ae.push_back(Ae);
			vec_Ap.push_back(Ap);
			vec_y.push_back(y);

			// Nee(imi*6-5:imi*6,imi*6-5:imi*6) += Ae_i' * Ae_i / std_IP ^ 2

			M = Ae.Trans() * Ae * (1.0 / m_std_IP / m_std_IP);

			for (r = 0; r < N_UNKNOWN_EO; r++)
			{
				for (c = 0; c < N_UNKNOWN_EO; c++)
				{
					d = Nee.Get(idxIM*N_UNKNOWN_EO + r, idxIM*N_UNKNOWN_EO + c) + M.Get(r, c);
					Nee.Set(idxIM*N_UNKNOWN_EO + r, idxIM*N_UNKNOWN_EO + c, d);
				}
			}

			// Npp(gpi*3-2:gpi*3,gpi*3-2:gpi*3) += Ap_i' * Ap_i / std_IP ^ 2

			M = Ap.Trans() * Ap * (1.0 / m_std_IP / m_std_IP);

			for (r = 0; r < N_UNKNOWN_GP; r++)
			{
				for (c = 0; c < N_UNKNOWN_GP; c++)
				{
					d = Npp.Get(idxGP*N_UNKNOWN_GP + r, idxGP*N_UNKNOWN_GP + c) + M.Get(r, c);
					Npp.Set(idxGP*N_UNKNOWN_GP + r, idxGP*N_UNKNOWN_GP + c, d);
				}
			}

			// Nep(imi*6-5:imi*6,gpi*3-2:gpi*3) = Ae_i' * Ap_i / std_IP ^ 2;

			M = Ae.Trans() * Ap * (1.0 / m_std_IP / m_std_IP);
			
			for (r = 0; r < N_UNKNOWN_EO; r++)
			{
				for (c = 0; c < N_UNKNOWN_GP; c++)
				{
					Nep.Set(idxIM*N_UNKNOWN_EO + r, idxGP*N_UNKNOWN_GP + c, M.Get(r, c));
				}
			}			

			// Ce(imi*6-5:imi*6,1) += Ae_i' * y_i / std_IP ^ 2;

			M = Ae.Trans() * y * (1.0 / m_std_IP / m_std_IP);

			for (r = 0; r < N_UNKNOWN_EO; r++)
			{
				d = Ce.Get(idxIM*N_UNKNOWN_EO + r, 0) + M.Get(r, 0);
				Ce.Set(idxIM*N_UNKNOWN_EO + r, 0, d);
			}

			// Cp(gpi*3-2:gpi*3,1) += Ap_i' * y_i / std_IP ^ 2;

			M = Ap.Trans() * y * (1.0 / m_std_IP / m_std_IP);

			for (r = 0; r < N_UNKNOWN_GP; r++)
			{
				d = Cp.Get(idxGP*N_UNKNOWN_GP + r, 0) + M.Get(r, 0);
				Cp.Set(idxGP*N_UNKNOWN_GP + r, 0, d);
			}

		}

		// ***** Consider EO constraints : EO's Perspective Center ***** //

		// Get the number of images in which their EO's perspective centers will be included.
		n_idx_eo_pc = vec_idx_eo_pc.size();

		// Iterate through each image according to the index in 'vec_idx_eo_pc'.
		for (idx = 0; idx < n_idx_eo_pc; idx++)
		{
			Matrix	M, N, ye(3,1);	// Temporary
			double	d;				// Temporary
			EO		eo_m, eo_adj;

			// Get 'Image ID' from the vector.
			n_IM = vec_idx_eo_pc[idx];

			// Get index of the Image ID
			idxIM = objEO.GetEOIndex(n_IM);

			// Nee(ni*6-5:ni*6-3,ni*6-5:ni*6-3) += eye(3) / std_GPS ^ 2;
			M = N.Eye(3) * (1.0 / m_std_GPS / m_std_GPS);

			for (r = 0; r < N_UNKNOWN_GP; r++)
			{
				for (c = 0; c < N_UNKNOWN_GP; c++)
				{
					d = Nee.Get(idxIM*N_UNKNOWN_EO + r, idxIM*N_UNKNOWN_EO + c) + M.Get(r, c);
					Nee.Set(idxIM*N_UNKNOWN_EO + r, idxIM*N_UNKNOWN_EO + c, d);
				}
			}

			// yepi{n} = EOi(imi,2:4)' - KTe(imi*6-5:imi*6-3,1);
			eo_m	= objGP.GetEO(vec_eo, n_IM);		// EO from the direct measurement.
			eo_adj	= objEO.GetEO(n_IM);				// EO from approximation process.	

			ye.Set(0, 0, eo_m.dXc - eo_adj.dXc);
			ye.Set(1, 0, eo_m.dYc - eo_adj.dYc);
			ye.Set(2, 0, eo_m.dZc - eo_adj.dZc);

			// Record ye into the resulting vector.
			vec_yepi.push_back(ye);

			// Ce(imi*6-5:imi*6-3,1) += yepi{n} / stdGPS ^ 2;
			for (r = 0; r < N_UNKNOWN_GP; r++)
			{
				d = Ce.Get(idxIM*N_UNKNOWN_EO + r, 0) + ye.Get(r, 0)/m_std_GPS/m_std_GPS;
				Ce.Set(idxIM*N_UNKNOWN_EO + r, 0, d);
			}
		}

		// ***** Consider EO constraints : EO's Attitude ***** //

		// Get the number of images in which their EO's attitudes will be included.
		n_idx_eo_at = vec_idx_eo_at.size();

		// Iterate through each image according to the index in 'vec_idx_eo_at'.
		for (idx = 0; idx < n_idx_eo_at; idx++)
		{
			Matrix	M, N, ye(3,1);	// Temporary
			double	d;				// Temporary
			EO		eo_m, eo_adj;

			// Get 'Image ID' from the vector.
			n_IM = vec_idx_eo_at[idx];

			// Get index of the Image ID
			idxIM = objEO.GetEOIndex(n_IM);

			// Nee(ni*6-2:ni*6,ni*6-2:ni*6) += eye(3) / std_INS ^ 2;
			M = N.Eye(3) * (1.0 / m_std_INS / m_std_INS);

			for (r = 3; r < N_UNKNOWN_EO; r++)
			{
				for (c = 3; c < N_UNKNOWN_EO; c++)
				{
					d = Nee.Get(idxIM*N_UNKNOWN_EO + r, idxIM*N_UNKNOWN_EO + c) + M.Get(r-N_UNKNOWN_GP, c-N_UNKNOWN_GP);
					Nee.Set(idxIM*N_UNKNOWN_EO + r, idxIM*N_UNKNOWN_EO + c, d);
				}
			}

			// yepi{n} = EOi(imi,2:4)' - KTe(imi*6-5:imi*6-3,1);
			eo_m	= objGP.GetEO(vec_eo, n_IM);		// EO from the direct measurement.
			eo_adj	= objEO.GetEO(n_IM);				// EO from approximation process.	

			ye.Set(0, 0, eo_m.dOmega - eo_adj.dOmega);
			ye.Set(1, 0, eo_m.dPhi - eo_adj.dPhi);
			ye.Set(2, 0, eo_m.dKappa - eo_adj.dKappa);

			// Record ye into the resulting vector.
			vec_yeai.push_back(ye);

			// Ce(imi*6-2:imi*6,1) += yeai{n} / stdINS ^ 2;
			for (r = 3; r < N_UNKNOWN_EO; r++)
			{
				d = Ce.Get(idxIM*N_UNKNOWN_EO + r, 0) + ye.Get(r - N_UNKNOWN_GP, 0)/m_std_INS/m_std_INS;
				Ce.Set(idxIM*N_UNKNOWN_EO + r, 0, d);
			}
		}

		// Iterate through each ground point index.
		for (gpID = 0; gpID < n_total_GP; gpID++)
		{
			Matrix	M(N_UNKNOWN_GP, N_UNKNOWN_GP);	// Temporary
			double	d;								// Temporary

			// Npp(n*3-2:n*3,n*3-2:n*3)
			for (r = 0; r < N_UNKNOWN_GP; r++)	
			{
				for (c = 0; c < N_UNKNOWN_GP; c++)	
				{
					d = Npp.Get(gpID*N_UNKNOWN_GP + r, gpID*N_UNKNOWN_GP + c);
					M.Set(r, c, d);
				}
			}

			// Evaluate the condition of partial Npp matrix.
			thCnd = M.EuclideanNorm() * (M.Inv()).EuclideanNorm();
			if DBL_GT(thCnd,THRESHOLD_CND)
			{
				return false;
			}
			
			M = M.Inv();

			// iNpp = inv(Npp(n*3-2:n*3,n*3-2:n*3))
			for (r = 0; r < N_UNKNOWN_GP; r++)	
			{
				for (c = 0; c < N_UNKNOWN_GP; c++)	
				{
					d = M.Get(r, c);
					iNpp.Set(gpID*N_UNKNOWN_GP + r, gpID*N_UNKNOWN_GP + c, d);
				}
			}
		}

		Matrix	Nr, R, iR;
		Matrix	kt_e, kt_g;

		// Q1 = Nep * iNpp;
		Q = Nep * iNpp;

		// Nr = ( Nee - Q1 * Nep' );
		Nr = Nee - Q * Nep.Trans();

		// Evaluate the condition of Nr matrix.
		thCnd = Nr.EuclideanNorm() * (Nr.Inv()).EuclideanNorm();
		if DBL_GT(thCnd,THRESHOLD_CND)
		{
			return false;
		}

		// Qee = inv(Nr)
		Qee = Nr.Inv();
		
		// kt_e = Qee * ( Ce - Q1 * Cp );
		kt_e = Qee * (Ce - Q * Cp);

		// kt_g = iNpp * Cp - Q1' * kt_e;
		kt_g = iNpp * Cp - Q.Trans() * kt_e;

		// Update the adjusted EOs.		
		for (idxIM = 0; idxIM < n_run_img; idxIM++)
		{
			// Kt_e = Kt_e + kt_e
			objEO.m_vec_eo[idxIM].dXc		+= kt_e.Get(idxIM*N_UNKNOWN_EO + 0, 0);
			objEO.m_vec_eo[idxIM].dYc		+= kt_e.Get(idxIM*N_UNKNOWN_EO + 1, 0);
			objEO.m_vec_eo[idxIM].dZc		+= kt_e.Get(idxIM*N_UNKNOWN_EO + 2, 0);
			objEO.m_vec_eo[idxIM].dOmega	+= kt_e.Get(idxIM*N_UNKNOWN_EO + 3, 0);
			objEO.m_vec_eo[idxIM].dPhi		+= kt_e.Get(idxIM*N_UNKNOWN_EO + 4, 0);
			objEO.m_vec_eo[idxIM].dKappa	+= kt_e.Get(idxIM*N_UNKNOWN_EO + 5, 0);
		}

		// Update the adjusted GPs.		
		for (gpID = 0; gpID < n_total_GP; gpID++)
		{
			// Kt_p = Kt_p + kt_g
			objGP.m_vec_initGP2[gpID].X += kt_g.Get(gpID*N_UNKNOWN_GP + 0, 0);
			objGP.m_vec_initGP2[gpID].Y += kt_g.Get(gpID*N_UNKNOWN_GP + 1, 0);
			objGP.m_vec_initGP2[gpID].Z += kt_g.Get(gpID*N_UNKNOWN_GP + 2, 0);
		}
		
		d_error = kt_g.Norm(0, false)/n_total_GP/3.0;
		//printf("\t\tUpdated error %lf\n", d_error);
		if ( DBL_LT( d_error, m_dDelta ) )
			break;
	}

	// no_IT = k;
	n_iterations = k;

	if (n_iterations >= m_max_iter)
		return false;

	// Qep = -Qee * Q1;
	Qep = (Qee * (-1.0)) * Q;

	// Qpp = iNpp + Q1' * Qee * Q1;
	Qpp = iNpp + (Q.Trans()* Qee * Q);

		#ifdef __LSM_DEBUGGING

			// Print out the adjusted EOs and GPs.

			printf("\tTotal running iteration : %d\n", n_iterations);		
			printf("\tThe adjusted EO parameters:\n");
			for (idxIM = 0; idxIM < n_run_img; idxIM++)
			{
				printf("\t[%d] %lf, %lf, %lf, %lf, %lf, %lf\n",
								objEO.m_vec_eo[idxIM].iImg,
								objEO.m_vec_eo[idxIM].dXc,
								objEO.m_vec_eo[idxIM].dYc,
								objEO.m_vec_eo[idxIM].dZc,
								objEO.m_vec_eo[idxIM].dOmega,
								objEO.m_vec_eo[idxIM].dPhi,
								objEO.m_vec_eo[idxIM].dKappa);
			}
			printf("\tThe adjusted GPs:\n");
			for (gpID = 0; gpID < n_total_GP; gpID++)
			{
				printf("\t[%d] %lf, %lf, %lf\n",
								objGP.m_vec_initGPID[gpID],
								objGP.m_vec_initGP2[gpID].X,
								objGP.m_vec_initGP2[gpID].Y,
								objGP.m_vec_initGP2[gpID].Z);
			}
	
		#endif	//__LSM_DEBUGGING


	// The function returns 'true' if the error is less than delta.
	return DBL_LT(d_error, m_dDelta)?true:false;
 }


/********************************************************************
 * Description:	Compute Ae, Ap and y matrices.
 ********************************************************************/
void CSeqSimAT::ComputeAy(	/*in*/	IP const &ip_f,
							/*in*/	IO const &io, 
							/*in*/	EO const &eo, 
							/*in*/	GP const &gp, 
							/*in*/	Matrix const &R,
							/*in*/	Matrix const &dROm,
							/*in*/	Matrix const &dRPh,
							/*in*/	Matrix const &dRKp,
							/*out*/	Matrix &Ae,
							/*out*/	Matrix &Ap,
							/*out*/	Matrix &y)
{
	// Variable declaration
	Matrix	GC(3,1), ND(3,1), F0(2,1), dND(3,9);
	Matrix	M, N;	// Temporary
	double dTemp;	// Temporary
	
	// GC = GP - PC
	GC.Set(0, 0, gp.X - eo.dXc);
	GC.Set(1, 0, gp.Y - eo.dYc);
	GC.Set(2, 0, gp.Z - eo.dZc);

	ND = R * GC;

	// F0 = IO(1:2) - IO(3) / ND(3) * ND(1:2,2)
	dTemp = io.dF / ND.Get(2, 0);
	F0.Set(0, 0, io.dPPX - (dTemp *  ND.Get(0, 0)));
	F0.Set(1, 0, io.dPPY - (dTemp *  ND.Get(1, 0)));

	// dND(:, 1:3) = -R
	dND.Set(0, 0, -R.Get(0, 0)); 
	dND.Set(1, 0, -R.Get(1, 0)); 
	dND.Set(2, 0, -R.Get(2, 0)); 

	dND.Set(0, 1, -R.Get(0, 1)); 
	dND.Set(1, 1, -R.Get(1, 1)); 
	dND.Set(2, 1, -R.Get(2, 1)); 

	dND.Set(0, 2, -R.Get(0, 2)); 
	dND.Set(1, 2, -R.Get(1, 2)); 
	dND.Set(2, 2, -R.Get(2, 2)); 

	// dND(:, 4) = dROm * GC
	M = dROm * GC;
	dND.Set(0, 3, M.Get(0, 0));
	dND.Set(1, 3, M.Get(1, 0));
	dND.Set(2, 3, M.Get(2, 0));

	// dND(:, 5) = dRPh * GC
	M = dRPh * GC;
	dND.Set(0, 4, M.Get(0, 0));
	dND.Set(1, 4, M.Get(1, 0));
	dND.Set(2, 4, M.Get(2, 0));

	// dND(:, 6) = dRKp * GC
	M = dRKp * GC;
	dND.Set(0, 5, M.Get(0, 0));
	dND.Set(1, 5, M.Get(1, 0));
	dND.Set(2, 5, M.Get(2, 0));

	// dND(:, 7:9) = R
	dND.Set(0, 6, R.Get(0, 0)); 
	dND.Set(1, 6, R.Get(1, 0)); 
	dND.Set(2, 6, R.Get(2, 0)); 

	dND.Set(0, 7, R.Get(0, 1)); 
	dND.Set(1, 7, R.Get(1, 1)); 
	dND.Set(2, 7, R.Get(2, 1)); 

	dND.Set(0, 8, R.Get(0, 2)); 
	dND.Set(1, 8, R.Get(1, 2)); 
	dND.Set(2, 8, R.Get(2, 2)); 

	// Ae = IO(3) / ND(3)^2 * [-ND(3) 0 ND(1); 0 -ND(3) ND(2)] * dND(:,1:6)
	M.SetDim(2, 3);		// For [-ND(3) 0 ND(1); 0 -ND(3) ND(2)]

	M.Set(0, 0, -ND.Get(2, 0));
	M.Set(0, 1, 0.0);
	M.Set(0, 2, ND.Get(0, 0));

	M.Set(1, 0, 0.0);
	M.Set(1, 1, -ND.Get(2, 0));
	M.Set(1, 2, ND.Get(1, 0));

	N.SetDim(3, 6);	// For dND(:,1:6)

	N.Set(0, 0, dND.Get(0, 0));
	N.Set(0, 1, dND.Get(0, 1));
	N.Set(0, 2, dND.Get(0, 2));
	N.Set(0, 3, dND.Get(0, 3));
	N.Set(0, 4, dND.Get(0, 4));
	N.Set(0, 5, dND.Get(0, 5));

	N.Set(1, 0, dND.Get(1, 0));
	N.Set(1, 1, dND.Get(1, 1));
	N.Set(1, 2, dND.Get(1, 2));
	N.Set(1, 3, dND.Get(1, 3));
	N.Set(1, 4, dND.Get(1, 4));
	N.Set(1, 5, dND.Get(1, 5));

	N.Set(2, 0, dND.Get(2, 0));
	N.Set(2, 1, dND.Get(2, 1));
	N.Set(2, 2, dND.Get(2, 2));
	N.Set(2, 3, dND.Get(2, 3));
	N.Set(2, 4, dND.Get(2, 4));
	N.Set(2, 5, dND.Get(2, 5));	

	// Ae = IO(3) / ND(3)^2 * [-ND(3) 0 ND(1); 0 -ND(3) ND(2)] * dND(:,1:6)
	dTemp = io.dF / (ND.Get(2, 0) * ND.Get(2, 0));
	Ae = M * N * dTemp;

	// Ap = IO(3) / ND(3)^2 * [-ND(3) 0 ND(1); 0 -ND(3) ND(2)] * dND(:,7:9)
	N.SetDim(3, 3);	// For dND(:,7:9)

	N.Set(0, 0, dND.Get(0, 6));
	N.Set(0, 1, dND.Get(0, 7));
	N.Set(0, 2, dND.Get(0, 8));

	N.Set(1, 0, dND.Get(1, 6));
	N.Set(1, 1, dND.Get(1, 7));
	N.Set(1, 2, dND.Get(1, 8));

	N.Set(2, 0, dND.Get(2, 6));
	N.Set(2, 1, dND.Get(2, 7));
	N.Set(2, 2, dND.Get(2, 8));
	
	Ap = M * N * dTemp;

	// y = IP - F0
	y.Set(0, 0, ip_f.dX - F0.Get(0, 0));
	y.Set(1, 0, ip_f.dY - F0.Get(1, 0));

}

/********************************************************************
 * Description:	Debug the AT results.
 ********************************************************************/
void CSeqSimAT::DebugSimATResults()
{
	int round;
	int	iGP, iIM;
	int	nGP, nIM;
	int id_start_img, id_last_img;
	int	idxGP;
	GP	gp_a, gp_i;
	EO	eo;

	// Open a file to write the summary results.
	FILE * fResult = _fsopen("Summary.txt", "w", _SH_DENYWR);

	if (fResult == NULL)
	{
		printf("CSeqSimAT::DebugSimATResults()\n");
		printf("\tFailed to create a resulting file\n");
		return;
	}

	// Iterate through each round result.
	for (round = 0; round < m_nRound; round++)
	{
		fprintf(fResult, "============================================================\n");
		fprintf(fResult, "                           AT Result                        \n", round);
		fprintf(fResult, "============================================================\n");

		// Report if the operation failed.
		if (m_vec_result[round] == false)
		{
			fprintf(fResult, "****** The operation failed. ******\n");		
			// continue;
		}

		// Get the number of image proceeded in this round.
		nIM = m_vec_objEO[round].m_nImages;

		// Obtain the starting Image ID.
		id_start_img	= 1;
		if (nIM > m_nmax_sim_img)
			id_start_img = nIM - m_nmax_sim_img + 1;

		// Obtain the last Image ID.
		id_last_img	= id_start_img + nIM - 1;

		fprintf(fResult, "The total number of images in computation: %d [ID%d ~ ID%d].\n", nIM, id_start_img, id_last_img);		
		fprintf(fResult, "The total number of running iterations: %d.\n", m_vec_iters[round]);
		fprintf(fResult, "The updated value with respect to the pre-defined threshold (%lf): %lf.\n", m_dDelta, m_vec_error[round]);
		fprintf(fResult, "The computational time for this round: %d milliseconds.\n", m_vec_elapsed_time[round]);

		// Print out the EO parameters.
		fprintf(fResult, "The adjusted EO parameters for each image:\n");
		fprintf(fResult, "\tImage\tXc\tYc\tZc\tOmega\tPhi\tKappa\n");

		for (iIM = 0; iIM < nIM; iIM++)
		{
			eo = m_vec_objEO[round].m_vec_eo[iIM];
		
			fprintf(fResult, "\t%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", 
							eo.iImg, eo.dXc, eo.dYc, eo.dZc, eo.dOmega, eo.dPhi, eo.dKappa);
			
		}

		// Print out the GP.
		fprintf(fResult, "The adjusted GPs compared with the initial GPs:\n");
		fprintf(fResult, "\tGP ID\tAdjusted GP [X,Y,Z]\t\tInitial GPs [X,Y,Z]\n");

		nGP = m_vec_objGP[round].m_vec_initGPID.size();
		for (iGP = 0; iGP < nGP; iGP++)
		{
			// Ground point index.
			idxGP = m_vec_objGP[round].m_vec_initGPID[iGP];

			// Adjusted GP by the Sim AT process.
			gp_a = m_vec_objGP[round].m_vec_initGP2[iGP];

			// Initial GP for comparison.
			gp_i = m_vec_objGPi[round].m_vec_initGP2[iGP];

			fprintf(fResult,	"\tGP %d\t[%lf,%lf,%lf]\t[%lf,%lf,%lf]\n",
								idxGP, gp_a.X, gp_a.Y, gp_a.Z, gp_i.X, gp_i.Y, gp_i.Z);
		}

		fprintf(fResult, "\n");
	}

	fclose(fResult);
}
