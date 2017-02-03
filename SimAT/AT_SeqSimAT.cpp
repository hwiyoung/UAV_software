/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "AT_Definition.h"
#include "AT_SeqSimAT.h"
#include "AT_Matrix.h"
#include "stdio.h"
#include <iostream>
#include <direct.h>   
#include <fstream>		
#include <strstream>
#include <vector>
#include <conio.h>
#include <math.h>
#include <iterator>

using namespace std;

/********************************************************************
 * Description:	Construction
 ********************************************************************/
AT_CSeqSimAT::AT_CSeqSimAT()
{}

/********************************************************************
 * Description:	Set the member parameters.
 ********************************************************************/
void AT_CSeqSimAT::Set(int nImages,int max_iter,int ninit_sim_img,
					double dDelta,double std_IP,double std_GPS,double std_INS)
{
	m_nImages		=	nImages;
	m_max_iter		=	max_iter;
	m_ninit_sim_img	=	ninit_sim_img;
	m_dDelta		=	dDelta;	
	m_std_IP		=	std_IP;
	m_std_GPS		=	std_GPS;
	m_std_INS		=	std_INS;
}

/********************************************************************
 * Description:	Estimate sequential AT simultaneous.
 ********************************************************************/
void AT_CSeqSimAT::SeqSimEstimation(AT_VEC_IP const &vec_ip_f, AT_VEC_INT const &vec_n_IP,
									AT_IO const &io, AT_VEC_EO const &vec_eo)
{
	int round;			// Round of iteration.
	int	n_run_img;		// Number of images used in processing.
	int	n_run_ip;		// Number of images points from Image#0 to current.
	DWORD start, end;	// Measuring computational time.
	int i;
	bool bSuccess;

	// Get the number of images points from the first image to the current image.
	n_run_ip = 0;

	// Sum up the number of image points.
	for (i = 1; i <= m_ninit_sim_img; i++)	// Image index starts at 1.
		n_run_ip = n_run_ip + vec_n_IP[i];

	// The round of processing.
	round = 0;

	// Perform simultaneous AT in which one image is added to each round of processing.
	for (n_run_img = m_ninit_sim_img+1; n_run_img <= m_nImages; n_run_img++, round++)//왜 이렇게 해놓은지 전혀 모르겠음
	{
		printf("====================================================\n");
		printf("          Simultaneous AT processing with %d images.\n", m_nImages);
		printf("====================================================\n");
		
		// Get the number of image points until this current image. 
		n_run_ip = n_run_ip + vec_n_IP[n_run_img];

		// Performing simultaneous AT.
		Matrix	Qee, Qep, Qpp;
		int		n_iterations = 0;
		double	d_error = -1.0;
		AT_CGP						objGP;	
		AT_CExteriorOrientation	objEO(vec_eo, n_run_img, false);
		AT_CRotationalMatrix		objRot;
		objRot.Compute3DRotMatrix(objEO.m_vec_eo, objEO.m_nImages);
		objRot.ComputeDerivative3DRotMatrix(objEO.m_vec_eo, objEO.m_nImages);

		// Get the tick of starting time.
		start = GetTickCount();

		// Perform Simultaneous AT
		bSuccess = SimultaneousAT(	/*in*/	n_run_img,
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

		// Get the tick of end time.
		end = GetTickCount();

		// Store the flag if the Sim AT operation success or not.
		//m_vec_result.push_back(bSuccess);

		//// Store the adjusted EO (CExteriorOrientation) into the vector.
		//m_vec_objEO.push_back(objEO);

		//// Store the adjusted GP (CGP object) into the vector.
		//m_vec_objGP.push_back(objGP);

		//// Store the no. iterations for this round to the vector.
		//m_vec_iters.push_back(n_iterations);

		//// Store the error wrt the pre-defined threshold to the vector.
		//m_vec_error.push_back(d_error);

		// Store the elapsed time into the vector.
		DWORD elapsed_time = end - start;
		m_vec_elapsed_time.push_back(elapsed_time);

		// Compute GP init for comparison.
		AT_CGP objGPi;
		objGPi.MakeGPList(vec_ip_f, n_run_img);		
		objGPi.GPInitApproximation(vec_ip_f, io, vec_eo, n_run_img, objRot, m_max_iter, m_dDelta, true);

		// Store the GP from initial approximation for comparison.
		m_vec_objGPi.push_back(objGPi);

		printf("AT의 처리 시간은 %d초 입니다\n\n", elapsed_time/1000);
	}

	m_nRound = round;
}

/********************************************************************
 * Description:	Perform Simultaneous AT tasks.
 ********************************************************************/
 bool AT_CSeqSimAT::SimultaneousAT(/*in*/	int n_run_img,
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
								/*out*/ double &d_error)
 {
	// Variable declaration
	Matrix	Nee, Nep, Npp, Ce, Cp, iNpp, Q;		// Normal matrices.
	int r, c, k, n_IP, n_total_GP, n_IM;
	int	imgID, gpID, idxGP;
	AT_VEC_IP vec_ip_f2;
	AT_VEC_IP vec_ip_f3;

	// Make the list of distinct ground point ID
	objGP.MakeGPList(vec_ip_f, n_run_img);

	// Initial approximation for GPs
	bool bFineApprox = true;
	objGP.GPInitApproximation(vec_ip_f, io, vec_eo, n_run_img, objRot, m_max_iter, m_dDelta, bFineApprox);

	/////////////////////////////////////////////////////////////////////////////////////////////////////////
	copy(vec_ip_f.begin(), vec_ip_f.end(), std::back_inserter(vec_ip_f2)); // IP가 const_iterator라서 수정이 안되므로 변경하여 작성

	// Erase GP 
	double sum = 0;
	double sum_devsq = 0;
	
	for (int i = 0; i < objGP.m_vec_initGP2.size(); i++)
	{
		sum = sum + objGP.m_vec_initGP2[i].Z;
	}
	double avg = sum / objGP.m_vec_initGP2.size();

	for (int j = 0; j < objGP.m_vec_initGP2.size(); j++)
	{
		sum_devsq = sum_devsq + pow((objGP.m_vec_initGP2[j].Z - avg),2);
	}
	double avg_devsq = sum_devsq / (objGP.m_vec_initGP2.size()-1);
	double std_Z = sqrt(avg_devsq);

	double high_Z = avg + std_Z;
	double low_Z = avg - std_Z;
	
	//GP 제거
	int a, k1, t;
	a = 0;
	for (vector<AT_GP>::iterator iter = objGP.m_vec_initGP2.begin(); iter != objGP.m_vec_initGP2.end(); ++iter) 
	{
		if (iter->Z < low_Z)
		{
			iter = objGP.m_vec_initGP2.erase(iter);
			//cout << a << endl;
			k1 = 0;
			for (vector<int>::iterator iter2 = objGP.m_vec_initGPID.begin(); iter2 != objGP.m_vec_initGPID.end(); ++iter2)
			{
				if (k1 == a)
				{
					iter2 = objGP.m_vec_initGPID.erase(iter2);
					a--;
					break;
				}
				k1++;
			}
		}
		else if (iter->Z > high_Z)
		{
			iter = objGP.m_vec_initGP2.erase(iter);
			t = 0;
			for (vector<int>::iterator iter3 = objGP.m_vec_initGPID.begin(); iter3 != objGP.m_vec_initGPID.end(); ++iter3)
			{
				if (t == a)
				{
					iter3 = objGP.m_vec_initGPID.erase(iter3);
					a--;
					break;
				}
				t++;
			}
		}
		a++;
	}

	int b;
	int asdf;
	asdf = 0;
	for (vector<AT_IP>::iterator iter4 = vec_ip_f2.begin(); iter4 != vec_ip_f2.end(); ++iter4)
	{		
		b = 0;
		for (vector<int>::iterator iter5 = objGP.m_vec_initGPID.begin(); iter5 != objGP.m_vec_initGPID.end(); ++iter5)
		{
			if (iter4->iObj == objGP.m_vec_initGPID[b])
			{
				vec_ip_f3.push_back(vec_ip_f2[asdf]);
			}						
			b++;				
		}		
		asdf++;
	}

	//FILE * ip_result;
	//FILE * gp_result;

	//

	//if ((ip_result = _fsopen("ip_result.txt", "w", _SH_DENYWR)) == NULL)
	//{
	//	cout << "지상점의 초기값 텍스트 파일을 생성할 수  없습니다" << endl;
	//	cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
	//	_getch();
	//	exit(0);
	//}

	//if ((gp_result = _fsopen("gp_result.txt", "w", _SH_DENYWR)) == NULL)
	//{
	//	cout << "지상점의 초기값 텍스트 파일을 생성할 수  없습니다" << endl;
	//	cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
	//	_getch();
	//	exit(0);
	//}

	//for (int ip_i = 0; ip_i < vec_ip_f3.size(); ip_i++)
	//{
	//	fprintf(ip_result, "%d\t%d\t%f\t%f\n", vec_ip_f3[ip_i].iImg, vec_ip_f3[ip_i].iObj, vec_ip_f3[ip_i].dX, vec_ip_f3[ip_i].dY);
	//}

	//for (int gp_i = 0; gp_i < objGP.m_vec_initGPID.size(); gp_i++)
	//{
	//	fprintf(gp_result, "%d\t%lf\t%lf\t%lf\n", objGP.m_vec_initGPID[gp_i], objGP.m_vec_initGP2[gp_i].X, objGP.m_vec_initGP2[gp_i].Y, objGP.m_vec_initGP2[gp_i].Z);
	//}
	//
	//fclose(ip_result);
	//fclose(gp_result);

	//_getch();

	/////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Debug out the inital approximated GPs.
	objGP.DebugInitGP(false);

	// Get the number of GPs that exists more than 2 images.

	n_total_GP = objGP.m_vec_initGP2.size();
	//n_total_GP = objGP.m_n_distinct_GP;

	if (n_total_GP == 0)
	{
		printf("Operation Failed!\n");
		printf("Not a single pair of tie points existed.\n");
		
		return false;
	}	

	// Iterate until the error is smaller than the delta threshold or maximum iterations is reached.
	for (k = 0; k < m_max_iter; k++)
	{
		printf("\tIteration %d\n", k+1);

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

		// Initialize the rotational matrix
		AT_CRotationalMatrix rot;
		rot.Compute3DRotMatrix(objEO.m_vec_eo, objEO.m_nImages);
		rot.ComputeDerivative3DRotMatrix(objEO.m_vec_eo, objEO.m_nImages);

		// ***** Compute the normal matrix ***** //
		// Iterate through each image point until the last point in of the current image.
		//cout << "n_run_ip : "<< n_run_ip<<endl;
		int nip=0;
		int n_run_ip2 = vec_ip_f3.size();
		//for (n_IP = 0; n_IP < 1; n_IP++)
		for (n_IP = 0; n_IP < n_run_ip2; n_IP++)
		{
			//cout << "n_IP = "<< n_IP+1 << endl;
			Matrix	Ae, Ap, y(2,1);
			Matrix	M;			// Temporary
			double	d;			// Temporary
			AT_EO		eo;			// Temporaryo
			AT_GP		gp;			// Temporary

			gpID	= vec_ip_f3[n_IP].iObj;			// GP ID
			
			// Continue if this GP ID exists less than two images.
			objGP.GetNumImagesGPExisted(gpID);
			if (objGP.GetNumImagesGPExisted(gpID) < 2)
				continue;
			
			imgID	= vec_ip_f3[n_IP].iImg;				// Image ID	
			
			int i;
			for ( i = 0; i < n_total_GP; i++)
			{			
				if(objGP.m_vec_initGPID[i] == gpID)
					idxGP=i;
					continue;		
			}	
			
			 
			gp		= objGP.m_vec_initGP2[idxGP];		// GP coordinates.			
			eo		= objEO.GetEO(imgID);				// EO of Image ID		
			

			//cout << "n_IP = " << n_IP << endl;
			// Compute the Ae, Ag and y matrices.
			ComputeAy(	/*in*/	vec_ip_f3[n_IP] /*vec_ip_f[nip]*/,
						/*in*/	io,
						/*in*/	eo,
						/*in*/	gp, 
						/*in*/	rot.m_vec_RotMatrix[imgID],
						/*in*/	rot.m_vec_om_dRotMatrix[imgID],
						/*in*/	rot.m_vec_ph_dRotMatrix[imgID],
						/*in*/	rot.m_vec_kp_dRotMatrix[imgID],
						/*out*/	Ae, /*out*/	Ap, /*out*/	y);
			//cout << Ae << endl;
			//cout << endl;
			//cout << "Compute Ae, Ag, y matrix end" << endl;
			// Nee(imi*6-5:imi*6,imi*6-5:imi*6) += Ae_i' * Ae_i / std_IP ^ 2
						
			M = Ae.Trans() * Ae * (1.0 / m_std_IP / m_std_IP);

			//ofstream M_res("Results\\AT_Results\\M_res.txt");

			//for(int i =0 ; i<6 ; i++) {
			//	for(int j=0; j<6 ; j++) {
			//		M_res << M[i+j] << "\t";
			//	} M_res << "\n";
			//}
			//M_res.close();

			for (r = 0; r < N_UNKNOWN_EO; r++)
			{
				for (c = 0; c < N_UNKNOWN_EO; c++)
				{
					d = Nee.Get((imgID-1)*N_UNKNOWN_EO + r, (imgID-1)*N_UNKNOWN_EO + c) + M.Get(r, c);
					Nee.Set((imgID-1)*N_UNKNOWN_EO + r, (imgID-1)*N_UNKNOWN_EO + c, d);
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
			//cout << " Ngg size = " << Npp.size() << endl;	

			// Neg(imi*6-5:imi*6,gpi*3-2:gpi*3) = Ae_i' * Ap_i / std_IP ^ 2;
			M = Ae.Trans() * Ap * (1.0 / m_std_IP / m_std_IP);
			
			for (r = 0; r < N_UNKNOWN_EO; r++)
			{
				for (c = 0; c < N_UNKNOWN_GP; c++)
				{
					Nep.Set((imgID-1)*N_UNKNOWN_EO + r, idxGP*N_UNKNOWN_GP + c, M.Get(r, c));
				}
			}			

			//cout << "idxGP : " << idxGP << endl;
			//cout << "imgID : " << imgID << endl << endl;

			// Ce(imi*6-5:imi*6,1) += Ae_i' * y_i / std_IP ^ 2;
			M = Ae.Trans() * y * (1.0 / m_std_IP / m_std_IP);
			
	
			for (r = 0; r < N_UNKNOWN_EO; r++)
			{
				d = Ce.Get((imgID-1)*N_UNKNOWN_EO + r, 0) + M.Get(r, 0);
			
				Ce.Set((imgID-1)*N_UNKNOWN_EO + r, 0, d);
			}
			
			// Cg(gpi*3-2:gpi*3,1) += Ap_i' * y_i / std_IP ^ 2;
			M = Ap.Trans() * y * (1.0 / m_std_IP / m_std_IP);

			for (r = 0; r < N_UNKNOWN_GP; r++)
			{
				d = Cp.Get(idxGP*N_UNKNOWN_GP + r, 0) + M.Get(r, 0);
				Cp.Set(idxGP*N_UNKNOWN_GP + r, 0, d);
			}	
			//cout << y << endl;
			//nip=nip+1;

		}

		// ***** Consider constraints ***** //

		// Iterate through each image (Image#ID starting from 1)
		for (n_IM = 1; n_IM <= n_run_img; n_IM++)
		{
			Matrix	M, N, ye(6,1);	// Temporary
			double	d;				// Temporary
			AT_EO		eo_m, eo_adj;

			// Nee(ni*6-5:ni*6-3,ni*6-5:ni*6-3) += eye(3) / std_GPS ^ 2;
			M = N.Eye(3) * (1.0 / m_std_GPS / m_std_GPS);
	
			for (r = 0; r < N_UNKNOWN_GP; r++)
			{
				for (c = 0; c < N_UNKNOWN_GP; c++)
				{	
					d = Nee.Get((n_IM-1)*N_UNKNOWN_EO + r, (n_IM-1)*N_UNKNOWN_EO + c) + M.Get(r, c);
					//cout << M.Get(r,c) << endl;
					Nee.Set((n_IM-1)*N_UNKNOWN_EO + r, (n_IM-1)*N_UNKNOWN_EO + c, d);
				}
			}

			// Nee(ni*6-2:ni*6,ni*6-2:ni*6) += eye(3) / std_INS ^ 2;
			M = N.Eye(3) * (1.0 / m_std_INS / m_std_INS);

			for (r = 3; r < N_UNKNOWN_EO; r++)
			{
				for (c = 3; c < N_UNKNOWN_EO; c++)
				{
					d = Nee.Get((n_IM-1)*N_UNKNOWN_EO + r, (n_IM-1)*N_UNKNOWN_EO + c) + M.Get(r-N_UNKNOWN_GP, c-N_UNKNOWN_GP);
					Nee.Set((n_IM-1)*N_UNKNOWN_EO + r, (n_IM-1)*N_UNKNOWN_EO + c, d);
				}
			}

			// ye_i = EO_c(ni,:)' - Kt_e(ni*6-5:ni*6,1);
			eo_m	= objGP.GetEO(vec_eo, n_IM);		// EO from the direct measurement.
			eo_adj	= objEO.GetEO(n_IM);				// EO from approximation process.	

			ye.Set(0, 0, eo_m.dXc - eo_adj.dXc);
			ye.Set(1, 0, eo_m.dYc - eo_adj.dYc);
			ye.Set(2, 0, eo_m.dZc - eo_adj.dZc);
			ye.Set(3, 0, eo_m.dOmega - eo_adj.dOmega);
			ye.Set(4, 0, eo_m.dPhi - eo_adj.dPhi);
			ye.Set(5, 0, eo_m.dKappa - eo_adj.dKappa);

			// Ce(ni*6-5:ni*6-3,1) += ye_i(1:3,1) / std_GPS ^ 2;
			for (r = 0; r < N_UNKNOWN_GP; r++)
			{
				d = Ce.Get((n_IM-1)*N_UNKNOWN_EO + r, 0) + ye.Get(r, 0)/m_std_GPS/m_std_GPS;
				Ce.Set((n_IM-1)*N_UNKNOWN_EO + r, 0, d);
			}

			// Ce(ni*6-2:ni*6,1) += ye_i(4:6,1) / std_INS ^ 2;
			for (r = 3; r < N_UNKNOWN_EO; r++)
			{
				d = Ce.Get((n_IM-1)*N_UNKNOWN_EO + r, 0) + ye.Get(r, 0)/m_std_INS/m_std_INS;
				Ce.Set((n_IM-1)*N_UNKNOWN_EO + r, 0, d);
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

		cout << "Multiply" << endl;
		// Q1 = Nep * iNpp;
		Q = Nep * iNpp;
		cout << "Multiply finish" << endl;

		// Nr = ( Nee - Q1 * Nep' );
		Nr = Nee - Q * Nep.Trans();

		// Qee = inv(Nr)
		Qee = Nr.Inv();

		// kt_e = Qee * ( Ce - Q1 * Cp );
		kt_e = Qee * (Ce - Q * Cp);
		//cout << kt_e << endl;
		// kt_g = iNpp * Cp - Q1' * kt_e;
		kt_g = iNpp * Cp - Q.Trans() * kt_e;

		// Update the adjusted EOs.		
		for (n_IM = 1; n_IM <= n_run_img; n_IM++)
		{
			int idxEO = objEO.GetEOIndex(n_IM);

			// Kt_e = Kt_e + kt_e		
			objEO.m_vec_eo[idxEO].dXc		+= kt_e.Get((n_IM-1)*N_UNKNOWN_EO + 0, 0);
			objEO.m_vec_eo[idxEO].dYc		+= kt_e.Get((n_IM-1)*N_UNKNOWN_EO + 1, 0);
			objEO.m_vec_eo[idxEO].dZc		+= kt_e.Get((n_IM-1)*N_UNKNOWN_EO + 2, 0);
			objEO.m_vec_eo[idxEO].dOmega	+= kt_e.Get((n_IM-1)*N_UNKNOWN_EO + 3, 0);
			objEO.m_vec_eo[idxEO].dPhi		+= kt_e.Get((n_IM-1)*N_UNKNOWN_EO + 4, 0);
			objEO.m_vec_eo[idxEO].dKappa	+= kt_e.Get((n_IM-1)*N_UNKNOWN_EO + 5, 0);	
		}

		// Update the adjusted GPs.		
		for (gpID = 0; gpID < n_total_GP; gpID++)
		{
			// Kt_p = Kt_p + kt_g
			objGP.m_vec_initGP2[gpID].X += kt_g.Get(gpID*N_UNKNOWN_GP + 0, 0);
			objGP.m_vec_initGP2[gpID].Y += kt_g.Get(gpID*N_UNKNOWN_GP + 1, 0);
			objGP.m_vec_initGP2[gpID].Z += kt_g.Get(gpID*N_UNKNOWN_GP + 2, 0);
		}

		// if norm(kt_g)/no_GP/3 < delta
	
		//d_error = kt_g.Norm(0, false)/objGP.m_n_distinct_GP/3.0;
		d_error = kt_g.Norm(0, false) / 1121 / 3.0;


		cout << "Error = " << d_error << endl;

		if (  d_error > 10000)
		{
			m_vec_result.push_back(false);
			objEO.m_vec_eo = vec_eo;
			m_vec_objEO.push_back(objEO);
			m_vec_objGP.push_back(objGP);
			m_vec_iters.push_back(n_iterations);
			m_vec_error.push_back(d_error);
			break;
		}

		if ( d_error < 0.001 )
		{
			m_vec_result.push_back(true);
			m_vec_objEO.push_back(objEO);
			m_vec_objGP.push_back(objGP);
			m_vec_iters.push_back(n_iterations);
			m_vec_error.push_back(d_error);
			break;

			//FILE * N_result;			

			//if ((N_result = _fsopen("N_result.txt", "w", _SH_DENYWR)) == NULL)
			//{
			//	cout << "지상점의 초기값 텍스트 파일을 생성할 수  없습니다" << endl;
			//	cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
			//	_getch();
			//	exit(0);
			//}			

			//for (int ip_i = 0; ip_i < vec_ip_f3.size(); ip_i++)
			//{
			//	fprintf(N_result, "%d\t%d\t%f\t%f\n", vec_ip_f3[ip_i].iImg, vec_ip_f3[ip_i].iObj, vec_ip_f3[ip_i].dX, vec_ip_f3[ip_i].dY);
			//}
			//
			//fclose(N_result);
		}

	}

	// no_IT = k;
	n_iterations = k;

	// Qep = -Qee * Q1;
	Qep = (Qee * (-1.0)) * Q;

	// Qpp = iNpp + Q1' * Qee * Q1;
	Qpp = iNpp + (Q.Trans()* Qee * Q);
 	
	return true;

 }


/********************************************************************
 * Description:	Compute Ae, Ag and y matrices.
 ********************************************************************/
void AT_CSeqSimAT::ComputeAy(	/*in*/	AT_IP const &ip_f,
							/*in*/	AT_IO const &io, 
							/*in*/	AT_EO const &eo, 
							/*in*/	AT_GP const &gp, 
							/*in*/	Matrix const &R,
							/*in*/	Matrix const &dROm,
							/*in*/	Matrix const &dRPh,
							/*in*/	Matrix const &dRKp,
							/*out*/	Matrix &Ae,
							/*out*/	Matrix &Ap,
							/*out*/	Matrix &y)
{
	// Variable declaration
	//cout << endl;
	//cout << "Compute Ae, Ag, and y matrix" << endl;
	Matrix	GC(3,1), ND(3,1), F0(2,1), dND(3,9);
	Matrix	dRe_om(3,3), dRe_ph(3,3), dRe_kp(3,3);
	//Matrix	R1t, R2t;
	Matrix	M, N;	// Temporary
	double dTemp;	// Temporary	
	double gpx, gpy, gpz;
	gpx=gp.X;
	gpy=gp.Y;
	gpz=gp.Z;

	// GC = GP - PC
/*	GC.Set(0, 0, gp.X - eo.dXc);
	GC.Set(1, 0, gp.Y - eo.dYc);
	GC.Set(2, 0, gp.Z - eo.dZc);*/
	GC.Set(0, 0, gpx - eo.dXc);
	GC.Set(1, 0, gpy - eo.dYc);
	GC.Set(2, 0, gpz - eo.dZc);
	
	//R1t=R.Trans();
	ND = R * GC;

	// F0 = IO(1:2) - IO(3) / ND(3) * ND(1:2,2)
	dTemp = io.dF / ND.Get(2, 0);
	//F0.Set(0, 0, /*io.dPPX*/ - (dTemp *  ND.Get(0, 0)));
	//F0.Set(1, 0, /*io.dPPY*/ - (dTemp *  ND.Get(1, 0)));

	F0.Set(0, 0, io.dPPX - (dTemp *  ND.Get(0, 0)));
	F0.Set(1, 0, io.dPPY - (dTemp *  ND.Get(1, 0)));

	//cout << "F0 = " << F0 << endl;
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
	//cout << -R.Get(2,2) << endl;

	// dND(:, 4) = dROm * GC
	//dRe_om=dROm.Trans();
	M = dROm * GC;
	dND.Set(0, 3, M.Get(0, 0));
	dND.Set(1, 3, M.Get(1, 0));
	dND.Set(2, 3, M.Get(2, 0));
	
	// dND(:, 5) = dRPh * GC
	//dRe_ph=dRPh.Trans();
	M = dRPh * GC;
	dND.Set(0, 4, M.Get(0, 0));
	dND.Set(1, 4, M.Get(1, 0));
	dND.Set(2, 4, M.Get(2, 0));

	// dND(:, 6) = dRKp * GC
	//dRe_kp=dRKp.Trans();
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

	//cout << "dND = " << dND << endl;

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
	
	//cout << "N = " << endl;
		//cout << N << endl;
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

	//cout << "GC" << endl;
	//cout << GC << endl;
	//cout << "ND" << endl;
	//cout << ND << endl;
	//cout << "F0" << endl;
	//cout << F0 << endl;
	//cout << "dND" << endl;
	//cout << dND << endl;
	//cout << "Ae" << endl;
	//cout << Ae << endl;
	//cout << "Ap" << endl;
	//cout << Ap << endl;
	//cout << "y" << endl;
	//cout << y << endl;
}

void AT_CSeqSimAT::GetLine(istream &in, char *buf, const int bufSize) const
{
	while(in.good())
	{
		in.getline(buf, bufSize);
		if(buf[0]!='%') break;
	}
}

/********************************************************************
 * Description:	Debug the AT results.
 ********************************************************************/
void AT_CSeqSimAT::DebugSimATResults()
{
	int round;
	int	iGP, iIM;
	int	nGP, nIM;
	int	idxGP;
	AT_GP	gp_a, gp_i;
	AT_EO	eo;
	int	TOTAL_IMAGE;
	string imagename[200];

	const int BUFFER_SIZE = 200;
	char buffer[BUFFER_SIZE+1] = {0};
		
	ifstream fileInput("Input_80_170110\\AerialConfig.txt");

	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf_s(buffer, "%d", &TOTAL_IMAGE);

	GetLine(fileInput, buffer, BUFFER_SIZE);
	GetLine(fileInput, buffer, BUFFER_SIZE);

	for (int i=0; i<TOTAL_IMAGE; i++)
	{
		string img_file;

		GetLine(fileInput, buffer, BUFFER_SIZE);
		istrstream(buffer) >> img_file;	
		imagename[i]=img_file;
	}

	fileInput.close();


	char strFolderPath[] = { "Results\\AT_Results" };
     
    int nResult = _mkdir( strFolderPath );
 
    if( nResult == 0 )
    {
        printf( "폴더 생성 성공" );
    }
    else if( nResult == -1 )
	{
        //printf( "error no : %d", errno );
    }  

	// Open a file to write the summary results.
	FILE * fResult;
	FILE * GPi_Result;
	FILE * GPe_Result;


	ofstream EO_Result("Results\\AT_Results\\Estimated_EO.txt");


	if ((fResult = _fsopen("Results\\AT_Results\\Summary.txt", "w", _SH_DENYWR)) == NULL)
	{
		cout <<"Summary 파일을 생성할 수  없습니다" << endl;
		cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
		_getch();			
		exit(0);
	}

	if ((GPi_Result = _fsopen("Results\\AT_Results\\GP_i.txt", "w", _SH_DENYWR)) == NULL)
	{
		cout <<"지상점의 초기값 텍스트 파일을 생성할 수  없습니다" << endl;
		cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
		_getch();			
		exit(0);
	}

	if ((GPe_Result = _fsopen("Results\\AT_Results\\GP_e.txt", "w", _SH_DENYWR)) == NULL)
	{
		cout <<"지상점의 조정된 값의 텍스트 파일을 생성할 수  없습니다" << endl;
		cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
		_getch();			
		exit(0);
	}
	
	for (round = 0; round < m_nRound; round++)
	{
		fprintf(fResult, "============================================================\n");
		fprintf(fResult, "                           AT_Summary                       \n");
		fprintf(fResult, "============================================================\n");

		fprintf(fResult, "The total number of images in computation: %d\n", m_ninit_sim_img + round + 1);		
		fprintf(fResult, "The total number of running iterations: %d\n", m_vec_iters[round]);
		fprintf(fResult, "The error with respect to the pre-defined threshold (%lf): %lf\n", m_dDelta, m_vec_error[round]);
		fprintf(fResult, "The computational time for this round: %d seconds\n", m_vec_elapsed_time[round]/1000);

		fprintf(fResult, "The adjusted EO parameters for each image:\n");
		fprintf(fResult, "Image\tXc\tYc\tZc\tOmega\tPhi\tKappa\n");

		nIM = m_vec_objEO[round].m_nImages;
		for (iIM = 0; iIM < nIM; iIM++)
		{
			eo = m_vec_objEO[round].m_vec_eo[iIM];	
			fprintf(fResult, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", eo.iImg, eo.dXc, eo.dYc, eo.dZc, eo.dOmega, eo.dPhi, eo.dKappa);
			EO_Result << imagename[iIM] << "\t"; 
			EO_Result.setf(ios_base::fixed, ios_base::floatfield);
			EO_Result << eo.dXc << "\t" << eo.dYc << "\t"<< eo.dZc << "\t";
			EO_Result << eo.dOmega << "\t" << eo.dPhi << "\t" << eo.dKappa << endl;
		}

		fprintf(GPi_Result, "\tGP ID\tInitial GP [X,Y,Z]\n");
		nGP = m_vec_objGP[round].m_vec_initGPID.size();
		for (iGP = 0; iGP < nGP; iGP++)
		{
			idxGP = m_vec_objGP[round].m_vec_initGPID[iGP];
			gp_i = m_vec_objGPi[round].m_vec_initGP2[iGP];

			fprintf(GPi_Result,	"%d\t%f\t%f\t%f\n",idxGP, gp_i.X, gp_i.Y, gp_i.Z);
		}

		fprintf(GPe_Result, "\tGP ID\tAdjusted GP [X,Y,Z]\n");
		for (iGP = 0; iGP < nGP; iGP++)
		{
			idxGP = m_vec_objGP[round].m_vec_initGPID[iGP];
			gp_a = m_vec_objGP[round].m_vec_initGP2[iGP];

			fprintf(GPe_Result,	"%d\t%lf\t%lf\t%lf\n",idxGP, gp_a.X, gp_a.Y, gp_a.Z);
		}
	}
	
	EO_Result.close();
	fclose(fResult);
	fclose(GPi_Result);
	fclose(GPe_Result);
}
