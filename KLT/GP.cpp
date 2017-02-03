/********************************************************
 * Project:			Automated AT
 * Last updated:	26 August 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#include <iterator>
#include "Definition.h"
#include "GP.h"
#include "Matrix.h"
#include "stdio.h"
#include <iostream>
#include <CRTDBG.H>	// To use _ASSERT

/********************************************************************
 * Description:	Construction
 ********************************************************************/
CGP::CGP()
{
	m_n_distinct_GP = 0;		// Number of distinct ground points
}

/********************************************************************
 * Description:	Constructor with member variables
 ********************************************************************/
CGP::CGP(VEC_INT const &vec_n_IP)
{
	// Store the (original) number of tie points exist for each image.
	copy(vec_n_IP.begin(), vec_n_IP.end(), std::back_inserter(m_vec_n_validGP));
}

/********************************************************************
 * Description:	Set the member variables
 ********************************************************************/
void CGP::Set(VEC_INT const &vec_n_IP)
{
	// Store the (original) number of tie points exist for each image.
	copy(vec_n_IP.begin(), vec_n_IP.end(), std::back_inserter(m_vec_n_validGP));
}

/********************************************************************
 * Description:	Reset all member variables
 ********************************************************************/
void CGP::ResetMember(VEC_INT const &vec_n_IP)
{
	m_n_distinct_GP = 0;		
	m_vec_distinct_GPID.clear();
	m_vec_nimg_GP.clear();
	m_vec_initGP.clear();
	m_vec_initGP2.clear();
	m_vec_initGPID.clear();
	m_vec_n_validGP.clear();
	m_vec_imagebits.clear();	
	m_base_img_id = 0;
	m_base_gp_id = 0;

	Set(vec_n_IP);
}

/********************************************************************
 * Description:	Clone all member variables.
 ********************************************************************/
CGP CGP::Clone(bool bAll)
{
	CGP obj;

	obj.m_n_distinct_GP = m_n_distinct_GP;
	
	if (m_vec_distinct_GPID.size() != 0)
		copy(m_vec_distinct_GPID.begin(), m_vec_distinct_GPID.end(), std::back_inserter(obj.m_vec_distinct_GPID));

	if (m_vec_nimg_GP.size() != 0)	
		copy(m_vec_nimg_GP.begin(), m_vec_nimg_GP.end(), std::back_inserter(obj.m_vec_nimg_GP));

	if (m_vec_initGP.size() != 0)	
		copy(m_vec_initGP.begin(), m_vec_initGP.end(), std::back_inserter(obj.m_vec_initGP));

	if (m_vec_initGPID.size() != 0)	
		copy(m_vec_initGPID.begin(), m_vec_initGPID.end(), std::back_inserter(obj.m_vec_initGPID));

	if (m_vec_imagebits.size() != 0)	
		copy(m_vec_imagebits.begin(), m_vec_imagebits.end(), std::back_inserter(obj.m_vec_imagebits));

	if (m_vec_n_validGP.size() != 0)
		copy(m_vec_n_validGP.begin(), m_vec_n_validGP.end(), std::back_inserter(obj.m_vec_n_validGP));

	if ( (bAll == true) && (m_vec_initGP2.size() != 0) )	
		copy(m_vec_initGP2.begin(), m_vec_initGP2.end(), std::back_inserter(obj.m_vec_initGP2));

	return obj;
}

/********************************************************************
 * Description:	Make all elements in m_vec_initGP1 and m_vec_initGP2
 *				to be equal according to the passed-in flag.
 * bRough = T	Both vector will store m_vec_initGP2
 * bRough = F	Both vector will store m_vec_initGP1
 ********************************************************************/
void CGP::ReplicateGP(bool bRough)
{
	if (false == bRough) // Set all elements as the 'rough' initial GP.
	{
		m_vec_initGP2.clear();
		copy(m_vec_initGP.begin(), m_vec_initGP.end(), std::back_inserter(m_vec_initGP2));
	}
	else // Set all elements as the 'fine' initial GP.
	{
		m_vec_initGP.clear();
		copy(m_vec_initGP2.begin(), m_vec_initGP2.end(), std::back_inserter(m_vec_initGP));
	}
}

/********************************************************************
 * Description:	Make the list of distinct ground point ID
 ********************************************************************/
void CGP::MakeGPList(	VEC_IP const &vec_ip_i,
						int idx_start_ip,
						int idx_last_ip,
						int start_img_id,
						int last_img_id,
						int nimages)
{
	IP		ip;
	int		i;
	SET_INT set_ip;	// To store the distince object ID
	SET_INT::iterator iter;

	// Set the based Image ID.
	m_base_img_id = start_img_id;

	// Set the based GP ID.
	m_base_gp_id = vec_ip_i[idx_start_ip].iObj;
	
	// Iterate through each point in the vector.
	for (i = idx_start_ip; i <= idx_last_ip; i++)
	{
		ip = vec_ip_i[i];

		// Skip if the image index is not in the image range.
		if ( (ip.iImg < start_img_id) || (ip.iImg > last_img_id) )
			continue;
		
		// Put the GP ID into the 'set' to produce distinct output.
		set_ip.insert(ip.iObj);

		// Fill the 'image bits' into the m_vec_imagebits vector.
		FillImageBits(ip);
	}

	// Get the number of distinct ground point.
	m_n_distinct_GP = set_ip.size();

	// Assign the distinct ground point into the vector list.
	for (iter = set_ip.begin(); iter != set_ip.end(); iter++)
	{
		m_vec_distinct_GPID.push_back(*iter);	
	}

}

/********************************************************************
 * Description:	Debug out the distinct object ID in the m_vec_distinct_GPID
 ********************************************************************/
void CGP::DebugDistinctGP()
{
	printf("Function: CGP::DebugDistinctGP()\n");
	printf("No. of distinct GP:\t%d\n", m_n_distinct_GP);

	for (int i = 0; i < m_n_distinct_GP; i++)
		printf("%d\n", m_vec_distinct_GPID[i]);

}

/********************************************************************
 * Description:	Debug out the intial approximation of GP.
 ********************************************************************/
void CGP::DebugInitGP(bool bRoughGP)
{
	VEC_GP * p_vec_GP;
	int i, nsize;

	// Set the pointer to the specified vector.
	if (bRoughGP == true)
		p_vec_GP = &m_vec_initGP;
	else
		p_vec_GP = &m_vec_initGP2;

	// Get the size of the GP.
	nsize = p_vec_GP->size();

	printf("Function: CGP::DebugInitGP()\n");
	printf("No. of initial approximated GP:\t%d\n", nsize);

	for (i = 0; i < nsize; i++)
		printf("[%d] : %lf, %lf, %lf\n", m_vec_initGPID[i], 
										 (*p_vec_GP)[i].X,
										 (*p_vec_GP)[i].Y, 
										 (*p_vec_GP)[i].Z);

}

/********************************************************************
 * Description:	Compute the 'rough' initial approximation of GP.
 *				This technique computes the ground point by considering
 *				the pair of tie points that have the largest distance
 *				between baseline.
 *				The approximated GPs stored in 'm_vec_initGP'.
 ********************************************************************/
void CGP::GPInitApproximation(	VEC_IP const &vec_ip_f, 
								int idx_start_ip,
								int idx_last_ip,
								IO const &io, 
								VEC_EO const &vec_eo,
								int start_img_id,
								int last_img_id,
								int const nimages,
								CRotationalMatrix const &rot,
								int const n_max_iter)
{
	int			iGP, nGP;
	int			GPID;
	PAIR_INT	pair_image;
	GP			ptGround;	

	// Get the size of m_vec_nimg_GP to find the largest GP ID. 
	// NOTE: The largest GP ID = nGP + m_base_gp_id
	nGP	= m_vec_nimg_GP.size();

	// Iterate through each index of GP (iGP+m_base_gp_id corresponding to GP ID)
	for (iGP = 0; iGP < nGP; iGP++)
	{
		// Get the corresponding GP ID.
		GPID = iGP + m_base_gp_id;

		// Find the pair of images with largest distance between PCs.
		if ( false == FindLargestDistancePair(GPID, vec_eo, start_img_id, last_img_id, nimages, pair_image) )
			continue;

		// Compute the 'rough' initial approximation of GP from the pair of tie points.
		ptGround = ComputeGPFromATiePoint(GPID, pair_image, vec_ip_f, vec_eo, idx_start_ip, idx_last_ip, io, rot);

		// Assign the initial approximation of GP into the vector.
		m_vec_initGPID.push_back(GPID);			// GP ID.
		m_vec_initGP.push_back(ptGround);		// GP approximated result.
	}

	// Assign the 'rough' init GP to the fine vector.
	copy(m_vec_initGP.begin(), m_vec_initGP.end(), std::back_inserter(m_vec_initGP2));

}

/********************************************************************
 * Description:	Compute the initial approximation of GP.
 *				The approximated GPs stored in 'm_vec_initGP2'.
 ********************************************************************/
void CGP::GPInitApproximation(	VEC_IP const &vec_ip_f, 
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
								bool b_fine_approx)
{
	int i, k, iGP, iIM;
	int	nGP, nIM;
	int	idxGP;
	int GP_ID, IM_ID;
	VEC_INT	vec_idx;
	GP	gp;

	// Compute the 'rough' initial approximation of GP.
	GPInitApproximation(vec_ip_f, 
						idx_start_ip,
						idx_last_ip,
						io, 
						vec_eo,
						start_img_id,
						last_img_id,
						nimages,
						rot,
						n_max_iter);

	if (b_fine_approx == false)
		return;

	// Compute the 'fine' initial approximation of GP.

	// Get the total number of GPs. 
	nGP	= m_vec_initGP.size();

	// Iterate through each GP.
	for (iGP = 0; iGP < nGP; iGP++)
	{	
		GP_ID = m_vec_initGPID[iGP];	// The GP ID.

		idxGP = GP_ID - m_base_gp_id;	// The index of GP ID.

		nIM = m_vec_nimg_GP[idxGP];		// The number of images for this GP.

		// Get the image ID in which this GP ID exists.
		i = -1;
		vec_idx.clear();
		while ( (i = m_CBits.GetNextBitSet(m_vec_imagebits[idxGP], i+1, nimages)) != -1 )
			vec_idx.push_back(i);

		// Assert to validate the accuracy of program.
		_ASSERT( nIM == vec_idx.size() );

		// Iterate through each round until reach the max. iterations.
		for (k = 0; k < n_max_iter; k++)
		{
			Matrix	Npp(3,3), Cp(3,1), kt;
			int		idxIP;

			idxIP = idx_start_ip;

			gp = m_vec_initGP2[iGP];		// The ground point.

			for (iIM = 0; iIM < nIM; iIM++)
			{
				// Get Image ID.
				IM_ID = vec_idx[iIM] + m_base_img_id;

				// Get the image point index.
				idxIP = GetImagePointIndex(vec_ip_f, IM_ID, GP_ID, idxIP);

				// Get the image point.
				IP ip_f = vec_ip_f[idxIP];

				// Get the EO of this image ID.
				EO eo = GetEO(vec_eo, IM_ID);

				// Compute the Ap and y matrix.
				Matrix Ap, y(2,1);
				Matrix rotM = rot.GetRotMatrix(IM_ID);	// Rotational matrix
				ComputeAy(ip_f, io, eo, gp, rotM, Ap, y);

				// Npp = Npp + Api' * Api
				Npp = Npp + Ap.Trans() * Ap;

				// Cp = Cp + Api' * yi
				Cp = Cp + Ap.Trans() * y;
			}

			// kt = inv(Npp) * Cp
			kt = Npp.Inv() * Cp;

			// GP_2(gi,:) = GP_2(gi,:) + kt'
			m_vec_initGP2[iGP].X = m_vec_initGP2[iGP].X + kt.Get(0, 0);
			m_vec_initGP2[iGP].Y = m_vec_initGP2[iGP].Y + kt.Get(1, 0);
			m_vec_initGP2[iGP].Z = m_vec_initGP2[iGP].Z + kt.Get(2, 0);

			// Break if its norm is less than the threshold.
			if ( DBL_LT( kt.Norm(0, false), threshold ) ) 
				break;
		}
	}
}

/********************************************************************
 * Description:	Fill the 'image bits' into the m_vec_imagebits vector.
 *				Count the number of images in which the GP exists.
 ********************************************************************/
void CGP::FillImageBits(IP const &ip)
{
	// Get the latest index in the vector.
	int iLatestGP = m_vec_imagebits.size() - 1;

	// Get the current index.
	int iCurrentGP = ip.iObj - m_base_gp_id;

	// Fill the blank GP between the latest and/include the current.
	for (int i = iLatestGP + 1; i <= iCurrentGP; i++)
	{
		m_vec_imagebits.push_back( DWORD_BITS() );
		m_vec_nimg_GP.push_back(0);
	}

	// Mark the 'image bits' in which this object appears.
	DWORD_BITS image_bits = m_vec_imagebits[iCurrentGP];
	m_CBits.SetBit(image_bits, (ip.iImg - m_base_img_id));
	m_vec_imagebits[iCurrentGP] = image_bits;

	// Increment the count of number of images the GP exists in.
	m_vec_nimg_GP[iCurrentGP] = m_vec_nimg_GP[iCurrentGP] + 1;
}


/********************************************************************
 * Description:	Find a pair of images that contributes largest Euclidean
 *				distance between PC of the pair.
 * Return:		true if able to find the largest distance pair.
 ********************************************************************/
bool CGP::FindLargestDistancePair(	int GPID, 
									VEC_EO const &vec_eo, 
									int start_img_id,
									int last_img_id,
									int const n_img_sequence, 
									/*out*/ PAIR_INT &pair_img)
{
	int			n_img_pt_exist;	// The number of images that this GP exists.
	double		dist, max_dist;	// The distance and maximum distance between an image pair.
	VEC_INT		vec_idx;		// The vector of image ID that this GP appears.
	int			i, j;			// The temporary image index.
	EO			EO_i, EO_j;		// The EO of image#i and image#j.
	int			iGP;			// Index of the GP ID in the vector.

	iGP = GPID - m_base_gp_id;

	// Consider only GP that exists in more than 1 image.
	if ((n_img_pt_exist = m_vec_nimg_GP[iGP]) < 2)
	{
		// Check if which image id (wrt m_base_img_id) contains this GP.
		i = m_CBits.GetNextBitSet(m_vec_imagebits[iGP], 0, n_img_sequence);

		// Decrease the number of shared TP/GP for this image by one.
		m_vec_n_validGP[i + m_base_img_id] = m_vec_n_validGP[i + m_base_img_id] - 1;

		return false;
	}
	
	// Simply return the pair if the point only appears on 2 images.
	if (2 == n_img_pt_exist)
	{
		pair_img.elm1 = m_CBits.GetNextBitSet(m_vec_imagebits[iGP], 0, n_img_sequence);
		pair_img.elm2 = m_CBits.GetNextBitSet(m_vec_imagebits[iGP], pair_img.elm1+1, n_img_sequence);

		if ( (-1 == pair_img.elm1) || (-1 == pair_img.elm2) )
			return false;
		else
		{
			// Update the actual image ID.
			pair_img.elm1 = pair_img.elm1 + m_base_img_id;
			pair_img.elm2 = pair_img.elm2 + m_base_img_id;

			return true;
		}
	}

	/*** GP appears on more than 2 image ***/

	// Initialization;
	max_dist = -1.0;
	pair_img.elm1 = -1;
	pair_img.elm2 = -1;

	// Obtain the image index in which this GP exists.
	i = -1;
	while ( (i = m_CBits.GetNextBitSet(m_vec_imagebits[iGP], i+1, n_img_sequence)) != -1 )
		vec_idx.push_back(i);

	// Iterate through each pair.
	for (i = 0; i < n_img_pt_exist; i++)
	{
		for (j = i+1; j < n_img_pt_exist; j++)
		{	
			// Obtain the actual image ID of index #i and #j.
			int IMG_i = vec_idx[i] + m_base_img_id;
			int IMG_j = vec_idx[j] + m_base_img_id;

			// Obtain the EO of Image#i and Image#j.
			EO_i = GetEO(vec_eo, IMG_i);
			EO_j = GetEO(vec_eo, IMG_j);

			// Calculate the euclidean distance between PCs of the image pair.
			dist = sqrt( (EO_i.dXc - EO_j.dXc)*(EO_i.dXc - EO_j.dXc) +
						 (EO_i.dYc - EO_j.dYc)*(EO_i.dYc - EO_j.dYc) +
						 (EO_i.dZc - EO_j.dZc)*(EO_i.dZc - EO_j.dZc) );

			// Update the pair of max distance if 'dist' is greater than 'max_distance'.
			if (DBL_GT(dist, max_dist))
			{
				max_dist = dist;
				pair_img.elm1 = IMG_i;
				pair_img.elm2 = IMG_j;
			}
		}
	}

	if ( (pair_img.elm1 == -1) || (pair_img.elm2 == -1) )
		return false;
	else
		return true;
}


/********************************************************************
 * Description:	Return the EO of each image.
 ********************************************************************/
EO CGP::GetEO(VEC_EO const &vec_eo, int imageID)
{
	int i;
	int vec_eoSize=0;
	vec_eoSize=vec_eo.size();

	for (i = 0; i < vec_eoSize; i++)
	{
		if (imageID == vec_eo[i].iImg)
			break;
	}

	if (i == vec_eo.size())
	{
		cout << "Error! The EO of the specified image does not exist." << endl;
		exit(0);
	}

	return vec_eo[i];
}

/********************************************************************
 * Description:	Get the image point vector of the specified GP ID
 *				and specified image ID in the form of 3x1 matrix:
 *				[x, y, -c]
 ********************************************************************/
Matrix CGP::GetImagePointVector(	int GPID,
									int imageID,
									VEC_IP const &vec_ip_f, 
									IO const &io)
{
	int		i;
	int vec_ip_fSize=0;
	vec_ip_fSize=vec_ip_f.size();
	Matrix	p(3,1);

	for (i = 0; i < vec_ip_fSize; i++)
	{
		if ( (vec_ip_f[i].iObj == GPID) && (vec_ip_f[i].iImg == imageID) )
		{
			// Form a vector of image point [x,y,-c]
			p.Set(0, 0, vec_ip_f[i].dX);
			p.Set(1, 0, vec_ip_f[i].dY);
			p.Set(2, 0, -io.dF);

			return p;
		}
	}
}

Matrix CGP::GetImagePointVector(	int GPID,
									int imageID,
									VEC_IP const &vec_ip_f, 
									int idx_start_ip,
									int idx_last_ip, 
									IO const &io)
{
	int		i;
	Matrix	p(3,1);

	for (i = idx_start_ip; i <= idx_last_ip; i++)
	{
		if ( (vec_ip_f[i].iObj == GPID) && (vec_ip_f[i].iImg == imageID) )
		{
			// Form a vector of image point [x,y,-c]
			p.Set(0, 0, vec_ip_f[i].dX);
			p.Set(1, 0, vec_ip_f[i].dY);
			p.Set(2, 0, -io.dF);

			return p;
		}
	}
}

/********************************************************************
 * Description:	Compute the GP from a pair of tie points.
 * Concept:		Both image points correspond to the same GP -> P.
 *				P = (L1)(R1)'p1 + C1
 *				P = (L2)(R2)'p2 + C2
 *				(L1)(R1)'p1 + C1 = (L2)(R2)'p2 + C2
 *				(L1)(R1)'p1 - (L2)(R2)'p2 = C2 - C1
 *				
 *				[(R1)'p1	-(R2)'p2] * [L1 L2]' = [C2 - C1]
 *				|--------A----------|	|--x--|	   |----y---|	
 *
 ********************************************************************/
GP CGP::ComputeGPFromATiePoint(	int GPID,
								PAIR_INT &pair_image,
								VEC_IP const &vec_ip_f, 
								VEC_EO const &vec_eo,
								int idx_start_ip,
								int idx_last_ip, 
								IO const &io, 
								CRotationalMatrix const &rot)
{
	Matrix	R1t, R2t;				// (R1)' and (R2)'
	Matrix	p1(3,1), p2(3,1);
	Matrix	A11, A12, A(3,2);
	Matrix	y(3,1);
	Matrix	C1(3,1), C2(3,1);
	Matrix	x, N, GP1, GP2;
	EO		EO1, EO2;
	GP		gp;
	
	// Construct A11 = (R1)'p1
	R1t	= (rot.GetRotMatrix(pair_image.elm1)).Trans();
	p1	= GetImagePointVector(GPID, pair_image.elm1, vec_ip_f, idx_start_ip, idx_last_ip, io);
	A11 = R1t * p1;

	// Construct A12 = (R2)'p2
	R2t = (rot.GetRotMatrix(pair_image.elm2)).Trans();
	p2	= GetImagePointVector(GPID, pair_image.elm2, vec_ip_f, idx_start_ip, idx_last_ip, io);
	A12	= R2t * p2;

	// Construct A matrix : A = [A11 -A12]
	A.Set(0, 0, A11.Get(0, 0));
	A.Set(1, 0, A11.Get(1, 0));
	A.Set(2, 0, A11.Get(2, 0));
	A.Set(0, 1, A12.Get(0, 0) * (-1.0));
	A.Set(1, 1, A12.Get(1, 0) * (-1.0));
	A.Set(2, 1, A12.Get(2, 0) * (-1.0));

	// Construct the y matrix
	EO1 = GetEO(vec_eo, pair_image.elm1);
	EO2 = GetEO(vec_eo, pair_image.elm2);
	y.Set(0, 0, EO2.dXc - EO1.dXc);
	y.Set(1, 0, EO2.dYc - EO1.dYc);
	y.Set(2, 0, EO2.dZc - EO1.dZc);

	// Obtain the scale: Ax = y --> x = inv(A'A)*A'*y
	N = A.Trans() * A;
	x = N.Inv() * A.Trans() * y;

	// Calculate the ground point.
	C1.Set(0, 0, EO1.dXc);
	C1.Set(1, 0, EO1.dYc);
	C1.Set(2, 0, EO1.dZc);
	C2.Set(0, 0, EO2.dXc);
	C2.Set(1, 0, EO2.dYc);
	C2.Set(2, 0, EO2.dZc);

	GP1 = A11 * x.Get(0,0) + C1;
	GP2 = A12 * x.Get(1,0) + C2;

	// Get the average of the two GPs.
	gp.X = (GP1.Get(0, 0) + GP2.Get(0, 0))/2;
	gp.Y = (GP1.Get(1, 0) + GP2.Get(1, 0))/2;
	gp.Z = (GP1.Get(2, 0) + GP2.Get(2, 0))/2;

	return gp;
}

/********************************************************************
 * Description:	Return the index of the specified GP ID in the passed-in vector.
 ********************************************************************/
int CGP::GetIndexGP(VEC_INT const &vec_id, int gpID)
{
	int i;
	int vec_idSize=0;
	vec_idSize=vec_id.size();

	for (i = 0; i < vec_idSize; i++)
	{
		if (gpID == vec_id[i])
			return i;
	}

	return -1;
}

/********************************************************************
 * Description:	Get the number of images that this GP exists.
 ********************************************************************/
int CGP::GetNumImagesGPExisted(int gpID)
{
	int idxGP = gpID - m_base_gp_id;

	if (m_vec_nimg_GP.size() != 0)
		return m_vec_nimg_GP[idxGP];
	else
		return -1;
}

/********************************************************************
 * Description:	Return the index of the specified GP ID in the distinct GP vector.
 ********************************************************************/
int CGP::GetIndexDistinctGP(int gpID)
{
	return GetIndexGP(m_vec_distinct_GPID, gpID);
}

/********************************************************************
 * Description:	Return the image point which has the specified Img ID and GP ID.
 ********************************************************************/
IP CGP::GetImagePoint(VEC_IP const &vec_ip, int ImageID, int GPID, int iStartHint)
{
	int i, nSize;
	
	nSize = vec_ip.size();

	for (i = iStartHint; i < nSize; i++)
	{
		if ( (vec_ip[i].iImg != ImageID) || (vec_ip[i].iObj != GPID) )
			continue;
		
		return vec_ip[i];
	}

	// Falling through this point is not found.
	return IP();
}

/********************************************************************
 * Description:	Return the index of the specified GP ID from the list of GP Init.
 ********************************************************************/
int CGP::GetIndexInitGP(int gpID)
{
	return GetIndexGP(m_vec_initGPID, gpID);
}

/********************************************************************
 * Description:	Return the image point which has the specified Img ID and GP ID.
 ********************************************************************/
int CGP::GetImagePointIndex(VEC_IP const &vec_ip, int ImageID, int GPID, int iStartHint)
{
	int i, nSize;
	
	nSize = vec_ip.size();

	for (i = iStartHint; i < nSize; i++)
	{
		if ( (vec_ip[i].iImg != ImageID) || (vec_ip[i].iObj != GPID) )
			continue;
		
		return i;
	}

	// Not found.
	return -1;
}

/********************************************************************
 * Description:	Compute the Ap and y matrices.
 ********************************************************************/
void CGP::ComputeAy(/*in*/	IP const &ip_f,
					/*in*/	IO const &io, 
					/*in*/	EO const &eo, 
					/*in*/	GP const &gp, 
					/*in*/	Matrix const &R, // Rotational Matrix	
					/*out*/	Matrix &Ap,
					/*out*/	Matrix &y)
{
	// Variable declaration
	Matrix	GC(3,1), ND(3,1), F0(2,1), dND(3,3);
	Matrix	M;		// Temporary
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

	dND = R;

	// Ap = IO(3) / ND(3)^2 * [-ND(3) 0 ND(1); 0 -ND(3) ND(2)] * dND
	M.SetDim(2, 3);		// For [-ND(3) 0 ND(1); 0 -ND(3) ND(2)]

	M.Set(0, 0, -ND.Get(2, 0));
	M.Set(0, 1, 0.0);
	M.Set(0, 2, ND.Get(0, 0));

	M.Set(1, 0, 0.0);
	M.Set(1, 1, -ND.Get(2, 0));
	M.Set(1, 2, ND.Get(1, 0));

	dTemp = io.dF / (ND.Get(2, 0) * ND.Get(2, 0));
	Ap = M * dND * dTemp;

	// y = IP - F0
	y.Set(0, 0, ip_f.dX - F0.Get(0, 0));
	y.Set(1, 0, ip_f.dY - F0.Get(1, 0));

}
