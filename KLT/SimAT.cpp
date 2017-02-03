/********************************************************
 * Project:			Automated AT
 * Last updated:	17 Jaunary 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "SimAT.h"
#include "Global.h"
#include <iostream>
#include <fstream>		// To utilize ifstream for dealing with text file input.
#include <strstream>
#include "stdio.h"

/********************************************************************
 * Description:	Construction
 ********************************************************************/
CSimAT::CSimAT(const string &_CFG_File)
{
	// Store the configuration file into the member variable.
	m_CFG_File = _CFG_File;

	// Initialize the epoch ID.
	m_EPOCH_ID = 0;

	// Initialize the number of images in the sequence.
	m_nImages = 0;

	// Initialize the sum of image points.
	m_nIP_raw = 0;

	// Initialize the number of points for the first image ID.
	// Image ID: 0 (not existence) always has 0 image point.
	m_vec_n_IP.push_back(0);	

}

/********************************************************************
 * Description:	Main function to perform Simultaneous AT
 ********************************************************************/
void CSimAT::RunSimAT()
{
	// Initialize all neccessary parameters.
	Initialization();

	// Perform Simultaneous AT estimation.
	SimATEstimation();
}

/********************************************************************
 * Description:	Call necessary functions for initialization prior to
 *				perform the simultaneous AT.
 ********************************************************************/
void CSimAT::Initialization()
{
	// Read the configuration file.
	ReadConfig();

	// Read the IO file.
	ReadIOFile();

	// Initialize the output tie point file.
	if (true == m_bOutputTP)
		InitializeTPFile();
}

/********************************************************************
 * Description:	Read the configuration file
 ********************************************************************/
void CSimAT::ReadConfig()
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

	// Read the flag whether to run this process after TP extraction.
	//GetLine(fileInput, buffer, BUFFER_SIZE);
	//sscanf(buffer, "%d", &m_bRun);
	m_bRun=0;

	// Read the flag whether to out put tie points into a file.
	//GetLine(fileInput, buffer, BUFFER_SIZE);
	//sscanf(buffer, "%d", &m_bOutputTP);
	m_bOutputTP=1;

	// Read the IO file name.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	istrstream(buffer) >> m_IO_File;

	// Read the rotational matrix order.
	//GetLine(fileInput, buffer, BUFFER_SIZE);
	//sscanf(buffer, "%d", &ROT_ORDER);
	ROT_ORDER=2;

	// Read the rotational matrix type.
	//GetLine(fileInput, buffer, BUFFER_SIZE);
	//sscanf(buffer, "%d", &ROT_TYPE);
	ROT_TYPE=2;

	// Read the max. number of iterations for performing AT.
	//GetLine(fileInput, buffer, BUFFER_SIZE);
	//sscanf(buffer, "%d", &m_max_iter);
	m_max_iter=15;

	// Read the initial number of images for starting seq. AT sim.
	//GetLine(fileInput, buffer, BUFFER_SIZE);
	//sscanf(buffer, "%d", &m_ninit_sim_img);
	m_ninit_sim_img=2;

	// Read the max. number of images for seq. AT computation.
	//GetLine(fileInput, buffer, BUFFER_SIZE);
	//sscanf(buffer, "%d", &m_nmax_sim_img);
	m_nmax_sim_img=50;

	GetLine(fileInput, buffer, BUFFER_SIZE);
	GetLine(fileInput, buffer, BUFFER_SIZE);

	// Read the 'Delta' stop criteria.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf_s(buffer, "%lg", &m_dDelta);

	// Read the standard deviation of the image point errors (pixels).
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf_s(buffer, "%lf", &m_std_IP);

	// Read the standard deviation of the GPS errors (cm).
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf_s(buffer, "%lf", &m_std_GPS);

	// Read the standard deviation of the INS errors (degrees).
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf_s(buffer, "%lf", &m_std_INS);

	// Close the file after finish reading.
	fileInput.close();
}

/********************************************************************
 * Description:	Read the IO file
 ********************************************************************/
void CSimAT::ReadIOFile()
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
	sscanf_s(buffer, "%lf,%lf,%lf", &m_IO.dPPX, &m_IO.dPPY, &m_IO.dF);

	// Read the radial distortion.
	GetLine(fIO, buffer, BUFFER_SIZE);
	sscanf_s(buffer, "%lg,%lg,%lg", &m_RDC.dK1, &m_RDC.dK2, &m_RDC.dK3);
	
	// Read the resolution / image dimension
	GetLine(fIO, buffer, BUFFER_SIZE);
	sscanf_s(buffer, "%d,%d", &m_IS.iWidth, &m_IS.iHeight);

	// Read the pixel size (mm)
	GetLine(fIO, buffer, BUFFER_SIZE);
	sscanf_s(buffer, "%lf,%lf", &m_PS.dX, &m_PS.dY);

	// Close the file after complete reading.
	fIO.close();
}

/********************************************************************
 * Description:	Read the EO file
 *
 * NOTE:		This function is substituted by the GetNewImage()
 *				function which will store the EO parameters of 
 *				each image into the m_EO_m buffer.
 ********************************************************************/
void CSimAT::ReadEOFile()
{
	EO		eo;

	// If the file exist, read the EO data from the file.
	FILE * fEO;

	if ((fEO = _fsopen(m_EO_File.c_str(), "r", _SH_DENYWR)) != NULL)
	{
		cout << "Error! Failed to open the EO file." << endl;
		exit(0);
	}

	const int BUFFER_SIZE = 200;
	char buffer[BUFFER_SIZE+1] = {0};

	while(fgets(buffer, BUFFER_SIZE, fEO) != NULL)
	{
		// Skip the comment
		if (buffer[0] == '%')
			continue;

		sscanf_s(buffer, "%d\t%lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", &eo.iImg, &eo.dXc, &eo.dYc, &eo.dZc, &eo.dOmega, &eo.dPhi, &eo.dKappa);

		m_EO_m.push_back(eo);
	}

	fclose(fEO);
}

/********************************************************************
 * Description:	Read the IP/Tiepoints file
 ********************************************************************/
void CSimAT::ReadIPFile()
{
	IP	ip;
	int nImg, numIP, imgID, nIP;

	// If the file exist, read the tiepoints from the file.
	FILE * fIP;

	if((fIP = _fsopen(m_IP_File.c_str(), "r", _SH_DENYWR)) != NULL)
	{
		cout << "Error! Failed to open the IP file." << endl;
		exit(0);
	}

	// Read the number of images in the sequence.
	fscanf_s(fIP, "No. Images : %d\n", &m_nImages);

	if (m_nImages <= 0)
	{
		cout << "Error! Number of image to process is invalid." << endl;
		exit(0);
	}

	// Image ID: 0 (not existence) always has 0 image point.
	m_vec_n_IP.push_back(0);	

	// Initialize the sum of image points.
	m_nIP_raw = 0;

	// Iterate through each image in the sequence.
	for (nImg = 1; nImg <= m_nImages; nImg++)
	{
		// Read the Image ID and number of image points in that image.
		fscanf_s(fIP, "Image ID : %d\n", &imgID);
		fscanf_s(fIP, "No. IP : %d\n", &numIP);

		// Add the number of image points into the vector.
		m_vec_n_IP.push_back(numIP);

		// Increment the sum of image points.
		m_nIP_raw += numIP;

		// Store each point into the vector.
		for (nIP = 0; nIP < numIP; nIP++)
		{
			// Obain the 'ImageID', 'ObjectID', 'X-coord', 'Y-coord' in order.
			fscanf_s(fIP, "%d\t%d\t%lf\t%lf\n", &ip.iImg, &ip.iObj, &ip.dX, &ip.dY);

			// Put into the list of points. 
			m_IP_i.push_back(ip);
		}
	}
	
	// Close the file after finish reading.
	fclose(fIP);
}

/********************************************************************
 * Description:	Read the IP/Tiepoints file generated from Simulator
 ********************************************************************/
void CSimAT::ReadIPFile2()
{
	IP	ip;
	int nImg, numIP, imgID, nIP;

	// If the file exist, read the tiepoints from the file.
	FILE * fIP;
	
	if((fIP = _fsopen(m_IP_File.c_str(), "r", _SH_DENYWR)) != NULL)
	{
		cout << "Error! Failed to open the IP file." << endl;
		exit(0);
	}

	// Read the number of images in the sequence.
	fscanf_s(fIP, "No. Images : %d\n", &m_nImages);

	if (m_nImages <= 0)
	{
		cout << "Error! Number of image to process is invalid." << endl;
		exit(0);
	}

	// Initialize the number of image points in the m_vec_n_IP vector.
	for (nImg = 0; nImg <= m_nImages; nImg++)
		m_vec_n_IP.push_back(0);	

	// Initialize the sum of image points.
	m_nIP_raw = 0;

	// Iterate through each image in the sequence starting from Image#1.
	for (nImg = 1; nImg <= m_nImages; nImg++)
	{
		// Read the Image ID and number of image points in that image.
		fscanf_s(fIP, "Image ID : %d\n", &imgID);
		fscanf_s(fIP, "No. IP : %d\n", &numIP);

		// Increment the sum of image points.
		m_nIP_raw += numIP;

		// Store each point into the vector.
		for (nIP = 0; nIP < numIP; nIP++)
		{
			// Obain the 'ImageID', 'ObjectID', 'X-coord', 'Y-coord' in order.
			fscanf_s(fIP, "%d\t%d\t%lf\t%lf\n", &ip.iImg, &ip.iObj, &ip.dX, &ip.dY);

			// Put into the list of points. 
			m_IP_i.push_back(ip);

			// Increment the number of image points in the vector.
			m_vec_n_IP[ip.iImg] = m_vec_n_IP[ip.iImg] + 1;
		}
	}
	
	// Close the file after finish reading.
	fclose(fIP);
}

/********************************************************************
 * Description:	Initialize EO.
 ********************************************************************/
void CSimAT::InitializeEO()
{
	// Compute the 3D rotational matrix.
	m_RotMatrix.Compute3DRotMatrix(m_EO_m, m_nImages, ROT_ORDER, ROT_TYPE);
	
	// Compute the derivative 3D rotational matrix.
	m_RotMatrix.ComputeDerivative3DRotMatrix(m_EO_m, m_nImages, ROT_ORDER, ROT_TYPE);
}

/********************************************************************
 * Description:	Initialize ground points.
 ********************************************************************/
void CSimAT::InitializeGP()
{
/*
	// Make the list of distinct ground point ID
	m_GP.MakeGPList(m_IP_i, m_nImages);

	// Debug out the list of distinct GP.
	m_GP.DebugDistinctGP();

	// Initial approximation for GPs
	m_GP.GPInitApproximation(m_IP_f, m_IO, m_EO_m, m_nImages, m_RotMatrix, m_max_iter, m_dDelta, true);

	// Debug out the inital approximated GPs.
	m_GP.DebugInitGP(false);
*/
}

/********************************************************************
 * Description:	Initialize member variables prior to processing.
 ********************************************************************/
void CSimAT::InitializeMembers()
{
	// Convert the units of members to appropriate scale.
	// UnitConversion();

	// Initialize the photo coordinate object.
	// m_PhotoCoord.Initialization(m_IO, m_IS, m_PS);

	// Convert the units of std. deviation of the error.
    // StdUnitConversion();

	// Convert the image points from 'Pixel Coord' to 'Photo Coord'
	// m_IP_f = m_PhotoCoord.ConvertPixel2PhotoCoords(m_IP_i);

	// Initialize EO parameters.
	//InitializeEO();

}


/********************************************************************
 * Description:	Convert the units of input parameters to appropriate scale.
 ********************************************************************/
void CSimAT::UnitConversion()
{
	// Principle point - convert from mm to m
	m_IO.dPPX	=	m_IO.dPPX/1000.0;
	m_IO.dPPY	=	m_IO.dPPY/1000.0;
	m_IO.dF		=	m_IO.dF/1000.0;

	// Pixel size - convert from mm to m
	m_PS.dX		=	m_PS.dX/1000.0;
	m_PS.dY		=	m_PS.dY/1000.0;
}


/********************************************************************
 * Description:	Convert the units of std. deviation of error to appropriate scale.
 ********************************************************************/
void CSimAT::StdUnitConversion()
{
	// Std. of IP error - convert the unit from pixels to m
	PS ps = m_PhotoCoord.ConvertNumPixel2Metre(m_std_IP);
	m_std_IP = (ps.dX + ps.dY)/2;		// Average value

	// Std. of GPS error - convert the unit from cm to m
	m_std_GPS = m_std_GPS/100.0;

	// Std. of INS error - convert the unit from degree to radian
	m_std_INS = deg2rad(m_std_INS);
}

/********************************************************************
 * Description:	Wrapper function to get each line from the input file.
 ********************************************************************/
void CSimAT::GetLine(istream &in, char *buf, const int bufSize) const
{
	while(in.good())
	{
		in.getline(buf, bufSize);
		if(buf[0]!='%') break;
	}
}

/********************************************************************
 * Description:	This function will be called whenever a new image
 *				is captured. 
 ********************************************************************/
void CSimAT::GetNewImage(string image_name, EO eo)
{
	// Store the EO for future use while discard the unused image file name.
	m_EO_m.push_back(eo);

	// Increment the number of images in the sequence.
	m_nImages++;

	// Get the tie points from the KLT object.
	if (m_EO_m.size() >= 2)
	{
		GetTiePoints();

		// At this point, the TP is in the form of photo-coordinates.
	}

	// Perform the Simultaneous AT if the number of images is greater
	// than the minimum number of the required images to perform.
	if (true == m_bRun)
	{
		int m_EO_mSize=0;
		m_EO_mSize=m_EO_m.size();
		if (m_EO_mSize >= m_ninit_sim_img)
			SimATEstimation();
	}
}

/********************************************************************
 * Description:	Perform Simultaneous AT estimation.
 ********************************************************************/
void CSimAT::SimATEstimation()
{
	// Set the parameters to the Sim AT object's member variables.
	m_SeqSimAT.Set(	m_nImages, 
					m_max_iter, 
					m_ninit_sim_img, 
					m_nmax_sim_img,
					m_dDelta,
					m_std_IP,
					m_std_GPS,
					m_std_INS,
					ROT_ORDER,
					ROT_TYPE);

	// Get the latest image ID.
	int latest_image_id = m_EO_m[m_nImages-1].iImg;

	// Perform the sequential AT simultaneous estimation.
	m_SeqSimAT.SeqSimEstimation(m_vec_n_IP,
								m_IP_f, 
								m_vec_n_IP, 
								m_IO, 
								m_EO_m,
								latest_image_id);

	// The results from Sim AT are its member variable.

}

/********************************************************************
 * Description:	Print out the result of the Simultaneous AT estimation.
 ********************************************************************/
void CSimAT::PrintOutResults()
{

	IP	ip;
	int nImg, numIP, imgID, nIP, all_Imgs;

	FILE * fIP = _fsopen("Results\\Matching_Results\\TiePoint_ics.txt", "r", _SH_DENYWR);
	fscanf_s(fIP, "No. Images : %d\n", &all_Imgs);

	Matrix	image_ip;
	image_ip.SetDim(all_Imgs, 1);

	r_vec_n_IP.push_back(0);
	m_nIP_raw = 0;

	for (nImg = 1; nImg < all_Imgs; nImg++)
	{
		fscanf_s(fIP, "Image ID : %d\n", &imgID);
		fscanf_s(fIP, "No. IP : %d\n", &numIP);
	
		r_vec_n_IP.push_back(numIP);


		m_nIP_raw += numIP;


		for (nIP = 0; nIP < numIP; nIP++)
		{			
			fscanf_s(fIP, "%d\t%d\t%lf\t%lf\n", &ip.iImg, &ip.iObj, &ip.dX, &ip.dY);
			
			if (ip.iImg == nImg)
			{
				image_ip.Set(nImg-1, image_ip.Get(nImg-1) + 1);
			}	
			else
			{
				image_ip.Set(nImg, image_ip.Get(nImg) + 1);
			}
			r_IP_i.push_back(ip);
		}
	}
	fclose(fIP);
	
	r_IP_f = m_PhotoCoord.ConvertPixel2PhotoCoords(r_IP_i);

	FILE * iFpoch = _fsopen("Results\\Matching_Results\\TiePoint_ics_AT.txt", "w", _SH_DENYWR);
	fclose(iFpoch);
	FILE * fFpoch = _fsopen("Results\\Matching_Results\\TiePoint_fcs_AT.txt", "w", _SH_DENYWR);
	fclose(fFpoch);
	FILE * fmatlab_fcs = _fsopen("Results\\Matching_Results\\TiePoint_MATLAB_fcs.txt", "w", _SH_DENYWR);
	fclose(fmatlab_fcs);
	FILE * fmatlab_ics = _fsopen("Results\\Matching_Results\\TiePoint_MATLAB_ics.txt", "w", _SH_DENYWR);
	fclose(fmatlab_ics);

	iFpoch = _fsopen("Results\\Matching_Results\\TiePoint_ics_AT.txt", "a+", _SH_DENYWR);
	fprintf(iFpoch, "No. Images : %d\n", all_Imgs); 

	fFpoch = _fsopen("Results\\Matching_Results\\TiePoint_fcs_AT.txt", "a+", _SH_DENYWR);
	fprintf(fFpoch, "No. Images : %d\n", all_Imgs); 

	fmatlab_fcs = _fsopen("Results\\Matching_Results\\TiePoint_MATLAB_fcs.txt", "a+", _SH_DENYWR);
	fmatlab_ics = _fsopen("Results\\Matching_Results\\TiePoint_MATLAB_ics.txt", "a+", _SH_DENYWR);

	int k=0;
	int noip=0;
	for (int i = 0; i < all_Imgs; i++)
	{
		noip=(int)image_ip[i];
		fprintf(iFpoch, "Image ID : %d\n", i+1);
		fprintf(iFpoch, "No. IP : %d\n", noip);	
		fprintf(fFpoch, "Image ID : %d\n", i+1);
		fprintf(fFpoch, "No. IP : %d\n", noip);	

		for(int j=k; j < k+image_ip[i]; j++)
		{
			fprintf(iFpoch, "%d\t%d\t%lf\t%lf\n", r_IP_i[j].iImg, r_IP_i[j].iObj, r_IP_i[j].dX, r_IP_i[j].dY);
			fprintf(fFpoch, "%d\t%d\t%lf\t%lf\n", r_IP_f[j].iImg, r_IP_f[j].iObj, r_IP_f[j].dX, r_IP_f[j].dY);
			fprintf(fmatlab_fcs, "%d\t%d\t%lf\t%lf\n", r_IP_f[j].iImg, r_IP_f[j].iObj, r_IP_f[j].dX, r_IP_f[j].dY);
			fprintf(fmatlab_ics, "%d\t%d\t%lf\t%lf\n", r_IP_i[j].iImg, r_IP_i[j].iObj, r_IP_i[j].dX, r_IP_i[j].dY);
		}
		k=k+(int)image_ip[i];
	}

	fclose(iFpoch);
	fclose(fFpoch);
	fclose(fmatlab_fcs);
	fclose(fmatlab_ics);
	
}

/********************************************************************
 * Description:	Store the pointer to the singleton CKLTTracker object.
 ********************************************************************/
void CSimAT::RegisterKLT(CKLTTracker * pKLT)
{
	m_pKLT = pKLT;
}


void CSimAT::GetTiePoints()
{
	int npoints, npoints1, npoints2;	// The number of tie points.
	int	IMG_ID;							// Image ID
	int	last_element_idx;				// Last element index.
	int last_stored_GP_ID;				// Latest stored GP ID
	int	idxIMG;							// Index of Image ID
	
	// Start the first epoch by requesting the tie points for the first two images.
	if (0 == m_EPOCH_ID)
	{
		// Get the tie points from the first image in the EO vector : m_EO_m[0].iImg
		npoints1 = m_pKLT->GetTiePoints(m_EO_m[0].iImg, &m_IP_f);

		// Store the number of tie points into the vector of tie-point number.
		m_vec_n_IP.push_back(npoints1);

		// Get the tie points from the second image in the EO vector : m_EO_m[1].iImg
		npoints2 = m_pKLT->GetTiePoints(m_EO_m[1].iImg, &m_IP_f);

		// Store the number of tie points into the vector of tie-point number.
		m_vec_n_IP.push_back(npoints2);

		// Calculate the total number of tie points for EPOCH ID#0
		npoints = npoints1 + npoints2;

	}
	else
	{
		// Get the first set of EPOCH for the 'first' image.

			// Get the image ID by computing from the epoch ID.
			IMG_ID = m_EO_m[m_EPOCH_ID].iImg;

			// Get the latest stored GP ID.
			last_element_idx = m_IP_f.size() - 1;
			last_stored_GP_ID = m_IP_f[last_element_idx].iObj;

			// Get the tie points.
			npoints1 = m_pKLT->GetTiePoints(IMG_ID, &m_IP_f, last_stored_GP_ID);

			// Update the number of tie points into the vector of tie-point number.
			idxIMG = m_vec_n_IP.size() - 1;
			m_vec_n_IP[idxIMG] = m_vec_n_IP[idxIMG] + npoints1;

		// Get the second set of EPOCH for the 'second' image.

			// Get the image ID by computing from the epoch ID.
			IMG_ID = m_EO_m[m_EPOCH_ID + 1].iImg;

			// Get the tie points.
			npoints2 = m_pKLT->GetTiePoints(IMG_ID, &m_IP_f);

			// Store the number of tie points into the vector of tie-point number.
			m_vec_n_IP.push_back(npoints2);
	}

	// Calculate the total number of tie points for EPOCH.
	npoints = npoints1 + npoints2;

	// Record the number of tie points for each EPOCH
	m_vec_n_EPOCH.push_back(npoints);	

	// Increment the sum of image points.
	m_nIP_raw += npoints;

	// Increment the epoch ID.
	m_EPOCH_ID++;	

	// To print out the immediate TP.
	if (true == m_bOutputTP)
	{
		FILE * fEpoch = _fsopen("Results\\Matching_Results\\TiePoint_fcs.txt", "a+", _SH_DENYWR);
		fprintf(fEpoch, "Image ID : %d\n", m_EPOCH_ID); 
		fprintf(fEpoch, "No. IP : %d\n", npoints);
		int idxIP = m_nIP_raw - npoints;
		for (int i = idxIP; i < m_nIP_raw; i++)
		{
			fprintf(fEpoch, "%d\t%d\t%lf\t%lf\n", 
				m_IP_f[i].iImg, m_IP_f[i].iObj, m_IP_f[i].dX, m_IP_f[i].dY);
		}
		fclose(fEpoch);
	}

	m_PhotoCoord.Initialization(m_IO, m_IS, m_PS);
	m_IP_i = m_PhotoCoord.ConvertPhoto2PixelCoords(m_IP_f);
	
	if (true == m_bOutputTP)
	{
		FILE * fFpoch = _fsopen("Results\\Matching_Results\\TiePoint_ics.txt", "a+", _SH_DENYWR);
		fprintf(fFpoch, "Image ID : %d\n", m_EPOCH_ID); 
		fprintf(fFpoch, "No. IP : %d\n", npoints);
		int idxIP = m_nIP_raw - npoints;
		for (int i = idxIP; i < m_nIP_raw; i++)
		{
			fprintf(fFpoch, "%d\t%d\t%lf\t%lf\n", 
				m_IP_i[i].iImg, m_IP_i[i].iObj, m_IP_i[i].dX, m_IP_i[i].dY);
		}
		fclose(fFpoch);
	}

}

/********************************************************************
 * Description:	Initialize the output tie point file.
 ********************************************************************/
void CSimAT::InitializeTPFile()
{
	//FILE * fEpoch = fopen("Matching_Results\\IP_s.txt", "w");
	//fclose(fEpoch);
}
