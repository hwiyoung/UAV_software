/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "AT_SimAT.h"
#include "AT_Global.h"
#include <iostream>
#include <fstream>		// To utilize ifstream for dealing with text file input.
#include <strstream>
#include "stdio.h"
#include <conio.h>

/********************************************************************
 * Description:	Construction
 ********************************************************************/
AT_CSimAT::AT_CSimAT(const string &_CFG_File)
{
	m_CFG_File = _CFG_File;

}

/********************************************************************
 * Description:	Main function to perform Simultaneous AT
 ********************************************************************/
void AT_CSimAT::RunSimAT()
{
	cout << "********************* " << endl;
	cout << "Simultaneous AT Start " << endl;
	cout << "********************* " << endl;
	// Initialize all neccessary parameters.
	Initialization();

	// Perform Simultaneous AT estimation.
	SimATEstimation();
}

/********************************************************************
 * Description:	Call necessary functions for initialization prior to
 *				perform the simultaneous AT.
 ********************************************************************/
void AT_CSimAT::Initialization()
{

	// Read the configuration file.
	ReadConfig();

	// Read the IO file.
	ReadIOFile();

	// Read the EO file.
	ReadEOFile();

	// Read the IP/Tiepoints file.
	ReadIPFile();

	// Initialize member variables prior to processing.
 	InitializeMembers();

}

/********************************************************************
 * Description:	Read the configuration file
 ********************************************************************/
void AT_CSimAT::ReadConfig()
{
	// Open the configuration file to read.
	ifstream fileInput(m_CFG_File.c_str());

	if (!fileInput)
	{
		cout <<"Input 폴더의 Config 파일을 확인할 수 없습니다" << endl;
		cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
		_getch();			
		exit(0);
	}
	
	const int BUFFER_SIZE = 200;
	char buffer[BUFFER_SIZE+1] = {0};

	// Read the IO file name.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	istrstream(buffer) >> m_IO_File;

	// Read the EO file name.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	istrstream(buffer) >> m_EO_File;

	// Read the TP file name.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	istrstream(buffer) >> m_IP_File;

	// Read the max. number of iterations for performing AT.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf_s(buffer, "%d", &m_max_iter);

	// Read the initial number of images for starting seq. AT sim.
	//GetLine(fileInput, buffer, BUFFER_SIZE);
	//sscanf(buffer, "%d", &m_ninit_sim_img);
	
	
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
void AT_CSimAT::ReadIOFile()
{
	// Open the IO file to read.
	ifstream fIO(m_IO_File.c_str());

	if (!fIO)
	{
		cout <<"Input 폴더의 IO 파일을 확인할 수 없습니다" << endl;
		cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
		_getch();			
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
 ********************************************************************/
void AT_CSimAT::ReadEOFile()
{
	AT_EO	eo;

	// If the file exist, read the EO data from the file.
	FILE * fEO;
	fEO = _fsopen(m_EO_File.c_str(), "r", _SH_DENYWR);
		
	if(!fEO)
	{
		cout <<"Results 폴더의 EO 파일을 확인할 수 없습니다" << endl;
		cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
		_getch();
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
void AT_CSimAT::ReadIPFile()
{
	AT_IP	ip;
	int nImg, numIP, imgID, nIP;
	
	// If the file exist, read the tiepoints from the file.
	FILE * fIP;
	fIP = _fsopen(m_IP_File.c_str(), "r", _SH_DENYWR);

	if(!fIP)
	{
		cout <<"Results\\Matching_Results 폴더의 TiePoint_fcs_AT 파일을 확인할 수 없습니다" << endl;
		cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
		_getch();			
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
 * Description:	Initialize EO.
 ********************************************************************/
void AT_CSimAT::InitializeEO()
{
	// Compute the 3D rotational matrix.
	m_RotMatrix.Compute3DRotMatrix(m_EO_m, m_nImages);
	
	// Compute the derivative 3D rotational matrix.
	m_RotMatrix.ComputeDerivative3DRotMatrix(m_EO_m, m_nImages);
}

/********************************************************************
 * Description:	Initialize ground points.
 ********************************************************************/
void AT_CSimAT::InitializeGP()
{
	// Make the list of distinct ground point ID
	m_GP.MakeGPList(m_IP_i, m_nImages);

	// Debug out the list of distinct GP.
	m_GP.DebugDistinctGP();

	// Initial approximation for GPs
	m_GP.GPInitApproximation(m_IP_f, m_IO, m_EO_m, m_nImages, m_RotMatrix, m_max_iter, m_dDelta, true);

	// Debug out the inital approximated GPs.
	m_GP.DebugInitGP(false);
}

/********************************************************************
 * Description:	Initialize member variables prior to processing.
 ********************************************************************/
void AT_CSimAT::InitializeMembers()
{
	// Initialize the photo coordinate object.
	//m_PhotoCoord.Initialization(m_IO, m_IS, m_PS);
	m_IP_f=m_IP_i;
}


/********************************************************************
 * Description:	Wrapper function to get each line from the input file.
 ********************************************************************/
void AT_CSimAT::GetLine(istream &in, char *buf, const int bufSize) const
{
	while(in.good())
	{
		in.getline(buf, bufSize);
		if(buf[0]!='%') break;
	}
}

/********************************************************************
 * Description:	Perform Simultaneous AT estimation.
 ********************************************************************/
void AT_CSimAT::SimATEstimation()
{
	// Set the parameters to the Sim AT object's member variables.
	m_ninit_sim_img=m_nImages-1;
	//m_ninit_sim_img=1;
	m_SeqSimAT.Set(	m_nImages, m_max_iter, m_ninit_sim_img,	m_dDelta, m_std_IP, m_std_GPS, m_std_INS);

	// Perform the sequential AT simultaneous estimation.
	m_SeqSimAT.SeqSimEstimation(m_IP_f,m_vec_n_IP,m_IO,m_EO_m);

	// The results from Sim AT are its member variable.

	// Debug out the results into a file named "Summary.txt".
	m_SeqSimAT.DebugSimATResults();
}
