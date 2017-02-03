/********************************************************
 * Project:			Automated AT
 * Last updated:	25 Jaunary 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "AerialSimulator.h"
#include "Global.h"
#include <iostream>
#include <fstream>		// To utilize ifstream for dealing with text file input.
#include <strstream>
#include "stdio.h"
#include <direct.h>       //mkdir
#include <conio.h>

/********************************************************************
 * Description:	Construction
 ********************************************************************/
CAerialSimulator::CAerialSimulator(const string &_CFG_File)
{
	m_CFG_File = _CFG_File; 

	RUNNING_INDEX = 0;			// Initialize the image index.
}

/********************************************************************
 * Description:	Set member variables
 ********************************************************************/
void CAerialSimulator::SetMembers()
{

}

/********************************************************************
 * Description:	Main function to start the Aerial section.  
 ********************************************************************/
void CAerialSimulator::StartAerialSimulator()
{
	// Initialize all neccessary parameters.
	Initialization();
}

/********************************************************************
 * Description:	Call necessary functions for initialization.
 ********************************************************************/
void CAerialSimulator::Initialization()
{
	_mkdir("Results");

	ReadConfig();

	ReadEOFile();

}

/********************************************************************
 * Description:	Return the period of time to capture a new image (second)
 ********************************************************************/
double CAerialSimulator::GetExposurePeriod()
{
	return CAPTURE_PERIOD;
}

/********************************************************************
 * Description:	Store the pointer to the singleton CKLTTracker object.
 ********************************************************************/
void CAerialSimulator::RegisterKLTClient(CKLTTracker * pKLT)
{
	m_pKLT = pKLT;
}

/********************************************************************
 * Description:	Store the pointer to the singleton CSimAT object.
 ********************************************************************/
void CAerialSimulator::RegisterSimAT(CSimAT * pSimAT)
{
	m_pSimAT = pSimAT;
}

/********************************************************************
 * Description:	This function is to simulated that the system has 
 *				capture a new image. It then manages to notify its
 *				registered clients that a new image is captured.
 ********************************************************************/
void CAerialSimulator::RunImageExposure()
{
	int nIMG;		// Total number of images.
	int idxIMG;		// Index of image and EO.

	// Get the total number of images.
	nIMG = m_deq_IMG_File.size();

	// Iterate through the deque of the image file and eo.	
	for (idxIMG = 0; idxIMG < nIMG; idxIMG++)
	{
		// Assume that a new image is captured.
		if (1 == ZA_FILE_EXIST)
			OnExposureEventFired(m_deq_IMG_File[idxIMG], m_deq_EO[idxIMG], m_deq_Za[idxIMG]);		
		else
			OnExposureEventFired(m_deq_IMG_File[idxIMG], m_deq_EO[idxIMG]);		
	}

}

/********************************************************************
 * Description:	Handle the process once the exposure event fired.
 ********************************************************************/
void CAerialSimulator::OnExposureEventFired(string image_name, EO eo)
{
	// Notify the KLT tracker that a new image is captured.
	m_pKLT->GetNewImage(image_name, eo);

	// Notify the Simultaneous AT that a new image is captured.
	m_pSimAT->GetNewImage(image_name, eo);
}

void CAerialSimulator::OnExposureEventFired(string image_name, EO eo, double za)
{
	// Notify the KLT tracker that a new image is captured.
	m_pKLT->GetNewImage(image_name, eo, za);

	// Notify the Simultaneous AT that a new image is captured.
	m_pSimAT->GetNewImage(image_name, eo);
}

/********************************************************************
 * Description:	Wrapper function to get each line from the input file.
 ********************************************************************/
void CAerialSimulator::GetLine(istream &in, char *buf, const int bufSize) const
{
	while(in.good())
	{
		in.getline(buf, bufSize);
		if(buf[0]!='%') break;
	}
}

/*******************************************************************
 * Description:	Read the configuration file for the aerial section.
 ********************************************************************/
void CAerialSimulator::ReadConfig()
{
	// Open the configuration file to read.
	ifstream fileInput(m_CFG_File.c_str());

	if (!fileInput)
	{
		cout <<"Input 폴더의 Configuration file 파일을 확인할 수 없습니다" << endl;
		cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
		_getch();			
		exit(0);
	}
	
	const int BUFFER_SIZE = 200;
	char buffer[BUFFER_SIZE+1] = {0};

	// Read the time period to capture a new image(second)
	//GetLine(fileInput, buffer, BUFFER_SIZE);
	//sscanf(buffer, "%lf", &CAPTURE_PERIOD);
	CAPTURE_PERIOD=1.0;

	// Read the number of images in the sequence.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	sscanf_s(buffer, "%d", &TOTAL_IMAGE);

	// Read the EO file name.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	istrstream(buffer) >> m_EO_File;

	// Read the image directory.
	GetLine(fileInput, buffer, BUFFER_SIZE);
	istrstream(buffer) >> m_img_directory;

	// Read the image file names.
	for (int i=0; i<TOTAL_IMAGE; i++)
	{
		string img_file;

		GetLine(fileInput, buffer, BUFFER_SIZE);
		istrstream(buffer) >> img_file;	
		
		// Store the image file name into the deque.
		m_deq_IMG_File.push_back(img_file);
	}


	fileInput.close();

	char strFolderPath[] = { "Results\\Matching_Results" };
     
    int nResult = _mkdir( strFolderPath );
 
    if( nResult == 0 )
    {
        cout <<"Matching_Results 폴더 생성" << endl;
    }
    else if( nResult == -1 )
	{
        //cout << "error no : " <<  errno << endl;;
    }  

	FILE * fEpoch = _fsopen("Results\\Matching_Results\\TiePoint_fcs.txt", "w", _SH_DENYWR);
	fclose(fEpoch);

	fEpoch = _fsopen("Results\\Matching_Results\\TiePoint_fcs.txt", "a+", _SH_DENYWR);
	fprintf(fEpoch, "No. Images : %d\n", TOTAL_IMAGE); // To compliance with MATLAB
	fclose(fEpoch);

	fEpoch = _fsopen("Results\\Matching_Results\\TiePoint_ics.txt", "w", _SH_DENYWR);
	fclose(fEpoch);

	fEpoch = _fsopen("Results\\Matching_Results\\TiePoint_ics.txt", "a+", _SH_DENYWR);
	fprintf(fEpoch, "No. Images : %d\n", TOTAL_IMAGE); // To compliance with MATLAB
	fclose(fEpoch);	


}

/********************************************************************
 * Description:	Read the EO file
 ********************************************************************/
void CAerialSimulator::ReadEOFile()
{
	EO		eo;

	// If the file exist, read the EO data from the file.
	FILE * fEO = _fsopen(m_EO_File.c_str(), "r", _SH_DENYWR);

	if (!fEO)
	{
		cout <<"Input 폴더의 EO 파일을 확인할 수 없습니다" << endl;
		cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
		_getch();			
		exit(0);
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
	
		m_deq_EO.push_back(eo);
	}

	fclose(fEO);
}
