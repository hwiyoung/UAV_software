/********************************************************
 * Project:			Config  父甸扁
 * Last updated:	2014. 05. 28
 * Developer:		傈狼劳
 ********************************************************/

#include <vector>
#include <stdio.h>
#include <io.h>
#include <conio.h>
#include <iostream> 
#include <fstream>

using namespace std;

void AerialConfig(vector<string> imageList);
void SimATConfig();
void KLTConfig();

void file_find()
{
	int i=0;
		
    _finddata_t fd;
    long handle;
    int result = 1;
    handle = _findfirst("Input_image\\*.jpg", &fd);
 
    if (handle == -1)
    {
        printf("There were no files.\n");
        return;
    }
 
    while (result != -1)
    {
        printf("File: %s\n", fd.name);
		result = _findnext(handle, &fd);

		i=i+1;
    }
	
	result = 1;
	handle = _findfirst("Input_image\\*.jpg", &fd); 

	vector<string> imageList(i); 
	i=0;

	while (result != -1)
    {
		imageList[i]=string(fd.name);
        result = _findnext(handle, &fd);	
		i=i+1;

    }
    _findclose(handle);

	AerialConfig(imageList);
	SimATConfig();
	KLTConfig();
    
}
void AerialConfig(vector<string> imageList)
{
	cout << "Aerial Config 颇老 积己" << endl;
	ofstream config;
	config.open("Input\\AerialConfig.txt");


	config << "% Number of images in the sequence"<< endl;
	config << imageList.size() << endl;

	config << "% EO File" << endl;	
	config << "Input\\EO.txt" <<endl;
	config << "% Image directory"<< endl;
	config << "Input_Image\\ " << endl;
	config << "% List of image file names" << endl;

	int imageSizeRow=0;
	int imageSizeCol=0;
	imageSizeRow = imageList.size();
	imageSizeCol = imageList[0].size();

	for (int i=0; i<imageSizeRow; i++ )
	{		
		for (int j=0; j<imageSizeCol; j++)
		{
			config << imageList[i][j];
		}
		config << endl;
	}		

	config.close();
}

void SimATConfig()
{
	cout << "SimAT Config 颇老 积己" << endl;
	ofstream config;
	config.open("Input\\SimATConfig.txt");

	config << "% IO File"<< endl;
	config << "Input\\IO.txt" <<endl;

	config << "% EO File"<< endl;
	config << "Input\\EO.txt" <<endl;

	config << "% Tiepoint File" << endl;
	config << "Results\\Matching_Results\\TiePoint_fcs_AT.txt" << endl;

	config << "% Maximum number of iterations for AT"<< endl;
	config << 10 << endl;

	config << "% Delta: stop criteria " << endl;
	config << 1e-3 << endl;

	config << "% std_IP: standard deviation of the image point errors (unit: mm) " << endl;
	config << 0.00477 << endl;
	
	config << "% std_GPS: standard deviation of the GPS errors (unit: m) " << endl;
	config << 1.5 << endl;
	
	config << "% std_INS: standard deviation of the INS point errors (unit: rad) " << endl;
	config << 0.5*3.14/180 << endl;

	config.close();
}



void KLTConfig()
{
	cout << "kLT Config 颇老 积己" << endl;
	ofstream config;
	config.open("Input\\KLTConfigModel.txt");

	config << "% IO File " << endl;
	config << "Input\\IO.txt"<<endl;
	config <<"% Image directory" << endl;
	config <<"Input_Image\\" << endl;
	config << "% SCALE FACTOR: 1(100%) or 2(50%) or 4(25%) is recommended" << endl; 
	config << 4 << endl;
	config <<  "% Rotational matrix order - 1 : Rx*Ry*Rz = 1 or 2 : Rz*Ry*Rx " << endl;
	config << 2 << endl; 
	config << "% Rotational matrix type: 1 or 2" << endl;
	config << "% Type 1" << endl;
	config << "% 	Rx = [1 0 0; 0 cos(om) -sin(om); 0 sin(om) cos(om)] " << endl;
	config << "% 	Ry = [cos(ph) 0 sin(ph); 0 1 0; -sin(ph) 0 cos(ph)]" << endl;
	config << "% 	Rz = [cos(kp) -sin(kp) 0; sin(kp) cos(kp) 0; 0 0 1]" << endl;
	config << "% Type 2 " << endl;
	config << "%	Rx = [1 0 0; 0 cos(om) sin(om); 0 -sin(om) cos(om)]" << endl;
	config << "%	Ry = [cos(ph) 0 -sin(ph); 0 1 0; sin(ph) 0 cos(ph)]" << endl;
	config << "%	Rz = [cos(kp) sin(kp) 0; -sin(kp) cos(kp) 0; 0 0 1]" << endl;
	config << 2 << endl;
	config << "% MAX_CORNERS: Limit the number of points to track" << endl;
	config << 4 << endl;
	config << "% MAX_BLOCK_TP: Limit the maximum number of tie point per block (0 for all possible)" << endl;
	config << 0 << endl;
	config << "% TRACKER: 1 for KLT/ 2 for FAST" << endl;
	config << 1 << endl;
	config << "% KLT::TRACKING_MODEL: Tracking model for KLT tracker" << endl;
	config << "% 1 - Translation" << endl;
	config << "% 2 - Translation with rotated image" << endl;
	config << "% 3 - Translation with rotated window" << endl;
	config << "% 4 - Affine" << endl;
	config << 1 << endl;
	config << "% KLT::MINIMUM_EIGENVALUE: Minimum eigen value" << endl;
	config << 0.20 << endl;
	config << "% KLT/FAST::MIN_DISTANCE: Minimum distance between adjacent good features point" << endl;
	config << 30.0 << endl;
	config << "% KLT::BLOCK_AUTOCORR: Window size of autocorrelation" << endl;
	config << 15 << endl;
	config << "% KLT::WIN_SIZE_SUBPIXEL: Window size for subpixel" << endl;
	config << 7 << endl;
	config << "% KLT::WIN_SIZE_TRACK: Window size for tracking" << endl;
	config << 25 << endl;
	config << "% KLT/FAST::WIN_SIZE_CCC: Window size for CCC" << endl;
	config << 15 << endl;
	config << "% KLT::SUBPATTERN: Pattern of the overlapping P x P e.g. 5x5, 4x4, 3x3" << endl;
	config << 5 << endl;
	config << "% KLT::CRI_NUM_ITERATION: Stop Criteria - iteration" << endl;
	config << 1000 << endl;
	config << "% KLT::CRI_EPSILON_LIMIT: Stop Criteria - distance" << endl;
	config << 0.001 << endl;
	config << "% KLT::DEPTH_LEVELS: Pyramid depth level" << endl;
	config << 5 << endl;
	config << "% KLT::FEATURE_ERROR: Tracking error" << endl;
	config << 10000 << endl;
	config << "% FAST::FAST_THRESHOLD: Threshold for FAST feature extractor" << endl;
	config << 30 << endl;
	config << "% FAST::NON_MAX_SUPPRESSION: 1 for non-maximum suppression, 0 otherwise" << endl;
	config << 1 << endl;
	config << "% FAST::NUM_FAST_PIXEL: between 9-12" << endl;
	config << 9 << endl;
	config << "% Minimum angle to rotate image (degree)" << endl;
	config << 15.0 << endl;
	config << "% Perform Calculating Initial Guessed Tie-points (1 = Yes, 0 = No)" << endl;
	config << 0 << endl;
	config << "% Consider out-of-scope tie-points (0 = remove permanently, 1 = previous pos)" << endl;
	config << 0 << endl;
	config << "% Perform Cross Corrrelation Coefficient checking (1 = Yes, 0 = No)" << endl;
	config << 0 << endl;
	config << "% CCC Accepting threshold (0.0 - 1.0)" << endl;
	config << 0.90 << endl;
	config << "% Perform outlier removal by epipolar line (1 = Yes, 0 = No)" << endl;
	config << 0 << endl;
	config << "% Epipolar distance reject ratio (Mean + R*Std)" << endl;
	config << 2.0 << endl;
	config << "% Perform outlier removal by affine (1 = Yes, 0 = No)" << endl;
	config << 1 << endl;
	config << "% Affine distance reject ratio (Mean + R*Std)" << endl;
	config << 1.5 << endl;
	config << "% ZA: Average terrain elevation" << endl;
	config << 0 << endl;
	config << "% Debug out intermediate tie points into TP.txt file. 1:Yes, 0:No" << endl;
	config << 1 << endl;

	config.close();
}


