#pragma once

#include "../KLT/AerialSimulator.h"	// To use the CAerialSimulator object
#include "../KLT/KLTTracker.h"			// To use the CKLTTracker object
#include <iostream>
#include <algorithm>
#include <cv.h>
#include <cxcore.h>
#include <highgui.h>

typedef struct tagData
{
	string imagename;
	double eo[6];
}Data;

void ortho();
istream& ReadInputFile(istream& in, vector<Data>& vec);
IplImage* orthophoto_m(IplImage* srcImage, int* pixel_cnt, double pixel_size,
	double focal_length, double* EO, double gsd, double ground_height);
CvMat* rot3D(double* EO);
void vertex_g(int* pixel_cnt, double pixel_size, double focal_length, double* EO, CvMat* R, double ground_height, double* g_c);
CvMat* xy_g_min(double* EO, CvMat* R, CvMat* xy_image, double ground_height);
CvMatND* dem_m(double* g_c, int row_s, int col_s, double gsd, double ground_height);
CvMatND* image_coordinate(CvMatND* dem, CvMat* R, double* EO, double focal_length, double pixel_size, int* pixel_cnt);
IplImage* pixel_color(int* pixel_cnt, CvMatND* xy_fcs, IplImage* srcImage);
double getMax(double* n, int size);
double getMin(double* n, int size);