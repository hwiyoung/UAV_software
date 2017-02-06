#include <direct.h>
#include <fstream>
#include <vector>
#include <string>
#include <iostream>
#include <memory>
#include <conio.h>
#include "ortho.h"
using namespace std;

void ortho()
{
	cout << endl;
	cout << "*************************** " << endl;
	cout << "Orthoimage Generation Start " << endl;
	cout << "*************************** " << endl;

	_mkdir("Results\\Orthoimage_Results");

	IplImage *srcImage;

	double pixel_size = 6E-6;		// unit : m
	double focal_length = 0.035;	// unit : m
	double gsd = 0.0;

	string input_foldername = "Input_Image\\";
	string output_foldername = "Results\\Orthoimage_Results\\";
	string img_input, img_output, img_name;

	int pixel_cnt[2];

	double ground_height = 0.0;

	ifstream in("Results\\AT_Results\\Estimated_EO.txt");

	vector<Data> vec;
	if (!in.is_open())
	{
		cout << "EO 파일을 확인할 수 없습니다" << endl;
		cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
		_getch();
		exit(0);
	}
	ReadInputFile(in, vec);

	//ofstream eo_ortho("Results\\Orthoimage_Results\\eo_ortho.txt");
	//eo_ortho << "%%Image Name\tGSD(m)\tX\tY\trow\tcol\n" << endl;
	
	int siz = 0;
	siz = vec.size();
	for (int k = 0; k < siz; k++)
	{
		img_name = vec[k].imagename;
		img_input = input_foldername + img_name;
		cout << endl;
		cout << img_name << " 처리 시작" << endl;

		srcImage = cvLoadImage(img_input.c_str());
		
		pixel_cnt[0] = srcImage->height;	// height
		pixel_cnt[1] = srcImage->width;		// width

		if (srcImage == NULL) {
			cout << img_name << "에 해당하는 사진이 없습니다" << endl;
			cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
			_getch();
			exit(0);
		}

		if (vec[k].eo[2] < 0)
		{
			cout << "Z 값 : " << vec[k].eo[2] << endl;
			cout << k << "번째의 Estimate EO의 Z값이 0보다 작으니 확인하시기 바랍니다" << endl;
			cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
			_getch();
			exit(0);
		}

		double EO[6] = { vec[k].eo[0], vec[k].eo[1], vec[k].eo[2], vec[k].eo[3], vec[k].eo[4], vec[k].eo[5] };
		gsd = (pixel_size * EO[2]) / focal_length;

		IplImage* dstImage;
		dstImage = orthophoto_m(srcImage, pixel_cnt, pixel_size, focal_length, EO, gsd, ground_height);

		img_output = output_foldername + img_name;
		cvSaveImage(img_output.c_str(), dstImage);
		
		//cvReleaseImage(&srcImage);
		cvReleaseImage(&dstImage);

		cout << img_name << "저장 완료" << endl;
	}	
}

IplImage* orthophoto_m(IplImage* srcImage, int* pixel_cnt, double pixel_size,
	double focal_length, double* EO, double gsd, double ground_height)
{
	// rotation matrix
	CvMat* R = rot3D(EO);
	/*R = cvCreateMat(3, 3, CV_32FC1);
	R = rot3D(EO);*/

	// ground coverage
	double g_c[4];
	vertex_g(pixel_cnt, pixel_size, focal_length, EO, R, ground_height, g_c);

	double col_s = (g_c[1] - g_c[0]) / gsd;
	double row_s = (g_c[3] - g_c[2]) / gsd;

	// construct grid point
	//CvMatND* dem = dem_m(g_c, ceil(row_s), ceil(col_s), gsd, ground_height);
	CvMatND* dem = dem_m(g_c, floor(row_s), floor(col_s), gsd, ground_height);

	// re-calculate image coordinate
	CvMatND* xy_fcs = image_coordinate(dem, R, EO, focal_length, pixel_size, pixel_cnt);

	// save image pixel value
	IplImage* dstImage = pixel_color(pixel_cnt, xy_fcs, srcImage);

	cvReleaseData(R);
	cvReleaseData(dem);
	cvReleaseData(xy_fcs);
	//cvReleaseImage(&srcImage);

	return dstImage;
}

CvMat* rot3D(double* EO)
{
	double	x;				// Angle		
	double	cos_x, sin_x;	// cos(Angle), sin(Angle)
	CvMat	*Rx, *Ry, *Rz, *R1, *R2;	// Rotational matrix


	// 1. Initialize the rotational matrix.
	Rx = cvCreateMat(3, 3, CV_32FC1);
	Ry = cvCreateMat(3, 3, CV_32FC1);
	Rz = cvCreateMat(3, 3, CV_32FC1);

	// 2. Manage R1 rotaion matrix
	R1 = cvCreateMat(3, 3, CV_32FC1);
	R2 = cvCreateMat(3, 3, CV_32FC1);

	// Construct rotational matrix for Image#1
	
	//		|	1			0		0		|
	// Rx =	|	0		 cos(Om)	sin(Om)	|
	//		|	0		 -sin(Om)	cos(Om)	|		

	x = EO[3];
	cos_x = cos(x);
	sin_x = sin(x);

	cvmSet(Rx, 0, 0, 1.0);
	cvmSet(Rx, 0, 1, 0.0);
	cvmSet(Rx, 0, 2, 0.0);
	cvmSet(Rx, 1, 0, 0.0);
	cvmSet(Rx, 1, 1, cos_x);
	cvmSet(Rx, 1, 2, sin_x);
	cvmSet(Rx, 2, 0, 0.0);
	cvmSet(Rx, 2, 1, -sin_x);
	cvmSet(Rx, 2, 2, cos_x);

	//		|	cos(Ph)		0		-sin(Ph)|
	// Ry =	|	  0			1		   0	|
	//		|	sin(Ph)		0		cos(Ph)	|

	x = EO[4];
	cos_x = cos(x);
	sin_x = sin(x);

	cvmSet(Ry, 0, 0, cos_x);
	cvmSet(Ry, 0, 1, 0.0);
	cvmSet(Ry, 0, 2, -sin_x);
	cvmSet(Ry, 1, 0, 0.0);
	cvmSet(Ry, 1, 1, 1.0);
	cvmSet(Ry, 1, 2, 0.0);
	cvmSet(Ry, 2, 0, sin_x);
	cvmSet(Ry, 2, 1, 0.0);
	cvmSet(Ry, 2, 2, cos_x);

	//		|	cos(Ka)		sin(Ka)		0	|
	// Rz =	|	-sin(Ka)	cos(Ka)		0	|
	//		|	  0			   0		1	|

	x = EO[5];	// Kappa
	cos_x = cos(x);
	sin_x = sin(x);

	cvmSet(Rz, 0, 0, cos_x);
	cvmSet(Rz, 0, 1, sin_x);
	cvmSet(Rz, 0, 2, 0.0);
	cvmSet(Rz, 1, 0, -sin_x);
	cvmSet(Rz, 1, 1, cos_x);
	cvmSet(Rz, 1, 2, 0.0);
	cvmSet(Rz, 2, 0, 0.0);
	cvmSet(Rz, 2, 1, 0.0);
	cvmSet(Rz, 2, 2, 1.0);

	// R = Rz * Ry * Rx
	cvMatMul(Rz, Ry, R1);	// R = Rz * Ry
	cvMatMul(R1, Rx, R2);	// R = (Rz * Ry) * Rx

	cvReleaseData(Rx);
	cvReleaseData(Ry);
	cvReleaseData(Rz);
	cvReleaseData(R1);

	return R2;
}

void vertex_g(int* pixel_cnt, double pixel_size, double focal_length, double* EO, CvMat* R, double ground_height, double* g_c)
{
	CvMat* image_xy;
	image_xy = cvCreateMat(4, 3, CV_32FC1);

	cvmSet(image_xy, 0, 0, pixel_cnt[1] * pixel_size / 2);
	cvmSet(image_xy, 0, 1, pixel_cnt[0] * pixel_size / 2);
	cvmSet(image_xy, 0, 2, -focal_length);
	cvmSet(image_xy, 1, 0, pixel_cnt[1] * pixel_size / 2);
	cvmSet(image_xy, 1, 1, -pixel_cnt[0] * pixel_size / 2);
	cvmSet(image_xy, 1, 2, -focal_length);
	cvmSet(image_xy, 2, 0, -pixel_cnt[1] * pixel_size / 2);
	cvmSet(image_xy, 2, 1, -pixel_cnt[0] * pixel_size / 2);
	cvmSet(image_xy, 2, 2, -focal_length);
	cvmSet(image_xy, 3, 0, -pixel_cnt[1] * pixel_size / 2);
	cvmSet(image_xy, 3, 1, pixel_cnt[0] * pixel_size / 2);
	cvmSet(image_xy, 3, 2, -focal_length);

	CvMat* xy_ground;
	xy_ground = cvCreateMat(4, 2, CV_32FC1);

	// each row matrix of image coordinate matrix
	CvMat* image_xy_row;
	image_xy_row = cvCreateMat(1, 3, CV_32FC1);

	// ground corner coordinates from image corner coorinates
	CvMat* xy_g_coord;
	xy_g_coord = cvCreateMat(1, 2, CV_32FC1);
	
	for (int i = 0; i < 4; i++) {
		for (int j = 0; j < 2; j++) {
			cvGetRow(image_xy, image_xy_row, i);
			xy_g_coord = xy_g_min(EO, R, image_xy_row, ground_height);
			cvmSet(xy_ground, i, j, cvmGet(xy_g_coord, 0, j));
		}
	}

	double edge_x[4] = { cvmGet(xy_ground, 0, 0), cvmGet(xy_ground, 1, 0) , cvmGet(xy_ground, 2, 0) ,cvmGet(xy_ground, 3, 0) };
	double edge_y[4] = { cvmGet(xy_ground, 0, 1), cvmGet(xy_ground, 1, 1) , cvmGet(xy_ground, 2, 1) ,cvmGet(xy_ground, 3, 1) };

	// ground coverage
	g_c[0] = getMin(edge_x, 4);
	g_c[1] = getMax(edge_x, 4);
	g_c[2] = getMin(edge_y, 4);
	g_c[3] = getMax(edge_y, 4);

	cvReleaseData(image_xy);
	cvReleaseData(xy_ground);
	cvReleaseData(image_xy_row);
	cvReleaseData(xy_g_coord);
}

CvMat* xy_g_min(double* EO, CvMat* R, CvMat* xy_image, double ground_height)
{
	// inverse matrix
	CvMat* inv_R;
	inv_R = cvCreateMat(3, 3, CV_32FC1);
	cvTranspose(R, inv_R);
	/*cvInvert(R, inv_R, CV_SVD_SYM);*/

	double scale = (ground_height - EO[2]) / (cvmGet(inv_R, 2, 0)*cvmGet(xy_image, 0, 0) + 
		cvmGet(inv_R, 2, 1)*cvmGet(xy_image, 0, 1) + cvmGet(inv_R, 2, 2)*cvmGet(xy_image, 0, 2));

	CvMat* sc;
	sc = cvCreateMat(3, 1, CV_32FC1);

	// transpose matrix
	CvMat* xy_image_t;
	xy_image_t = cvCreateMat(3, 1, CV_32FC1);
	cvTranspose(xy_image, xy_image_t);

	// sc = inv_R * xy_image' 3 by 1
	cvMatMul(inv_R, xy_image_t, sc);

	CvMat* xy_ground;
	xy_ground = cvCreateMat(1, 2, CV_32FC1);
	cvmSet(xy_ground, 0, 0, EO[0] + scale*cvmGet(sc, 0, 0));
	cvmSet(xy_ground, 0, 1, EO[1] + scale*cvmGet(sc, 1, 0));

	cvReleaseData(inv_R);
	cvReleaseData(sc);
	cvReleaseData(xy_image_t);

	return xy_ground;
}

CvMatND* dem_m(double* g_c, int row_s, int col_s, double gsd, double ground_height)
{
	// grid point
	int size[] = { row_s, col_s, 3 };
	CvMatND* dem;
	dem = cvCreateMatND(3, size, CV_32FC1);
	cvSetZero(dem);

	for (int i = 0; i < dem->dim[0].size; i++) {
		for (int j = 0; j < dem->dim[1].size; j++) {
			cvSetReal3D(dem, i, j, 0, g_c[3] - i*gsd);
			cvSetReal3D(dem, i, j, 1, g_c[0] + j*gsd);
			cvSetReal3D(dem, i, j, 2, ground_height);
		}
	}

	return dem;
}

CvMatND* image_coordinate(CvMatND* dem, CvMat* R, double* EO, double focal_length, double pixel_size, int* pixel_cnt)
{
	// fiducial point
	int size[] = { dem->dim[0].size, dem->dim[1].size, dem->dim[2].size };
	CvMatND* xy_fcs;
	xy_fcs = cvCreateMatND(3, size, CV_32FC1);
	cvSetZero(xy_fcs);

	double x_fcs = 0, y_fcs = 0;
	
	// ground point of each grid point
	CvMat* xyz_ground;
	xyz_ground = cvCreateMat(1, 3, CV_32FC1);

	double x_ground, y_ground, z_ground;

	//// eo convert double* to scalar
	//CvScalar eo_scalar;
	//eo_scalar = cvScalar(EO[0], EO[1], EO[2]);

	// substract matrix
	CvMat* subMat;
	subMat = cvCreateMat(1, 3, CV_32FC1);

	// transpose matrix
	CvMat* xyz_t;
	xyz_t = cvCreateMat(3, 1, CV_32FC1);

	// multiply matrix
	CvMat* coll_pcs;
	coll_pcs = cvCreateMat(3, 1, CV_32FC1);

	CvMat* xy_image;
	xy_image = cvCreateMat(3, 1, CV_32FC1);

	for (int i = 0; i < xy_fcs->dim[0].size; i++) {
		for (int j = 0; j < xy_fcs->dim[1].size; j++) {
			x_ground = cvGetReal3D(dem, i, j, 1);
			y_ground = cvGetReal3D(dem, i, j, 0);
			z_ground = cvGetReal3D(dem, i, j, 2);

			cvmSet(xyz_ground, 0, 0, x_ground - EO[0]);
			cvmSet(xyz_ground, 0, 1, y_ground - EO[1]);
			cvmSet(xyz_ground, 0, 2, z_ground - EO[2]);
				
			/*cvSubS(xyz_ground, eo_scalar, subMat, NULL);
			cvTranspose(subMat, xyz_t);*/
			cvTranspose(xyz_ground, xyz_t);

			// coll_pcs = R * xyz_t' 3 by 1
			cvMatMul(R, xyz_t, coll_pcs);

			double scale = - cvmGet(coll_pcs, 2, 0) / focal_length;

			cvmSet(xy_image, 0, 0, cvmGet(coll_pcs, 0, 0) / scale / pixel_size);
			cvmSet(xy_image, 1, 0, cvmGet(coll_pcs, 1, 0) / scale / pixel_size);
			cvmSet(xy_image, 2, 0, cvmGet(coll_pcs, 2, 0) / scale / pixel_size);

			x_fcs = pixel_cnt[0] / 2 + cvmGet(xy_image, 1, 0);
			y_fcs = pixel_cnt[1] / 2 - cvmGet(xy_image, 0, 0);

			cvSetReal3D(xy_fcs, i, j, 0, x_fcs);
			cvSetReal3D(xy_fcs, i, j, 1, y_fcs);
		}
	}
	cvReleaseData(subMat);
	cvReleaseData(xyz_t);
	cvReleaseData(coll_pcs);
	cvReleaseData(xy_image);
	cvReleaseData(xyz_ground);

	return xy_fcs;
}

IplImage* pixel_color(int* pixel_cnt, CvMatND* xy_fcs, IplImage* srcImage)
{
	int row = xy_fcs->dim[0].size;
	int col = xy_fcs->dim[1].size;

	 CvMat matHeader1, matHeader2, matHeader3, *pSrcMatB, *pSrcMatG, *pSrcMatR, *em_b, *em_g, *em_r;
	IplImage *dstB, *dstG, *dstR, *FB, *FG, *FR, *dstImage;

	// save pixel value
	em_b = cvCreateMat(pixel_cnt[0], pixel_cnt[1], CV_32FC1);		// cvCreateMat(row, col, *);
	em_g = cvCreateMat(pixel_cnt[0], pixel_cnt[1], CV_32FC1);
	em_r = cvCreateMat(pixel_cnt[0], pixel_cnt[1], CV_32FC1);

	// declare dstB, dstG, dstR
	dstB = cvCreateImage(cvSize(pixel_cnt[1], pixel_cnt[0]), IPL_DEPTH_8U, 1);		// cvSize(width, height);
	dstG = cvCreateImage(cvSize(pixel_cnt[1], pixel_cnt[0]), IPL_DEPTH_8U, 1);
	dstR = cvCreateImage(cvSize(pixel_cnt[1], pixel_cnt[0]), IPL_DEPTH_8U, 1);

	// consturct dstB, dstG, dstR
	cvSplit(srcImage, dstB, dstG, dstR, NULL);

	cvSetImageROI(dstB, cvRect(0, 0, pixel_cnt[1], pixel_cnt[0]));		// cvSize(leftOrigin, rightOrigin, width, height);
	cvSetImageROI(dstG, cvRect(0, 0, pixel_cnt[1], pixel_cnt[0]));
	cvSetImageROI(dstR, cvRect(0, 0, pixel_cnt[1], pixel_cnt[0]));

	// allign em_r, em_g, em_b memory
	pSrcMatB = cvGetMat(dstB, &matHeader1);
	pSrcMatG = cvGetMat(dstG, &matHeader2);
	pSrcMatR = cvGetMat(dstR, &matHeader3);

	// save pixel value to em_r, em_g, em_b
	for (int i = 0; i<pixel_cnt[0]; i++)
	{
		for (int j = 0; j<pixel_cnt[1]; j++)
		{
			cvmSet(em_b, i, j, cvGetReal2D(pSrcMatB, i, j));
			cvmSet(em_g, i, j, cvGetReal2D(pSrcMatG, i, j));
			cvmSet(em_r, i, j, cvGetReal2D(pSrcMatR, i, j));
		}
	}
	
	cvReleaseImage(&srcImage);
	cvReleaseData(pSrcMatR);
	cvReleaseData(pSrcMatG);
	cvReleaseData(pSrcMatB);

	cvReleaseImage(&dstB);
	cvReleaseImage(&dstG);
	cvReleaseImage(&dstR);

	FB = cvCreateImage(cvSize(col, row), IPL_DEPTH_8U, 1);
	FG = cvCreateImage(cvSize(col, row), IPL_DEPTH_8U, 1);
	FR = cvCreateImage(cvSize(col, row), IPL_DEPTH_8U, 1);

	cvSetImageROI(FB, cvRect(0, 0, col, row));		// cvRect(leftOrigin, rightOrigin, width, height);
	cvSetImageROI(FG, cvRect(0, 0, col, row));
	cvSetImageROI(FR, cvRect(0, 0, col, row));

	double outOfBoundsRow, outOfBoundsCol;
	for (int i = 0; i < row; i++) {
		for (int j = 0; j < col; j++) {
			outOfBoundsRow = cvGetReal3D(xy_fcs, row - (i + 1), col - (j + 1), 0);
			outOfBoundsCol = cvGetReal3D(xy_fcs, row - (i + 1), col - (j + 1), 1);

			if (outOfBoundsRow < 0 || outOfBoundsRow > pixel_cnt[0])
				continue;
			else if (outOfBoundsCol < 0 || outOfBoundsCol > pixel_cnt[1])
				continue;
			else {
				//cout << int(outOfBoundsRow) << "\t" << int(outOfBoundsCol) << endl;
				cvSetReal2D(FB, i, j, cvmGet(em_b, int(outOfBoundsRow), int(outOfBoundsCol)));
				cvSetReal2D(FG, i, j, cvmGet(em_g, int(outOfBoundsRow), int(outOfBoundsCol)));
				cvSetReal2D(FR, i, j, cvmGet(em_r, int(outOfBoundsRow), int(outOfBoundsCol)));
			}
		}
	}

	cvReleaseData(em_b);
	cvReleaseData(em_g);
	cvReleaseData(em_r);

	dstImage = cvCreateImage(cvSize(col, row), IPL_DEPTH_8U, 3);	// cvSize(width, height);
	cvSetZero(dstImage);
	cvMerge(FB, FG, FR, NULL, dstImage);
	
	cvReleaseImage(&FB);
	cvReleaseImage(&FG);
	cvReleaseImage(&FR);

	return dstImage;
}

//Get max value from n
double getMax(double* n, int size) {
	double max = n[0];

	for (int i = 1; i < size; i++)
		if (n[i] > max) max = n[i];

	return max;
}

//Get min value from n
double getMin(double* n, int size) {
	double min = n[0];

	for (int i = 1; i < size; i++)
		if (n[i] < min) min = n[i];

	return min;
}
