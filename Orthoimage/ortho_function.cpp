#include "rt_nonfinite.h"
#include "ortho_function.h"
#include "ortho_function_emxutil.h"
#include "min.h"
#include "max.h"
#include "mldivide.h"
#include "cv.h"
#include "highgui.h"

#include <iostream> 
#include <fstream> 
#include <string> 
#include <cstdlib> 
#include <stdio.h>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <direct.h> //mkdir
#include <conio.h>
using namespace std; 

static int x_fcs[100000000], y_fcs[100000000];
static IplImage *FR, *FG, *FB;
static string output_foldername="Results\\Orthoimage_Results\\";

typedef struct tagData
{ 
	string imagename; 
	double eo[6]; 
}Data;

istream& ReadInputFile(istream& in, vector<Data>& vec)
{ 
     if( in ) { 
            vec.clear(); 
            Data d; 
            while( in >> d.imagename >> d.eo[0]>> d.eo[1] >> d.eo[2] >> d.eo[3] >> d.eo[4] >> d.eo[5] ) 
                vec.push_back(d); 
            in.clear(); 
     } 
     return in; 
} 

void ortho_function(const real_T EO[6], const double em_r[], const double em_g[], const double em_b[], string img_name)
{
	real_T R[3], dv0[9], dv1[9], dv2[9], dv3[9];
	int32_T i0, i1, loop_ub;
	real_T scale1,col_s, row_s;
	real_T sc1[3], b_R[9];
	real_T edge_x[4], edge_y[4];
	double pixel_size=6E-6;			// unit : m
	double pixel_cnt[2];
	double focal_length=0.035;		// unit : m
	pixel_cnt[0]=4000;				// unit : px - row
	pixel_cnt[1]=6000;				// unit : px - col
	double start_point_x, start_point_y;	// coord. of upper left
	int i,j;
	int tmp_r, tmp_g, tmp_b, xvalue, yvalue;
	real_T gsd=pixel_size*EO[2]/focal_length;
	
	real_T xyz_ground[3], b_xyz_ground[3];
	real_T XY1[2];real_T XY2[2];real_T XY3[2];real_T XY4[2];

	const int8_T iv0[3] = { 0, 0, 1 };
	const int8_T iv1[3] = { 0, 1, 0 };
	const int8_T iv2[3] = { 1, 0, 0 };
	const real_T dv4[3] = {  pixel_cnt[1]*pixel_size/2,  pixel_cnt[0]*pixel_size/2, -focal_length };
	const real_T dv5[3] = {  pixel_cnt[1]*pixel_size/2, -pixel_cnt[0]*pixel_size/2, -focal_length };
	const real_T dv6[3] = { -pixel_cnt[1]*pixel_size/2, -pixel_cnt[0]*pixel_size/2, -focal_length };
	const real_T dv7[3] = { -pixel_cnt[1]*pixel_size/2,  pixel_cnt[0]*pixel_size/2, -focal_length };

	//cout << "EO: "<< EO[0]<<" "<<EO[1]<<" "<<EO[2]<<" "<<EO[3]<<" "<<EO[4]<<" "<<EO[5]<< endl;
	R[0] = EO[3] * 3.14159265358979 / 180.0;
	R[1] = EO[4] * 3.14159265358979 / 180.0;
	R[2] = EO[5] * 3.14159265358979 / 180.0;

	//cout << "R[0]: " << R[0] << endl;
	//cout << "R[1]: " << R[1] << endl;
	//cout << "R[2]: " << R[2] << endl;

	// z-axis
	dv0[0] = cos(R[2]);	dv0[3] = sin(R[2]);	dv0[6] = 0.0;
	dv0[1] = -sin(R[2]);dv0[4] = cos(R[2]); dv0[7] = 0.0;

	// y-axis
	dv1[0] = cos(R[1]);	dv1[3] = 0.0;dv1[6] = -sin(R[1]);
	dv1[2] = sin(R[1]);dv1[5] = 0.0;dv1[8] = cos(R[1]);

	for (i0 = 0; i0 < 3; i0++) 
	{
		dv0[2 + 3 * i0] = (real_T)iv0[i0];
		dv1[1 + 3 * i0] = (real_T)iv1[i0];
	}	

	for (i0 = 0; i0 < 3; i0++) 
	{
		for (loop_ub = 0; loop_ub < 3; loop_ub++) 
		{
			dv2[i0 + 3 * loop_ub] = 0.0;
			for (i1 = 0; i1 < 3; i1++)
			{
				dv2[i0 + 3 * loop_ub] += dv0[i0 + 3 * i1] * dv1[i1 + 3 * loop_ub];
			}
		}
		dv3[3 * i0] = (real_T)iv2[i0];
	}

	// x-axis
	dv3[1] = 0.0;dv3[4] = cos(R[0]);dv3[7] = sin(R[0]);
	dv3[2] = 0.0;dv3[5] = -sin(R[0]);dv3[8] = cos(R[0]);

	for (i0 = 0; i0 < 3; i0++)
	{
		for (loop_ub = 0; loop_ub < 3; loop_ub++) 
		{
			b_R[i0 + 3 * loop_ub] = 0.0;
			for (i1 = 0; i1 < 3; i1++)
			{
				b_R[i0 + 3 * loop_ub] += dv2[i0 + 3 * i1] * dv3[i1 + 3 * loop_ub];
			}
		}
	}
	/*cout<<endl;
	cout << "역행렬"<<endl;
	cout<< b_R[0] << " " << b_R[1] << " " << b_R[2] << endl;
	cout<< b_R[3] << " " << b_R[4] << " " << b_R[5] << endl;
	cout<< b_R[6] << " " << b_R[7] << " " << b_R[8] << endl;
	cout<<endl;
*/

	scale1 = (0.0 - EO[2]) / (b_R[6] * dv4[0] + b_R[7] * dv4[1] + b_R[8] * -focal_length);
	mldivide(b_R, dv4, sc1);
	XY1[0] = EO[0] + scale1 * sc1[0];
	XY1[1] = EO[1] + scale1 * sc1[1];

	//cout << "sc1[0]: " << sc1[0] << " " <<sc1[1] << endl;
	//cout << "scale1: " <<scale1 << endl;
	//cout <<"XY1: "<<XY1[0]<< ", " << XY1[1]<< endl;
	//cout<<endl;

	scale1 = (0.0 - EO[2]) / ((b_R[6] * dv5[0]) + (b_R[7] * dv5[1]) + (b_R[8] * -focal_length));
	mldivide(b_R, dv5, sc1);
	XY2[0] = EO[0] + scale1 * sc1[0];
	XY2[1] = EO[1] + scale1 * sc1[1];

	//cout << "sc2[0]: " << sc1[0] << " " <<sc1[1] << endl;
	//cout << "scale2: " <<scale1 << endl;
	//cout <<"XY2: "<<XY2[0]<< ", " << XY2[1]<< endl;
	//cout<<endl;

	scale1 = (0.0 - EO[2]) / ((b_R[6] * dv6[0] + b_R[7] * dv6[1]) + b_R[8] * -focal_length);
	mldivide(b_R, dv6, sc1);
	XY3[0] = EO[0] + scale1 * sc1[0];
	XY3[1] = EO[1] + scale1 * sc1[1];

	//cout << "sc3[0]: " << sc1[0] << " " <<sc1[1] << endl;
	//cout << "scale3: " <<scale1 << endl;
	//cout <<"XY3: "<<XY3[0]<< ", " << XY3[1]<< endl;
	//cout<<endl;

	scale1 = (0.0 - EO[2]) / ((b_R[6] * dv7[0] + b_R[7] * dv7[1]) + b_R[8] * -focal_length);
	mldivide(b_R, dv7, sc1);
	XY4[0] = EO[0] + scale1 * sc1[0];
	XY4[1] = EO[1] + scale1 * sc1[1];

	//cout << "sc4[0]: " << sc1[0] << " " <<sc1[1] << endl;
	//cout << "scale4: " <<scale1 << endl;
	//cout <<"XY4: "<<XY4[0]<< ", " << XY4[1]<< endl;
	//cout<<endl;

	edge_x[0] = XY1[0];	edge_x[1] = XY2[0];	edge_x[2] = XY3[0];	edge_x[3] = XY4[0];
	edge_y[0] = XY1[1];	edge_y[1] = XY2[1];	edge_y[2] = XY3[1];	edge_y[3] = XY4[1];	


	gsd=pixel_size*EO[2]/focal_length;

	start_point_x=b_min(edge_x);
	start_point_y=b_max(edge_y);
	col_s = (b_max(edge_x) - start_point_x) / gsd;
	row_s = (start_point_y - b_min(edge_y)) / gsd;
	col_s=floor(col_s);row_s=floor(row_s);
	
	cout << "GSD : "<<gsd<<endl;
	cout << "격자수 : " <<(int)(col_s*row_s) <<endl;


	i0=(int32_T)row_s*col_s;

	R[0] = EO[0];R[1] = EO[1];R[2] = EO[2];

	for (i = 0; i <= (int)row_s-1; i++)
	{
		//cout << i << "th row(s)" << endl;
		for (j = 0; j <= (int)col_s-1; j++) 
		{			
			xyz_ground[0] = start_point_x + j*gsd;
			xyz_ground[1] = start_point_y - i*gsd;
			xyz_ground[2] = 0.0;

			for (i0 = 0; i0 < 3; i0++)	
			{
				b_xyz_ground[i0] = xyz_ground[i0] - R[i0];
			}
			for (i0 = 0; i0 < 3; i0++)
			{
				sc1[i0] = 0.0;
				for (loop_ub = 0; loop_ub < 3; loop_ub++)
				{
					sc1[i0] += b_R[i0 + 3 * loop_ub] * b_xyz_ground[loop_ub];
				}
			}
			scale1 = -sc1[2] / focal_length;
			for (i0 = 0; i0 < 3; i0++)
			{
				sc1[i0] = sc1[i0] / scale1 / pixel_size;
			}

			x_fcs[(int32_T)(j+i*col_s)] = 2000.0 + sc1[1];
			y_fcs[(int32_T)(j+i*col_s)] = 3000.0 - sc1[0];
		}
	}

	FR = cvCreateImage(cvSize((int)col_s, (int)row_s), IPL_DEPTH_8U, 1);
	FG = cvCreateImage(cvSize((int)col_s, (int)row_s), IPL_DEPTH_8U, 1);
	FB = cvCreateImage(cvSize((int)col_s, (int)row_s), IPL_DEPTH_8U, 1);
	cvSetImageROI(FR, cvRect(0,0,(int)col_s,(int)row_s));	
	cvSetImageROI(FG, cvRect(0,0,(int)col_s,(int)row_s));
	cvSetImageROI(FB, cvRect(0,0,(int)col_s,(int)row_s));
	cvSetZero(FR);cvSetZero(FG);cvSetZero(FB);	
	
	for (i = 0; i < (int)row_s; i++)
	{
		for (j = 0; j < (int)col_s; j++)
		{
			xvalue=x_fcs[(int)(row_s*col_s-j-i*col_s)];
			yvalue=y_fcs[(int)(row_s*col_s-j-i*col_s)];

			if ((xvalue<0.0) || (xvalue>pixel_cnt[0]) || (yvalue<0.0) || (yvalue>pixel_cnt[1]))
			{
				continue;
			}
			else 
			{
				//cout << xvalue << "\t" << yvalue << endl;
				tmp_b=(int)em_b[(int)(xvalue*pixel_cnt[1]+yvalue)];
				tmp_g=(int)em_g[(int)(xvalue*pixel_cnt[1]+yvalue)];
				tmp_r=(int)em_r[(int)(xvalue*pixel_cnt[1]+yvalue)];	

				cvSetReal2D(FR,i,j,tmp_r);
				cvSetReal2D(FG,i,j,tmp_g);
				cvSetReal2D(FB,i,j,tmp_b);
			}
		}		
	}
	
	IplImage *dstImage;
	dstImage=cvCreateImage(cvSize((int)col_s, (int)row_s), IPL_DEPTH_8U,3);
	cvSetZero(dstImage);
	cvMerge(FB, FG, FR, NULL, dstImage);
	
	string img=output_foldername + img_name ;


	cvSaveImage(img.c_str(),dstImage);

	cvReleaseData(dstImage);
	cvReleaseData(FR);
	cvReleaseData(FG);
	cvReleaseData(FB);
	cout << img << "저장 완료" << endl; 
}

void orthoimage_start()
{
	cout << endl;
	cout << "*************************** " << endl;
	cout << "Orthoimage Generation Start " << endl;
	cout << "*************************** " << endl;

	clock_t start, end;
	int i, j, k;
	int pixel_cnt[2]={6000, 4000};
	IplImage *srcImage,*dstB,*dstG,*dstR;
	static real_T em_r[24000000],em_g[24000000],em_b[24000000];
	string input_foldername="Input_Image\\";
	string img, img_name;
	_mkdir("Results\\Orthoimage_Results");
	
	start = clock();

	ifstream in("Results\\AT_Results\\Estimated_EO.txt"); 

    vector<Data> vec; 
    if( !in.is_open() ) 
	{     	
		cout <<"EO 파일을 확인할 수 없습니다" << endl;
		cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
		_getch();			
		exit(0);
    } 
    ReadInputFile(in, vec);
	
	int vecSize=0;
	vecSize=vec.size();
	for (k=0;k<vecSize;k++)
	{	
		img_name=vec[k].imagename;
		img=input_foldername + img_name;
		cout<<endl;
		cout << img_name << " 처리 시작" << endl;

		srcImage=cvLoadImage(img.c_str());
		
		if(srcImage==NULL){
			cout << img_name << "에 해당하는 사진이 없습니다" <<endl;
			cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
			_getch();	
			exit(0);
		}		
		if(vec[k].eo[2] < 0)
		{
			cout << "Z 값 : " << vec[k].eo[2] << endl;
			cout << k << "번째의 Estimate EO의 Z값이 0보다 작으니 확인하시기 바랍니다" << endl;
			cout << "프로그램을 종료하겠습니다. 아무키나 누르세요" << endl;
			_getch();			
			exit(0);

		}

		real_T EO[6]={vec[k].eo[0], vec[k].eo[1], vec[k].eo[2], vec[k].eo[3]*180/3.14, vec[k].eo[4]*180/3.14, vec[k].eo[5]*180/3.14};

		CvMat matHeader1,matHeader2,matHeader3, *pSrcMatB,*pSrcMatG,*pSrcMatR;			
		
		dstB=cvCreateImage(cvSize(pixel_cnt[0],pixel_cnt[1]), IPL_DEPTH_8U, 1);
		dstG=cvCreateImage(cvSize(pixel_cnt[0],pixel_cnt[1]), IPL_DEPTH_8U, 1);
		dstR=cvCreateImage(cvSize(pixel_cnt[0],pixel_cnt[1]), IPL_DEPTH_8U, 1);
		cout<<"dstB,G,R 선언"<<endl;

		cvSplit(srcImage, dstB, dstG, dstR, NULL);
		cvSetImageROI(dstB, cvRect(0,0,pixel_cnt[0],pixel_cnt[1]));	
		cvSetImageROI(dstG, cvRect(0,0,pixel_cnt[0],pixel_cnt[1]));
		cvSetImageROI(dstR, cvRect(0,0,pixel_cnt[0],pixel_cnt[1]));

		cout << "dstB,G,R 생성"<<endl;
		pSrcMatB=cvGetMat(dstB,&matHeader1);
		pSrcMatG=cvGetMat(dstG,&matHeader2);
		pSrcMatR=cvGetMat(dstR,&matHeader3);

		memset((void *)&em_r[0], 0, 24000000U * sizeof(real_T));
		memset((void *)&em_g[0], 0, 24000000U * sizeof(real_T));
		memset((void *)&em_b[0], 0, 24000000U * sizeof(real_T));	
		
		cout<<"em_r,g,b 메모리 할당"<<endl;

		for(i=0;i<pixel_cnt[1];i++)
		{
			for(j=0;j<pixel_cnt[0];j++)
			{
				em_r[j+i*pixel_cnt[0]]=cvGetReal2D(pSrcMatR,i,j);
				em_g[j+i*pixel_cnt[0]]=cvGetReal2D(pSrcMatG,i,j);
				em_b[j+i*pixel_cnt[0]]=cvGetReal2D(pSrcMatB,i,j);
			}
		}

		cout<<"em_r,g,b에 화소값 저장"<<endl;
		cvReleaseImage(&srcImage);
		cvReleaseData(pSrcMatR);
		cvReleaseData(pSrcMatG);
		cvReleaseData(pSrcMatB);
	
		cvReleaseData(dstB);
		cvReleaseData(dstG);
		cvReleaseData(dstR);

		cout << k+1 <<"번 째 ortho_function 진입"<<endl;
		ortho_function(EO, em_r, em_g, em_b, img_name);

		end=clock();

		printf("처리 시간 = %.2f\n",(double)(end-start)/CLOCKS_PER_SEC); 	
	}
	
}