/********************************************************
 * Project:			UAV Software
 * Last updated:	2014. 05. 28
 * Developer:		전의익
 ********************************************************/

#include <stdio.h>
#include <windows.h>
#include "Config/file_find.h"
#include "KLT/SimAT.h"				// To use the CSimAT object
#include "KLT/AerialSimulator.h"	// To use the CAerialSimulator object
#include "KLT/KLTTracker.h"			// To use the CKLTTracker object
#include "SimAT/AT_SimAT.h"
#include "Orthoimage/ortho_function.h"
#include <time.h>
#include "Ortho_test/ortho.h"

void main(int argc, char *argv[])
{	
	argv[1]="Input_80_170110\\AerialConfig.txt";
	argv[2]="Input_80_170110\\KLTConfigModel.txt";
	argv[3]="Input_80_170110\\SimATConfig.txt";

	clock_t start, end;
	start=clock();

	///*cout << "Config 파일 생성" << endl;
	//file_find();*/

	ofstream SaveFile("log.txt");

	//CAerialSimulator AerialSimulator(argv[1]);
	//AerialSimulator.StartAerialSimulator();
	//CKLTTracker KLT(argv[2]);
	//KLT.Startup();
	//CSimAT SimAT(argv[3]);
	//SimAT.Initialization();
	//AerialSimulator.RegisterKLTClient(&KLT);
	//AerialSimulator.RegisterSimAT(&SimAT);
	//SimAT.RegisterKLT(&KLT);
	//AerialSimulator.RunImageExposure();
	//SimAT.PrintOutResults();		
	//end=clock();
	//cout << "영상 매칭까지의 처리 시간은 "<< (end-start)/CLOCKS_PER_SEC << "초 입니다"<< endl;	
	//SaveFile << "영상 매칭까지의 처리 시간은 "<< (end-start)/CLOCKS_PER_SEC << "초 입니다"<< endl;

	AT_CSimAT AT_SimAT(argv[3]);
	AT_SimAT.RunSimAT();
	//_getch();
	end=clock();
	cout << "영상 지오레퍼런싱까지의 매칭까지의 처리 시간은 "<< (end-start)/CLOCKS_PER_SEC << "초 입니다"<< endl;
	SaveFile << "영상 지오레퍼런싱까지의 처리 시간은 "<< (end-start)/CLOCKS_PER_SEC << "초 입니다"<< endl;
	//_getch();

	//orthoimage_start();
	//end=clock();
	//cout << "정사영상 생성까지의 매칭까지의 처리 시간은 "<< (end-start)/CLOCKS_PER_SEC << "초 입니다"<< endl;
	////SaveFile << "정사영상 생성까지의 처리 시간은 "<< (end-start)/CLOCKS_PER_SEC << "초 입니다"<< endl;
	////SaveFile.close();
	//_getch();
	//exit(0);

	/*ortho();
	end=clock();
	cout << "정사영상 생성까지의 매칭까지의 처리 시간은 "<< (end-start)/CLOCKS_PER_SEC << "초 입니다"<< endl;
	_getch();
	exit(0);*/
}