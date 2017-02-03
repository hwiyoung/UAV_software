/********************************************************
 * Project:			Automated AT
 * Last updated:	17 Jaunary 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#include <stdio.h>
#include <windows.h>
#include "SimAT.h"				// To use the CSimAT object.
#include "AerialSimulator.h"	// To use the CAerialSimulator object.
#include "KLTTracker.h"			// To use the CKLTTracker object.

void ShowUsage();

void main(int argc, char *argv[])
{	
	cout << "KLT Matching start" << endl;
	argv[1]="..\\Input\\AerialConfig.txt";
	argv[2]="..\\Input\\KLTConfigModel.txt";
	argv[3]="..\\Input\\SimATConfig.txt";

	// Release ¿ë
	//argv[1]="..\\Input\\AerialConfig.txt";
	//argv[2]="..\\Input\\KLTConfigModel.txt";
	//argv[3]="..\\Input\\SimATConfig.txt";
	
	// 2. Initialize the three main components to build the Automated AT system.
	try
	{
		cout << "Aerial Section" << endl;
		CAerialSimulator AerialSimulator(argv[1]);
		AerialSimulator.StartAerialSimulator();

		cout << "Image Matching" << endl;
		CKLTTracker KLT(argv[2]);
		KLT.Startup();

		cout << "Simultaneous AT " << endl;
		CSimAT SimAT(argv[3]);
		SimAT.Initialization();

		cout << "Register the clients of the Aerial simulator." << endl;
		AerialSimulator.RegisterKLTClient(&KLT);
		AerialSimulator.RegisterSimAT(&SimAT);

		cout << "Register the KLT object to the Simultaneous AT object. " << endl;
		SimAT.RegisterKLT(&KLT);

		cout << "Start the Aerial section to capture images. " << endl;
		AerialSimulator.RunImageExposure();
	
		cout << "Finally, print out the result." << endl;
		SimAT.PrintOutResults();
	}
	catch(...)
	{
		
	}

}

/*
 * Print out the help message.
 */
void ShowUsage()
{
	printf("Automated Aerial Triangulation System\n");
	printf("     Usage:     AutomatedAT <Aerial-config-file> <KLT-config-file> <AT-config-file>\n");
	printf("   Example:     AutomatedAT AerialConfig.txt KLTConfig.txt SimATConfig.txt\n");
}
