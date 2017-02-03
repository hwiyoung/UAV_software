/********************************************************
 * Project:			Automated AT
 * Last updated:	17 Jaunary 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _DEFINITION_H_
#define	_DEFINITION_H_


#pragma once
#pragma warning (disable:4786)

#include <deque>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <bitset>
#include <string>
#include <windows.h>
#include <cv.h>			// To use the CvPoint2D32f structure

using namespace std;

/*
 * struct IP
 * This structure is to store the image point ID and coordinate.
 */
struct IP
{
	int		iImg;
	int		iObj;
	double	dX;
	double	dY;

	IP()
	{
		iImg= -1;
		iObj= -1;
		dX	= 0.0;
		dY	= 0.0;
	}

	IP(int imgID, int objID, double x, double y)
	{
		iImg= imgID;
		iObj= objID;
		dX	= x;
		dY	= y;
	}
};

/*
 * struct EO
 * This structure is to store the exterior orientation parameters for each image.
 */
struct EO
{
	int		iImg;
	double	dXc;
	double	dYc;
	double	dZc;
	double	dOmega;	// Roll
	double	dPhi;	// Nick
	double	dKappa;	// Kappa

	EO()
	{
		iImg	= -1;
		dXc		= 0.0;
		dYc		= 0.0;
		dZc		= 0.0;
		dOmega	= 0.0;
		dPhi	= 0.0;
		dKappa	= 0.0;
	}

	EO(int imgID, double Xc, double Yc, double Zc, double omega, double phi, double kappa)
	{
		iImg	= -1;
		dXc		= Xc;
		dYc		= Yc;
		dZc		= Zc;
		dOmega	= omega;
		dPhi	= phi;
		dKappa	= kappa;
	}
};

/*
 * struct IO
 * This structure is to store the interior orientation parameters.
 */
struct IO
{
	double	dPPX;	// Principal point coordinate - x
	double	dPPY;	// Principal point coordinate - y
	double	dF;		// Focal length

	IO()
	{
		dPPX	= 0.0;
		dPPY	= 0.0;
		dF		= 0.0;
	}

	IO(double PPX, double PPY, double F)
	{
		dPPX	= PPX;
		dPPY	= PPY;
		dF		= F;
	}
};

/*
 * struct PS
 * This structure is to store pixel size.
 */
struct PS
{
	double	dX;		// Pixel size in x-direction
	double	dY;		// Pixel size in y-direction

	PS()
	{
		dX	= 0.0;
		dY	= 0.0;
	}

	PS(double x, double y)
	{
		dX	= x;
		dY	= y;
	}
};

/*
 * struct IS
 * This structure is to store the image dimension.
 */
struct IS
{
	int	iWidth;		// The width of an image
	int	iHeight;	// The height of an image

	IS()
	{
		iWidth	= 0;
		iHeight	= 0;
	}

	IS(int w, int h)
	{
		iWidth	= w;
		iHeight	= h;
	}
};

/*
 * struct GP
 * This structure is to store the ground point.
 */
struct GP
{
	double	X;
	double	Y;
	double	Z;

	GP()
	{
		X	= 0.0;
		Y	= 0.0;
		Z	= 0.0;
	}

	GP(double x, double y, double z)
	{
		X	= x;
		Y	= y;
		Z	= z;
	}
};

/*
 * struct RDC
 * This structure is to store the radial distortion.
 */
struct RDC
{
	double	dK1;	// K1
	double	dK2;	// K2
	double	dK3;	// K3

	RDC()
	{
		dK1	= 0.0;
		dK2	= 0.0;
		dK3	= 0.0;
	}

	RDC(double K1, double K2, double K3)
	{
		dK1	= K1;
		dK2	= K2;
		dK3	= K3;
	}
};

/*
 * struct POINT_XY
 * This structure is to store pixel coordinate.
 */
struct POINT_XY
{
	int	x;		// x-coordinate
	int	y;		// y-coordinate

	POINT_XY()
	{
		x	= 0;
		y	= 0;
	}

	POINT_XY(int xx, int yy)
	{
		x	= xx;
		y	= yy;
	}
};

#define NUMBITS	(32)	// The number of bits per DWORD

/*
 * struct DWORD_BITS
 * This structure is to store DWORD to be used as array of bits.
 * The structure defines 4 DWORD. It limits to maximum of 128 bits.
 */
struct DWORD_BITS
{
	DWORD	dword0 : 32;	// Bit#0 - 31
	DWORD	dword1 : 32;	// Bit#32 - 63
	DWORD	dword2 : 32;	// Bit#64 - 95
	DWORD	dword3 : 32;	// Bit#96 - 127

	DWORD_BITS()
	{
		dword0	= 0x0;
		dword1	= 0x0;
		dword2	= 0x0;
		dword3	= 0x0;
	}

	DWORD_BITS(DWORD d0, DWORD d1, DWORD d2, DWORD d3)
	{
		dword0	= d0;
		dword1	= d1;
		dword2	= d2;
		dword3	= d3;
	}
};

/*
 * struct PAIR_INT
 * This structure is to pair of integer.
 */
struct PAIR_INT
{
	int elm1;
	int	elm2;

	PAIR_INT()
	{
		elm1	= 0;
		elm2	= 0;
	}

	PAIR_INT(int e1, int e2)
	{
		elm1	= e1;
		elm2	= e2;
	}
};

/*
 * struct TIE_POINT
 * This structure is to store the tie point.
 */
struct TIE_POINT
{
	int	pointID;
	int imageID;
	int regionID;
	CvPoint2D32f coord;
};

/*
 * struct SIM2D_TRANS_PARAM
 * This structure is to store the 2D similarity transform parameters.
 */
struct SIM2D_TRANS_PARAM
{
	float	scale;
	float	angle;
	float	tx;		// translation vector in x-axis
	float	ty;		// translation vector in y-axis

	SIM2D_TRANS_PARAM()
	{
		scale	= 0.0;
		angle	= 0.0;
		tx		= 0.0;	
		ty		= 0.0;	
	}

	SIM2D_TRANS_PARAM(float p1, float p2, float p3, float p4)
	{
		scale	= p1;
		angle	= p2;
		tx		= p3;	
		ty		= p4;	
	}
};

// Defined constant numbers for performing AT
#define		N_UNKNOWN_GP	(3)
#define		N_UNKNOWN_EO	(6)

typedef vector< EO > VEC_EO;		// Vector of EO parameters / EO structure
typedef vector< IP > VEC_IP;		// Vector of image points / IP structure
typedef vector< GP > VEC_GP;		// Vector of ground point / GP structure
typedef vector< int > VEC_INT;		// Vector of integers
typedef vector< bool > VEC_BOOL;	// Vector of bool
typedef vector< DWORD > VEC_DWORD;	// Vector of DWORD
typedef vector< double > VEC_DOUBLE;// Vector of double
typedef vector< DWORD_BITS > VEC_DWORD_BITS;	// Vector of array of bits / DWORD_BITS
typedef vector< TIE_POINT > VEC_TP;	// Vector of tie point / TIE_POINT structure
typedef vector< CvPoint2D32f > VEC_CVPOINT;	// Vector of CvPoint2D32f
typedef set< int > SET_INT;			// Set of integers
typedef map< int, int > MAP_II;		// Map of integer as the key and integer as the value
typedef bitset< NUMBITS >	BITS32;	// A 32-bit bitset
typedef deque< EO > DEQ_EO;			// Deque of EO parameters / EO structure
typedef deque< string > DEQ_STRING;	// Deque of std::string
typedef deque< double > DEQ_DOUBLE;	// Deque of double

// Comparing two floating point number
#define SMALL_NON_ZERO 1e-16 /* or something else small */
#define DBL_GT(X,Y) ( (X) - (Y) > SMALL_NON_ZERO ) /* X > Y */
#define DBL_LT(X,Y) ( (Y) - (X) > SMALL_NON_ZERO ) /* Y > X */
#define DBL_EQ(X,Y) ( !DBL_GT((X),(Y)) && !DBL_LT((X),(Y)) ) /* X == Y */

// Defined constant numbers for KLT tracker
#define		MAX_TRACKING_IMG	(2)		// Maximum number of images for tracking at a time.
#define		IMG_A				(0)		// Image A identifier (the first image)
#define		IMG_B				(1)		// Image B identifier (the consecutive image)

// Define a minim angle to decide whether or not warp image
#define		MIN_ANGLE_NO_WARP	(10)		// Minimum angle (degrees) to not warp image

// Detector type
enum DETECTOR_TYPE 
	{	KLT		= 1,
		FAST	= 2	
	};

// KLT tracking model
enum TRACKING_MODEL 
	{	TRANSLATION					= 1,
		TRANSLATION_ROTATED_IMG		= 2,	
		TRANSLATION_ROTATED_WIN		= 3,
		AFFINE						= 4
	};

// Timer constant values
#define TIME_UNITS_PER_SECOND	((LONGLONG)10000000)

// Matrix threshold condition
const double THRESHOLD_CND = 1e12;

#endif
