/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _AT_DEFINITION_H_
#define	_AT_DEFINITION_H_


#pragma once
#pragma warning (disable:4786)

#include <deque>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <bitset>
#include <windows.h>

using namespace std;

/*
 * struct IP
 * This structure is to store the image point ID and coordinate.
 */
struct AT_IP
{
	int		iImg;
	int		iObj;
	double	dX;
	double	dY;

	AT_IP()
	{
		iImg= -1;
		iObj= -1;
		dX	= 0.0;
		dY	= 0.0;
	}

	AT_IP(int imgID, int objID, double x, double y)
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
struct AT_EO
{
	int		iImg;
	double	dXc;
	double	dYc;
	double	dZc;
	double	dOmega;	// Roll
	double	dPhi;	// Nick
	double	dKappa;	// Kappa

	AT_EO()
	{
		iImg	= -1;
		dXc		= 0.0;
		dYc		= 0.0;
		dZc		= 0.0;
		dOmega	= 0.0;
		dPhi	= 0.0;
		dKappa	= 0.0;
	}

	AT_EO(int imgID, double Xc, double Yc, double Zc, double omega, double phi, double kappa)
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
struct AT_IO
{
	double	dPPX;	// Principal point coordinate - x
	double	dPPY;	// Principal point coordinate - y
	double	dF;		// Focal length

	AT_IO()
	{
		dPPX	= 0.0;
		dPPY	= 0.0;
		dF		= 0.0;
	}

	AT_IO(double PPX, double PPY, double F)
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
struct AT_PS
{
	double	dX;		// Pixel size in x-direction
	double	dY;		// Pixel size in y-direction

	AT_PS()
	{
		dX	= 0.0;
		dY	= 0.0;
	}

	AT_PS(double x, double y)
	{
		dX	= x;
		dY	= y;
	}
};

/*
 * struct IS
 * This structure is to store the image dimension.
 */
struct AT_IS
{
	int	iWidth;		// The width of an image
	int	iHeight;	// The height of an image

	AT_IS()
	{
		iWidth	= 0;
		iHeight	= 0;
	}

	AT_IS(int w, int h)
	{
		iWidth	= w;
		iHeight	= h;
	}
};

/*
 * struct GP
 * This structure is to store the ground point.
 */
struct AT_GP
{
	double	X;
	double	Y;
	double	Z;

	AT_GP()
	{
		X	= 0.0;
		Y	= 0.0;
		Z	= 0.0;
	}

	AT_GP(double x, double y, double z)
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
struct AT_RDC
{
	double	dK1;	// K1
	double	dK2;	// K2
	double	dK3;	// K3

	AT_RDC()
	{
		dK1	= 0.0;
		dK2	= 0.0;
		dK3	= 0.0;
	}

	AT_RDC(double K1, double K2, double K3)
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
struct AT_POINT_XY
{
	int	x;		// x-coordinate
	int	y;		// y-coordinate

	AT_POINT_XY()
	{
		x	= 0;
		y	= 0;
	}

	AT_POINT_XY(int xx, int yy)
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
struct AT_DWORD_BITS
{
	DWORD	dword0 : 32;	// Bit#0 - 31
	DWORD	dword1 : 32;	// Bit#32 - 63
	DWORD	dword2 : 32;	// Bit#64 - 95
	DWORD	dword3 : 32;	// Bit#96 - 127

	AT_DWORD_BITS()
	{
		dword0	= 0x0;
		dword1	= 0x0;
		dword2	= 0x0;
		dword3	= 0x0;

	}

	AT_DWORD_BITS(DWORD d0, DWORD d1, DWORD d2, DWORD d3)
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
struct AT_PAIR_INT
{
	int elm1;
	int	elm2;

	AT_PAIR_INT()
	{
		elm1	= 0;
		elm2	= 0;
	}

	AT_PAIR_INT(int e1, int e2)
	{
		elm1	= e1;
		elm2	= e2;
	}
};

// Defined constant numbers for performing AT
#define		N_UNKNOWN_GP	(3)
#define		N_UNKNOWN_EO	(6)

typedef vector< AT_EO > AT_VEC_EO;		// Vector of EO parameters / EO structure
typedef vector< AT_IP > AT_VEC_IP;		// Vector of image points / IP structure
typedef vector< AT_GP > AT_VEC_GP;		// Vector of ground point / GP structure
typedef vector< int > AT_VEC_INT;		// Vector of integers
typedef vector< bool > AT_VEC_BOOL;	// Vector of bool
typedef vector< DWORD > AT_VEC_DWORD;	// Vector of DWORD
typedef vector< double > AT_VEC_DOUBLE;// Vector of double
typedef vector< AT_DWORD_BITS> AT_VEC_DWORD_BITS;	// Vector of array of bits / DWORD_BITS
typedef set< int > AT_SET_INT;			// Set of integers
typedef map< int, int > AT_MAP_II;		// Map of integer as the key and integer as the value
typedef bitset<NUMBITS>	AT_BITS32;			// A 32-bit bitset
typedef deque< string > AT_DEQ_STRING;

// Comparing two floating point number
//#define SMALL_NON_ZERO 1e-16 /* or something else small */
//#define DBL_GT(X,Y) ( (X) - (Y) > SMALL_NON_ZERO ) /* X > Y */
//#define DBL_LT(X,Y) ( (Y) - (X) > SMALL_NON_ZERO ) /* Y > X */
//#define DBL_EQ(X,Y) ( !DBL_GT((X),(Y)) && !DBL_LT((X),(Y)) ) /* X == Y */

#define SMALL_NON_ZERO 0.001 /* or something else small */
#define BIG_NUM 10000
#define DBL_GT(X,Y) ( (X) - (Y) > SMALL_NON_ZERO ) /* X > Y */
#define DBL_LT(X) ( X < SMALL_NON_ZERO ) 
#define DBL_KT(X) ( X > BIG_NUM ) 
#define DBL_EQ(X,Y) ( !DBL_GT((X),(Y)) && !DBL_LT((X),(Y)) ) /* X == Y */

#endif
