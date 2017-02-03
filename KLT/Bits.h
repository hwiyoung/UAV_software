/********************************************************
 * Project:			Automated AT
 * Last updated:	17 Jaunary 2011
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _BITS_H_
#define	_BITS_H_


#pragma once

#include "Definition.h"
#include <windows.h>
#include <math.h>

// This class is worked based on the DWORD_BITS struct.

class CBits
{
public:
	// Constructor
	CBits();

	// Return the value of 2 power x.
	DWORD Power(int x);

	// Check if the specified bit is set to 1.
	bool IsBitSet(DWORD dword_bits, int x);

	// Set the specified bit to 1.
	void SetBit(DWORD_BITS &dword_bits, int x);

	// Get the next bit which is setted to 1 based on 'DWORD'.
	int GetNextBitSet(DWORD const &dword_bits, int iStart, int iStop);

	// Get the next bit which is setted to 1 based on 'DWORD_BITS' struct.
	int GetNextBitSet(DWORD_BITS const &dword_bits, int iStart, int iStop);

	// Variables declaration.

	VEC_DWORD	m_vec_power2;	// Vector of power of 2

};


#endif // _BITS_H_
