/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#ifndef _AT_BITS_H_
#define	_AT_BITS_H_


#pragma once

#include "AT_Definition.h"
#include <windows.h>
#include <math.h>

// This class is worked based on the DWORD_BITS struct.

class AT_CBits
{
public:
	// Constructor
	AT_CBits();

	// Return the value of 2 power x.
	DWORD Power(int x);

	// Check if the specified bit is set to 1.
	bool IsBitSet(DWORD dword_bits, int x);

	// Set the specified bit to 1.
	void SetBit(AT_DWORD_BITS &dword_bits, int x);

	// Get the next bit which is setted to 1 based on 'DWORD'.
	int GetNextBitSet(DWORD const &dword_bits, int iStart, int iStop);

	// Get the next bit which is setted to 1 based on 'DWORD_BITS' struct.
	int GetNextBitSet(AT_DWORD_BITS const &dword_bits, int iStart, int iStop);

	// Variables declaration.

	AT_VEC_DWORD	m_vec_power2;	// Vector of power of 2

};


#endif // _BITS_H_
