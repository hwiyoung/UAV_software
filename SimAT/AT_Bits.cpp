/********************************************************
 * Project:			Simultaneous AT
 * Last updated:	27 October 2010
 * Developer:		Supannee Tanathong
 ********************************************************/

#include "AT_Bits.h"
#include <iostream>

/********************************************************************
 * Description:	Construction
 ********************************************************************/
AT_CBits::AT_CBits()
{
	// Initialize the vector containing power of 2.
	m_vec_power2.push_back(1);				// 2^0
	m_vec_power2.push_back(2);				// 2^1
	m_vec_power2.push_back(4);				// 2^2
	m_vec_power2.push_back(8);				// 2^3
	m_vec_power2.push_back(16);				// 2^4
	m_vec_power2.push_back(32);				// 2^5
	m_vec_power2.push_back(64);				// 2^6
	m_vec_power2.push_back(128);			// 2^7
	m_vec_power2.push_back(256);			// 2^8
	m_vec_power2.push_back(512);			// 2^9
	m_vec_power2.push_back(1024);			// 2^10
	m_vec_power2.push_back(2048);			// 2^11
	m_vec_power2.push_back(4096);			// 2^12
	m_vec_power2.push_back(8192);			// 2^13
	m_vec_power2.push_back(16384);			// 2^14
	m_vec_power2.push_back(32768);			// 2^15
	m_vec_power2.push_back(65536);			// 2^16
	m_vec_power2.push_back(131072);			// 2^17
	m_vec_power2.push_back(262144);			// 2^18
	m_vec_power2.push_back(524288);			// 2^19
	m_vec_power2.push_back(1048576);		// 2^20
	m_vec_power2.push_back(2097152);		// 2^21
	m_vec_power2.push_back(4194304);		// 2^22
	m_vec_power2.push_back(8388608);		// 2^23
	m_vec_power2.push_back(16777216);		// 2^24
	m_vec_power2.push_back(33554432);		// 2^25
	m_vec_power2.push_back(67108864);		// 2^26
	m_vec_power2.push_back(134217728);		// 2^27
	m_vec_power2.push_back(268435456);		// 2^28
	m_vec_power2.push_back(536870912);		// 2^29
	m_vec_power2.push_back(1073741824);		// 2^30
	m_vec_power2.push_back(2147483648);		// 2^31

}

/********************************************************************
 * Description:	Return the value of 2 power x.
 ********************************************************************/
DWORD AT_CBits::Power(int x)
{
	return m_vec_power2[x];
}

/********************************************************************
 * Description:	Check if the specified bit is set to 1.
 ********************************************************************/
bool AT_CBits::IsBitSet(DWORD dword_bits, int x)
{
	return AT_BITS32(dword_bits).test(x);	
}

/********************************************************************
 * Description:	Set the specified bit to 1.
 ********************************************************************/
void AT_CBits::SetBit(AT_DWORD_BITS &dword_bits, int x)
{
	int	index_dword, index_bit;

	// Check if the number of bits is greater than the limitation (128)
	if (x >= NUMBITS*4)
	{
		cout << x << endl;
		cout << NUMBITS << endl;	
		cout << "Error! The number of images is over the limitation." << endl;
		exit(0);
	}

	// Get the index of DWORD bits: 0,1,2,3
	index_dword = x/NUMBITS;

	// Get the index of bits: 0,1,2,...,31
	index_bit = x%NUMBITS;

	switch(index_dword)
	{
		case 0:
			dword_bits.dword0 = dword_bits.dword0 | (1 << index_bit); break;
		case 1:
			dword_bits.dword1 = dword_bits.dword1 | (1 << index_bit); break;
		case 2:
			dword_bits.dword2 = dword_bits.dword2 | (1 << index_bit); break;
		case 3:
			dword_bits.dword3 = dword_bits.dword3 | (1 << index_bit); break;
		default:
			cout << "Error! The number of images is over the limitation." << endl;
			exit(0);
	}
}

/********************************************************************
 * Description:	Get the next image ID which is setted to 1 starting
 *				from ImageID = iStart to Image ID = iStop.
 *				This function is applied for only 1 DWORD.
 ********************************************************************/
int AT_CBits::GetNextBitSet(DWORD const &dword_bits, int iStart, int iStop)
{
	int i;

	for (i = iStart; i <= iStop; i++)
	{
		if (IsBitSet(dword_bits, i))
			return i;
	}
	
	return -1;
}

/********************************************************************
 * Description:	Get the next image ID which is setted to 1 starting
 *				from ImageID = iStart to Image ID = iStop.
 *				This function is applied to all 4 DWORDs.
 ********************************************************************/
int AT_CBits::GetNextBitSet(AT_DWORD_BITS const &dword_bits, int iStart, int iStop)
{
	int index_dword_start, index_dword_stop;
	int index_bit_start, index_bit_stop;
	int iresult;

	// Validate the passed-in parameters.
	if ( (iStart >= NUMBITS*4) || (iStop >= NUMBITS*4) || (iStart < 0) || (iStop < 0) || (iStart > iStop) )
	{
		return -1;
	}

	// Get the index of DWORD bits: 0,1,2,3
	index_dword_start	= iStart/NUMBITS;
	index_dword_stop	= iStop/NUMBITS;

	// Get the index of bits: 0,1,2,...,31
	index_bit_start		= iStart%NUMBITS;
	index_bit_stop		= iStop%NUMBITS;

	iresult = -1;

	switch(index_dword_start)
	{
	case 0:

		if (index_dword_stop == 0)
		{
			iresult = GetNextBitSet(dword_bits.dword0, index_bit_start, index_bit_stop);
			return iresult;
		} 
		else
		{
			iresult = GetNextBitSet(dword_bits.dword0, index_bit_start, NUMBITS-1);

			if (iresult != -1)
				return iresult;
			// else continue
		}

	case 1:

		if (index_dword_stop == 1)
		{
			iresult = GetNextBitSet(dword_bits.dword1, (index_dword_start==1)?index_bit_start:0, index_bit_stop);

			if (iresult != -1)
				return NUMBITS + iresult;

			return -1;
		} 
		else if (index_dword_stop > 1)
		{
			iresult = GetNextBitSet(dword_bits.dword1, (index_dword_start==1)?index_bit_start:0, NUMBITS-1);

			if (iresult != -1)
				return NUMBITS + iresult;
			// else continue
		}

	case 2:

		if (index_dword_stop == 2)
		{
			iresult = GetNextBitSet(dword_bits.dword2, (index_dword_start==2)?index_bit_start:0, index_bit_stop);

			if (iresult != -1)
				return NUMBITS*2 + iresult;

			return -1;
		} 
		else if (index_dword_stop > 2)
		{
			iresult = GetNextBitSet(dword_bits.dword2, (index_dword_start==2)?index_bit_start:0, NUMBITS-1);

			if (iresult != -1)
				return NUMBITS*2 + iresult;
			// else continue
		}

	case 3:

		if (index_dword_stop == 3)
		{
			iresult = GetNextBitSet(dword_bits.dword3, (index_dword_start==3)?index_bit_start:0, index_bit_stop);

			if (iresult != -1)
				return NUMBITS*3 + iresult;

			return -1;
		} 

	}

	return -1;
}