// Array.h: interface for the Array class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_ARRAY_H__FD29E191_7B91_11D5_9447_00105A29FE9E__INCLUDED_)
#define AFX_ARRAY_H__FD29E191_7B91_11D5_9447_00105A29FE9E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <vector>

using namespace std;

template <class T>
class Array : public vector<T>
{
protected:
	int row;
	int col;
public:
	// Indexing
	int Idx(const int r, const int c) const { return(r*col+c); }

	// Set Attributes
	void SetDim(const int _row, const int _col) { row=_row; col=_col; resize(row*col); }
	void Set(const T& val) { for(size_t n=0; n<size(); n++) (*this)[n]=val; }
	void Set(const int r, const int c, const T& val) { (*this)[Idx(r,c)]=val; }
	void Set(const int n, const T& val) { (*this)[n]=val; }

	// Get Attributes
	void Dim(int &r, int &c) const { r=row; c=col; }
	int NumCol() const { return col; }
	int NumRow() const { return row; }
	const T& Get(const int r, const int c) const { return (*this)[Idx(r,c)]; }
	const T& Get(const int n) const { return (*this)[n]; }
};

#endif // !defined(AFX_ARRAY_H__FD29E191_7B91_11D5_9447_00105A29FE9E__INCLUDED_)
