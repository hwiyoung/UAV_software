// Point.h: interface for the Point class.
//
//////////////////////////////////////////////////////////////////////

#if !defined(AFX_POINT_H__D4C0D00E_6FC6_11D5_9446_00105A29FE9E__INCLUDED_)
#define AFX_POINT_H__D4C0D00E_6FC6_11D5_9446_00105A29FE9E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include <cmath>
#include <iostream>
#include <ostream>

using namespace std;

class Matrix;

class Point  
{
protected:
	int id;
	double x;
	double y;
	double z;

public:
	// Constructors
	Point() { Set(-1, 0.0, 0.0, 0.0); }
	Point(const double _x, const double _y, const double _z) { Set(-1, _x, _y, _z); }
	Point(const Point &pt) { Set(pt.id, pt.x, pt.y, pt.z); }
	Point(const int _id, const double _x, const double _y, const double _z) { Set(_id, _x, _y, _z); }
	Point(const Matrix &mat);
	virtual ~Point() { }

	// Set attributes
	void SetId(const int _id) { id = _id; }
	void SetXYZ(const double _x, const double _y, const double _z) { x = _x; y = _y; z = _z; }
	void SetX(const double _x) { x = _x; }
	void SetY(const double _y) { y = _y; }
	void SetZ(const double _z) { z = _z; }
	void Set(const int _id, const double _x, const double _y, const double _z) { SetId(_id); SetXYZ(_x, _y, _z); }
	void Set(const int index, const double val)
	{
		switch(index)
		{
		case 0:	x = val; break;
		case 1:	y = val; break;
		case 2:	z = val; break;
		}
	}

	// Get attributes
	int Id() const { return id; };
	double X() const { return x; }
	double Y() const { return y; }
	double Z() const { return z; }
	void GetXYZ(double &_x, double &_y, double &_z) const { _x = x; _y = y; _z = z; }
	double Get(const int index) const
	{
		switch(index)
		{
		case 0:	return x;
		case 1: return y;
		case 2: return z;
		default: return 0;
		}
	}

	// Operations
	double Norm() const { return sqrt(x*x+y*y+z*z); }
	double Distance(const Point &pt) const { return sqrt((pt.x-x)*(pt.x-x)+(pt.y-y)*(pt.y-y)+(pt.z-z)*(pt.z-z)); }

	// Comparisons
	bool operator <(const Point &pt) const { if(id<pt.Id()) return true; else return false; }
	Point& operator =(const Point &pt) {	if(this!=&pt) Set(pt.id, pt.x, pt.y, pt.z);	return *this; }
	Point operator +(const Point &pt) const { return Point(x+pt.x, y+pt.y, z+pt.z); }
	Point operator -(const Point &pt) const { return Point(x-pt.x, y-pt.y, z-pt.z); }
	Point& operator +=(const Point &pt) { x+=pt.x; y+=pt.y; z+=pt.z; return *this; }
	Point& operator -=(const Point &pt) { x-=pt.x; y-=pt.y; z-=pt.z; return *this; }
	Point operator *(const double val) const { return Point(x*val, y*val, z*val); }
	Point& operator *=(const double val) { x*=val; y*=val; z*=val; return *this; }
	Point operator /(const double val) const { return Point(x/val, y/val, z/val); }
	Point& operator /=(const double val) { x/=val; y/=val; z/=val; return *this; }

	ostream& Print(const char *name, ostream &out=cout) const;
	friend ostream& operator<<(ostream& os, const Point& pt);
};

#endif // !defined(AFX_POINT_H__D4C0D00E_6FC6_11D5_9446_00105A29FE9E__INCLUDED_)
