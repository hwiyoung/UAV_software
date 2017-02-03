// Point.cpp: implementation of the Point class.
//
//////////////////////////////////////////////////////////////////////

#include "Point.h"
#include "Matrix.h"
#include <ostream>
#include <iomanip>

Point::Point(const Matrix &mat)
{
	SetXYZ(mat[0], mat[1], mat[2]);
}

ostream& operator<<(ostream& os, const Point& pt)
{
	os << pt.X() << "\t" << pt.Y() << "\t" << pt.Z();
	return os;
}

