// Matrix.h: interface for the Matrix class.
//
//////////////////////////////////////////////////////////////////////

/* Note: This version is with additional functions listed below:
 *			Clear() : To reset the content of the matrix
 *			EuclideanNorm() : To compute the euclidean norm of matrix
 */

#if !defined(AFX_MATRIX_H__F187AE80_753B_11D5_9447_00105A29FE9E__INCLUDED_)
#define AFX_MATRIX_H__F187AE80_753B_11D5_9447_00105A29FE9E__INCLUDED_

#if _MSC_VER > 1000
#pragma once
#endif // _MSC_VER > 1000

#include "Array.h"

#include <vector>
#include <istream>
#include <ostream>

using namespace std;

class Point;

class Matrix : public Array<double>
{
public:
	// Constructors and Destructors
	Matrix() { 	SetDim(0,0); }

	Matrix(const int _row, const int _col) { SetDim(_row, _col); }
	Matrix(const int _row, const int _col, double val) { SetDim(_row, _col); Set(val); }
	Matrix(const Matrix &mat) { *this = mat; }
	Matrix(const Point &pt);
	virtual ~Matrix() { };

	// Matrix Operators
	Matrix& Copy(const Matrix &A);
	Matrix& operator= (const Matrix &A);
	Matrix operator+ (const Matrix &B) const;
	Matrix& operator+= (const Matrix &B);
	Matrix operator- (const Matrix &B) const;
	Matrix operator- (const double &a) const;
	Matrix& operator-= (const Matrix &B);
	Matrix operator* (const Matrix &B) const;
	Matrix operator* (const double val) const;
	Matrix& operator*= (const double val);
	Matrix Sqrt() const;
	Matrix Mult(const Matrix &A) const;
	Matrix& Mult(const Matrix &A, Matrix &T) const;

	// Format Change
	Matrix Trans() const;
	Matrix Floor() const;
	Matrix Round() const;
	Matrix Diag(const bool option=true) const;
	Matrix Subset(const int si, const int ei, const bool option=true) const;
	Matrix& Subset(Matrix &A, const int si, const int ei, const bool option=true) const;
	Matrix Subset(const int rs, const int re, const int cs, const int ce) const;
	Matrix& Subset(Matrix &A, const int rs, const int re, const int cs, const int ce) const;
	Matrix DiagSubset();
	Matrix& DiagSubset(Matrix &A);

	// Special Matrix
	Matrix& Eye(const int dim);
	Matrix& RotateMat(const double Om, const double Ph, const double Kp);

	// Element Functions
	void swap(double &x, double &y) const { double tmp=y; y=x; x=tmp; }

	// Elementary Vector Functions
	Matrix Sort() const;
	double Median() const;
	double Norm(const int index, const bool option = true) const;
	double EuclideanNorm() const;
	double SumSquare(const int index, const bool option = true) const;
	double SumSquare(const Matrix &A, const int index, const bool option = true) const;
	
	// Matrix Functions

	double Det() const;
	bool Inv2(Matrix &T) const;
	double InnerProduct(const Matrix &A) const;
	Matrix DiagInv() const;
	Matrix& DiagInv(Matrix &A) const;
	Matrix Inv() const;
	Matrix& Inv(Matrix &A) const;
	Matrix& Inv(Matrix &A, bool &flag) const;
	bool InvCheck(Matrix &A) const;
	bool GaussJordan(Matrix &A, Matrix &b) const;

	// Debug the matrix.
	void DebugMatrix(Matrix &A)
	{
		for (int r = 0; r < row; r++)
			for (int c = 0; c < col; c++)
				printf("[%d,%d] = %.10lf\n", r, c, A.Get(r,c));
	}

	// Reset the content of the matrix.
	void Clear();

	// Input and Ouput through a stream
	void ImportBinary(FILE* fin);
	void ExportBinary(FILE* fout) const;
	friend istream& operator>>(istream& is, Matrix& A);
	friend ostream& operator<<(ostream& os, const Matrix& A);

};

#endif // !defined(AFX_MATRIX_H__F187AE80_753B_11D5_9447_00105A29FE9E__INCLUDED_)
