// Matrix.cpp: implementation of the Matrix class.
//
//////////////////////////////////////////////////////////////////////

/* Note: This version is with additional functions listed below:
 *			Clear() : To reset the content of the matrix
 */

#include "AT_Matrix.h"
#include "AT_Point.h"

#include <istream>
#include <ostream>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Matrix::Matrix(const Point &pt)
{
	SetDim(3,1);
	(*this)[0] = pt.X();
	(*this)[1] = pt.Y();
	(*this)[2] = pt.Z();
}

Matrix& Matrix::Copy(const Matrix &A)
{
	if (this == &A) return *this;
	SetDim(A.row,A.col);
	for(register int n=0; n<size(); n++) (*this)[n] = A[n];
	return *this;
}

Matrix& Matrix::operator= (const Matrix &A)
{
	return Copy(A);
}

Matrix Matrix::operator+ (const Matrix &B) const
{
	Matrix T(row,col);
	for(register int n=0; n<size(); n++) T[n] = (*this)[n] + B[n];
	return T;
}

Matrix& Matrix::operator+= (const Matrix &B)
{
	for(register int n=0; n<size(); n++) (*this)[n] -= B[n];
	return *this;
}

Matrix Matrix::operator- (const Matrix &B) const
{
	Matrix T(row,col);
	for(register int n=0; n<size(); n++) T[n] = (*this)[n] - B[n];
	return T;
}

Matrix Matrix::operator- (const double &a) const
{
	Matrix T(row,col);
	for(register int n=0; n<size(); n++) T[n] = (*this)[n] - a;
	return T;
}

Matrix& Matrix::operator-= (const Matrix &B)
{
	for(register int n=0; n<size(); n++) (*this)[n] -= B[n];
	return *this;
}

Matrix Matrix::operator* (const Matrix &B) const
{
	Matrix T(row,B.col);
	int i,j,k;
	double sum;
	for(i=0; i<T.row; i++) for(j=0; j<T.col; j++)
	{
		sum=0.0;
		for(k=0; k<col; k++) sum += Get(i,k) * B.Get(k,j);
		T.Set(i,j,sum);
	}
	return T;
}

Matrix Matrix::operator* (const double val) const
{
	Matrix T(row,col);
	for(register int n=0; n<size(); n++) T[n] = (*this)[n] * val;
	return T;
}

Matrix& Matrix::operator*= (const double val)
{
	for(register int n=0; n<size(); n++) (*this)[n] *= val;
	return *this;
}

Matrix Matrix::Sqrt() const
{
	Matrix R(row,col);
	for(register int n=0; n<size(); n++) R[n] = sqrt((*this)[n]);
	return R;
}

Matrix Matrix::Mult(const Matrix &A) const
{
	return Mult(A, Matrix());
}

Matrix& Matrix::Mult(const Matrix &A, Matrix &T) const
{
	T.SetDim(row,col);
	for(register int n=0; n<size(); n++) T[n] = (*this)[n] * A[n];

	return T;
}

Matrix Matrix::Trans() const
{
	Matrix T(col,row);
	int i,j;
	for(i=0; i<T.row; i++) for(j=0; j<T.col; j++)
		T.Set(i,j,Get(j,i));
	return T;
}

Matrix Matrix::Floor() const
{
	Matrix T(row,col);
	int i,j;
	for(i=0; i<T.row; i++) for(j=0; j<T.col; j++)
		T.Set(i,j, int(Get(i,j)));
	return T;
}

Matrix Matrix::Round() const
{
	Matrix T(row,col);
	int i,j;
	for(i=0; i<T.row; i++) for(j=0; j<T.col; j++)
		T.Set(i,j, int(Get(i,j)+0.5));
	return T;
}

Matrix Matrix::Diag(const bool option) const
{
	Matrix T;
	if(option)
	{	
		T.Set(row,col,0.0);
		for (int m=0; m<row; m++) T.Set(m,m,Get(m,m));
	}
	else
	{
		T.SetDim(row,1);
		for (int m=0; m<row; m++) T.Set(m,1,Get(m,m));
	}
	return T;
}

Matrix Matrix::Subset(const int si, const int ei, const bool option) const
{
	if(option) return Subset(Matrix(), si, ei, 0, col-1);
	else return Subset(Matrix(), 0, row-1, si, ei);
}

Matrix& Matrix::Subset(Matrix &A, const int si, const int ei, const bool option) const
{
	if(option) return Subset(A, si, ei, 0, col-1);
	else return Subset(A, 0, row-1, si, ei);
}

Matrix Matrix::Subset(const int rs, const int re, const int cs, const int ce) const
{
	return Subset(Matrix(), rs, re, cs, ce);
}

Matrix& Matrix::Subset(Matrix &A, const int rs, const int re, const int cs, const int ce) const
{
	int nrow = (re-rs)+1;
	int ncol = (ce-cs)+1;
	A.SetDim(nrow,ncol);

	int r,c;
	for(r=0; r<nrow; r++) for(c=0; c<ncol; c++) A[Idx(r,c)] = Get(r+rs,c+cs);
	return A;
}

Matrix& Matrix::DiagSubset(Matrix &A)
{
	for (int m=0; m<row; m++) A[Idx(m,m)] = Get(m,m);
	return A;
}

Matrix& Matrix::Eye(const int dim)
{
	SetDim(dim,dim);
	Set(0);
	for(register int n=0; n<dim; n++) Set(n,n,1);
	return *this;
}


Matrix& Matrix::RotateMat(const double Om, const double Ph, const double Kp)
{
	SetDim(3,3);
	double cOm = cos(Om);
	double sOm = sin(Om);
	double cPh = cos(Ph);
	double sPh = sin(Ph);
	double cKp = cos(Kp);
	double sKp = sin(Kp);
	Set(0,0,  cPh * cKp);
	Set(0,1, -cPh * sKp);
	Set(0,2,  sPh);
	Set(1,0,  cOm * sKp + sOm * sPh * cKp);
	Set(1,1,  cOm * cKp - sOm * sPh * sKp);
	Set(1,2, -sOm * cPh);
	Set(2,0,  sOm * sKp - cOm * sPh * cKp);
	Set(2,1,  sOm * cKp + cOm * sPh * sKp);
	Set(2,2,  cOm * cPh);
	return *this;
}

Matrix Matrix::Sort() const
{
	Matrix M(*this);
	sort(M.begin(), M.end());
	return M;
}

double Matrix::Median() const
{
	Matrix M = Sort();
	double medVal;
	if((size()%2)==1)
	{
		int index = (size()-1)/2;
		medVal = M[index];
	}
	else
	{
		int index = size()/2;
		medVal = (M[index-1] + M[index])/2.0;
	}
	return medVal;
}

double Matrix::Norm(const int index, const bool option) const
{
	double sum = 0.;
	if(option)	// row
	{
		for(register int k=0; k<col; k++)
			sum += Get(index,k) * Get(index,k);
	}
	else		// column
	{
		for(register int k=0; k<row; k++)
			sum += Get(k,index) * Get(k,index);
	}
	return sqrt(sum);
}

double Matrix::SumSquare(const int index, const bool option) const
{
	double sum = 0.;
	if(option)	// row
	{
		for(register int k=0; k<col; k++)
			sum += Get(index,k) * Get(index,k);
	}
	else		// column
	{
		for(register int k=0; k<row; k++)
			sum += Get(k,index) * Get(k,index);
	}
	return sum;
}

double Matrix::SumSquare(const Matrix &A, const int index, const bool option) const
{
	double sum = 0.;
	if(option)	// row
	{
		for(register int k=0; k<col; k++)
			sum += Get(index,k) * A.Get(index,k);
	}
	else		// column
	{
		for(register int k=0; k<row; k++)
			sum += Get(k,index) * A.Get(k,index);
	}
	return sum;
}


double Matrix::Det() const
{
	Matrix T(*this);

	double tmp, pivot;
	int i,j,l;

	int num_change = 0;
	bool flag = false;
	for(i=0; i<row-1; i++)	  // Gauss Elli.
	{
		if(T.Get(i,i)==0.)
		{
			j = i+1;
			while((j<row) && (T.Get(j,i)==0.) ) j++; //i:col, j:row		
			if(j>=row)	//one row or one column is whole zero
			{ 
				flag = true;
				break;	
			}
			for(l=i; l<col; l++) swap(T[Idx(i,l)], T[Idx(j,l)]);
			num_change++;
		}
	    pivot = T.Get(i,i);
		for(j=i+1; j<row ;j++ )
		{
			if(T.Get(j,i)!=0.)
			{
				tmp = T.Get(j,i)/pivot;
				for(l=i; l<col; l++) T[Idx(j,l)]-=T.Get(i,l)*tmp;
			}
		}
	}
	
	if(flag) return(0.);	//sigular
	
	// determinant: multiply of diagonal elements
	double det=T[0];	//first elemets
	for(i=1 ;i<row ;i++) det *= T.Get(i,i);			
	
	//odd change of row elements
	if((num_change%2) != 0)	det *= -1;
	return det;
}



bool Matrix::Inv2(Matrix &T) const
{
	double det = Get(0,0)*Get(1,1) - Get(0,1)*Get(1,0);
	if(det==0.) return false;
	else
	{
		T.SetDim(2,2);
		T.Set(0,0,  Get(1,1)/det);
		T.Set(1,1,  Get(0,0)/det);
		T.Set(0,1, -Get(0,1)/det);
		T.Set(1,0, -Get(1,0)/det);
		return true;
	}
}


Matrix Matrix::Inv() const
{
	return Inv(Matrix());
}

Matrix& Matrix::Inv(Matrix &A) const
{
	Matrix I;
	I.Eye(row);
	GaussJordan(A,I);
	return A;
}

Matrix& Matrix::Inv(Matrix &A, bool &flag) const
{
	Matrix I;
	I.Eye(row);
	flag = GaussJordan(A,I);
	return A;
}

bool Matrix::InvCheck(Matrix &A) const
{
	Matrix I;
	I.Eye(row);
	return GaussJordan(A,I);
}

double Matrix::InnerProduct(const Matrix &A) const
{
	double val = 0.0;
	for(register int n=0; n<size(); n++) val += (*this)[n] * A[n];
	return val;
}

Matrix Matrix::DiagInv() const
{
	return DiagInv(Matrix());
}

Matrix& Matrix::DiagInv(Matrix &A) const
{
	A.SetDim(row,col);
	for (int m=0; m<row; m++) A[Idx(m,m)] = 1/Get(m,m);
	return A;
}


// Porting from Numerical Recipe in C, p39-40
bool Matrix::GaussJordan(Matrix &A, Matrix &b) const
/*
Linear equation solution by Gauss-Jordan elimination, equation (2.1.1) above. a[1..n][1..n]
is the input matrix. b[1..n][1..m] is input containing the m right-hand side vectors. On
output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution
vectors.
*/
{
	A = *this;
	int n = A.NumRow();
	int m = b.NumCol();

	vector<int> indxc(n), indxr(n), ipiv(n);
	// The integer arrays ipiv, indxr, and indxc are used for bookkeeping on the pivoting.

	int i,icol,irow,j,k,l,ll;
	double big,dum,pivinv;
	
	for (j=0;j<n;j++) ipiv[j]=0;
	for (i=0;i<n;i++) // This is the main loop over the columns to be reduced. 
	{ 
		big=0.0;
		for (j=0;j<n;j++) if (ipiv[j] != 1) for (k=0;k<n;k++)
		// This is the outer loop of the search for a pivot element.
		{
			if (ipiv[k] == 0)
			{
				if (fabs(A.Get(j,k)) >= big)
				{
					big=fabs(A.Get(j,k));
					irow=j;
					icol=k;
				}
			}
			else if (ipiv[k] > 1) return false;
		}
		++(ipiv[icol]);
		
		/*
		We now have the pivot element, so we interchange rows, if needed, to put the pivot
		element on the diagonal. The columns are not physically interchanged, only relabeled:
		indxc[i], the column of the ith pivot element, is the ith column that is reduced, while
		indxr[i] is the row in which that pivot element was originally located. If indxr[i]
		6 = indxc[i] there is an implied column interchange. With this form of bookkeeping, the
		solution b's will end up in the correct order, and the inverse matrix will be scrambled
		by columns.
		*/
		
		if (irow != icol)
		{
			for (l=0;l<n;l++) swap(A[Idx(irow,l)],A[Idx(icol,l)]);
			for (l=0;l<m;l++) swap(b[Idx(irow,l)],b[Idx(icol,l)]);
		}
		indxr[i]=irow;
		
		// We are now ready to divide the pivot row by the pivot element, located at irow and icol.
		
		indxc[i]=icol;
		if (A.Get(icol,icol) == 0.0) return false;
		pivinv=1.0/A.Get(icol,icol);
		A[Idx(icol,icol)]=1.0;
		for (l=0;l<n;l++) A[Idx(icol,l)] *= pivinv;
		for (l=0;l<m;l++) b[Idx(icol,l)] *= pivinv;
		// Next, we reduce the rows... except for the pivot one, of course.
		for (ll=0;ll<n;ll++) if (ll != icol)
		{
			dum=A.Get(ll,icol);
			A[Idx(ll,icol)]=0.0;
			for (l=0;l<n;l++) A[Idx(ll,l)] -= A.Get(icol,l)*dum;
			for (l=0;l<m;l++) b[Idx(ll,l)] -= b.Get(icol,l)*dum;
		}
	}
	
	/*
	This is the end of the main loop over columns of the reduction. It only remains to unscram-
	ble the solution in view of the column interchanges. We do this by interchanging pairs of
	columns in the reverse order that the permutation was built up.
	*/
	for (l=n-1;l>=0;l--) if (indxr[l] != indxc[l]) for (k=0;k<n;k++)
		swap(A[Idx(k,indxr[l])],A[Idx(k,indxc[l])]);
	return true;
}

void Matrix::ExportBinary(FILE* fout) const
{
	fwrite(&((*this)[0]), sizeof(double), row*col, fout);
}

void Matrix::ImportBinary(FILE* fin)
{
	fread(&((*this)[0]), sizeof(double), row*col, fin);
}

istream& operator>>(istream& is, Matrix& A)
{
	int m,n;
	for(m=0; m<A.NumRow(); m++) for(n=0; n<A.NumCol(); n++) is >> A[A.Idx(m,n)];
	return is;
}

ostream& operator<<(ostream& os, const Matrix& A)
{
	int m,n;
	for(m=0; m<A.NumRow(); m++)
	{
		for(n=0; n<A.NumCol(); n++)
		{
			os << A.Get(m,n);
			if (n != A.NumCol() -1 )
			{
				os << "\t";
			}
		}
		if (m != A.NumRow() -1 )
		{
			os << endl;
		}
	}
	return os;
}

/*
 * Reset the content of the matrix.
 */
void Matrix::Clear()
{
	int r,c;
	
	for (r = 0; r < row; r++)
		for (c = 0; c < col; c++)
			Set(r, c, 0.0);
}


