/*
Name: MatrixOperations.hpp
Author: Jay Gambetta

Dependences: Matrix.hpp
Brief Discription: A namespace contianing all my matrix operations 
Limitations: None

Version History
	v0: June 14, 2007
	v1:	June 27, 2008. Changed to support my new storage for hermitian matrices
	v2: October 13th, 20008 - added the trace support
	v3: January 3rd, 2009 - changed it to support my new matrix template -- removed hermitian and symmetric matrix
	v4: May 20, 2009 - Fixed some bugs

*/
#ifndef MatrixOperations_h
#define MatrixOperations_h
#include "Matrix.hpp"


namespace MOs { //All my matrix operators
	
	// Type of matrices
	template <class T>
	bool IsSquare (const matrix<T>& A) { return (A.GetRows() == A.GetColumns() ); } 
	template <class T>
	bool IsDiagonal (const matrix<T>& A);
	template <class T>
	bool IsScalar (const matrix<T>& A);
	template <class T>
	bool IsNull (const matrix<T>& A);
	template <class T>
	bool IsSymmetric (const matrix<T>& A);
	template <class T>
	bool IsSkewSymmetric (const matrix<T>& A);
	template <class T>
	bool IsUpperTriangular (const matrix<T>& A);
	template <class T>
	bool IsLowerTriangular (const matrix<T>& A);
	bool IsHermitian (const matrix<std::complex<double> >& A);
	
	//Common Matricies
	template <class T> 
	void Null(matrix<T>& A);
	template <class T>
	void Identity(matrix<T>& I);
	template <class T> 
	void Destroy(matrix<T>& a);
	template <class T> 
	void GenPauliX(matrix<T>& a, size_t level_i, size_t level_j);
	template <class T>
	void GenPauliMinus(matrix<T>& a, size_t level_i, size_t level_j);	
        template <class T>
	void GenPauliZ(matrix<T>& a, size_t level_i, size_t level_j);
        template <class T> 
	void numberSymetry(matrix<T>& n);
	template <class T> 
	void costheta(matrix<T>& ct);
	void Pauli(matrix<std::complex<double> >& sx, matrix<std::complex<double> >& sy, matrix<std::complex<double> >& sz);
	void GellMann(matrix<std::complex<double> >& l1, matrix<std::complex<double> >& l2, matrix<std::complex<double> >& l3, matrix<std::complex<double> >& l4, matrix<std::complex<double> >& l5, matrix<std::complex<double> >& l6, matrix<std::complex<double> >& l7, matrix<std::complex<double> >& l8);
	
	template <class T> 
	void RandomHermitian(matrix<T>& A);
	
	
	//operations on Matrics
	template <class T> 
	matrix<T> Transpose(const matrix<T>& A);
	matrix<std::complex<double> > Dagger(const matrix<std::complex<double> >& A);
	matrix<std::complex<double> > Conjugate(const matrix<std::complex<double> >& A);
	template <class T>
	T Trace(const matrix<T>& A);
	template <class T>
	matrix<T> TraceOutA(const matrix<T>&, size_t);
	template <class T>
	matrix<T> TraceOutB(const matrix<T>&, size_t);
	template <class T>
	matrix<T> TensorProduct(const matrix<T>& A, const matrix<T>& B);
	inline void diagmult(matrix<std::complex<double> >& dest, std::complex<double> *diagelem,  const matrix<std::complex<double> >& source);
	inline void multdiag(matrix<std::complex<double> >& dest, const matrix<std::complex<double> >& source, std::complex<double>* diagelem);	

}
//------------Member funtions------------------
// inline bool MOs::IsSingular () { 
// 	// Determine if the matrix is singular
//    if (_m->Row != _m->Col)
//       return false;
//    return (Det() == T(0));
// }
template <class T>
inline bool MOs::IsDiagonal (const matrix<T>& A){
	// Determine if the matrix is diagonal
   	if (A.GetRows() != A.GetColumns())
      	return false;
   	for (size_t i=0; i < A.GetRows(); i++)
     	for (size_t j=0; j < A.GetColumns(); j++)
			if (i != j && A(i,j) != T(0))
	  			return false;
   return true;
}
template <class T>
inline bool MOs::IsScalar (const matrix<T>& A){
	// Determine if the matrix is scalar
   	if (!IsDiagonal(A))
     	return false;
   	T v = A(0,0);
   	for (size_t i=1; i < A.GetRows(); i++)
     	if (A(i,i) != v)
			return false;
	return true;
}
template <class T>
inline bool MOs::IsNull (const matrix<T>& A){
	// Determine if this is a null matrix
   	for (size_t i=0; i < A.GetRows(); i++)
     	for (size_t j=0; j < A.GetColumns(); j++)
	 		if (A(i,j) != T(0))
	    		return false;
   return true;
}
template <class T>
inline bool MOs::IsSymmetric (const matrix<T>& A){
	// Determine if the matrix is symmetric
   	if (A.GetRows() != A.GetColumns())
      	return false;
   	for (size_t i=0; i < A.GetRows(); i++)
     	for (size_t j=0; j < A.GetColumns(); j++)
	 		if (A(i,j) != A(j,i))
	    		return false;
   return true;
}
template <class T>
inline bool MOs::IsSkewSymmetric (const matrix<T>& A){
	// Determine if the matrix is skew-symmetric
   	if (A.GetRows() != A.GetColumns())
      	return false;
   	for (size_t i=0; i < A.GetRows(); i++)
     	for (size_t j=0; j < A.GetColumns(); j++)
	 		if (A(i,j) != -A(j,i))
	    		return false;
   return true;
}
template <class T>
inline bool MOs::IsUpperTriangular (const matrix<T>& A){
	// Determine if the matrix is upper triangular
	if (A.GetRows() != A.GetColumns())
      	return false;
	for (size_t i=1; i < A.GetRows(); i++)
     	for (size_t j=0; j < i; j++)
	 		if (A(i,j) != T(0))
	    		return false;
	return true;
}
template <class T>
inline bool MOs::IsLowerTriangular (const matrix<T>& A){
	// Determine if the matrix is lower triangular
	if (A.GetRows() != A.GetColumns())
      	return false;

   	for (size_t j=1; j < A.GetColumns(); j++)
      	for (size_t i=0; i < j; i++)
	 		if (A(i,j) != T(0))
	    		return false;

	return true;
}
inline bool MOs::IsHermitian (const matrix<std::complex<double> >& A){
	// Determine if the matrix is Hermitian
   	if (A.GetRows() != A.GetColumns())
      	return false;
   	for (size_t i=0; i < A.GetRows(); i++)
     	for (size_t j=0; j < A.GetColumns(); j++)
	 		if (A(i,j) != conj(A(j,i)))
	    		return false;
   return true;
}
template <class T> 
inline void MOs::Null(matrix<T>& A){
	// set zero to all elements of this matrix
   	for (size_t i=0; i < A.GetRows(); i++)
     	for (size_t j=0; j < A.GetColumns(); j++)
		 	A(i,j)= T(0);
    return;
}
template <class T> 
inline void MOs::Identity(matrix<T>& I){
	//The Idenity operator
	size_t dim = I.GetRows();
	if (dim != I.GetColumns())
			UFs::MyError("Routine MOs::Identity: Matrix is not square");
	for(size_t i=0; i<dim; i++)
		for(size_t j=0; j<dim; j++)
				I(i,j) = i == j ? T(1) : T(0);
	return;	
}
template <class T> 
inline void MOs::RandomHermitian(matrix<T>& A){
	size_t dim = A.GetRows();
	
	/* initialize random seed: do it earlier or use clock() perhaps*/
	//srand ( time(NULL) );
		
	for(size_t i = 0; i < dim; i++){
		A(i,i)=(rand()%100)/100.0;
		for (size_t j=0; j < i; j++){
			A(i,j)=T((rand()%100)/100.0,(rand()%100)/100.0);
			A(j,i) = conj(A(i,j));
		}
	}
	return;	
}

template <class T> 
inline void  MOs::GenPauliX(matrix<T>& a, size_t level_i, size_t level_j)
{
	Identity(a);
	a(level_i,level_i)=a(level_j,level_j)=T(0);
	a(level_i,level_j)=a(level_j,level_i)=T(1);
}
template <class T> 
inline void  MOs::GenPauliMinus(matrix<T>& a, size_t level_i, size_t level_j)
{
	Identity(a);
	
	a(level_i,level_j)=T(1);
	a(level_j,level_j)=a(level_i, level_i)=a(level_j,level_i)=T(0);
}





//creates a paul Z matrix
template <class T> 
inline void  MOs::GenPauliZ(matrix<T>& a, size_t level_i, size_t level_j)
{
	Identity(a);
	a(level_i,level_i)=T(1);
        a(level_i,level_j)=a(level_j,level_i)=T(0);
	a(level_j,level_j)=T(-1);
}
template <class T> 
inline void MOs::Destroy(matrix<T>& a){
	// The Annilation operator
	size_t dim = a.GetRows();
	if (dim != a.GetColumns())
			UFs::MyError("Routine MOs::Destroy: Matrix is not square");
	for(size_t i=0; i<dim; i++)
		for(size_t j=0; j<dim; j++)
			a(i,j) = j == i+1 ? std::sqrt(T(j)) : T(0);	
	return;	
}

template <class T> 
inline void MOs::numberSymetry(matrix<T>& n){
	// The symetry number operator - (dim-1)/2 to (dim-1)/2
	size_t dim = n.GetRows();
	if (dim != n.GetColumns())
			UFs::MyError("Routine MOs::Destroy: Matrix is not square");
	for(size_t i=0; i<dim; i++)
		for(size_t j=0; j<dim; j++)
			n(i,j) = j == i ? T(double(j)-(double(dim-1))/2) : T(0);	
	return;	
}
template <class T> 
inline void MOs::costheta(matrix<T>& ct){
	// The cos(phi) operator
	size_t dim = ct.GetRows();
	if (dim != ct.GetColumns())
			UFs::MyError("Routine MOs::Destroy: Matrix is not square");
	for(size_t i=0; i<dim; i++)
		for(size_t j=0; j<dim; j++)
			ct(i,j) = j == i+1 ? T(0.5) : j == i-1 ? T(0.5) : T(0);	
	return;
}
inline void MOs::Pauli(matrix<std::complex<double> >& sx, matrix<std::complex<double> >& sy, matrix<std::complex<double> >& sz){
	// Pauli x
	sx(0,0)=std::complex<double>(0);
	sx(1,0)=std::complex<double>(1);
	sx(0,1)=std::complex<double>(1);
	sx(1,1)=std::complex<double>(0);
	
	// Pauli y
	sy(0,0)=std::complex<double>(0);
	sy(1,0)=std::complex<double>(0,1);
	sy(0,1)=std::complex<double>(0,-1);
	sy(1,1)=std::complex<double>(0);
	
	// Pauli z
	sz(0,0)=std::complex<double>(1);
	sz(1,0)=std::complex<double>(0);
	sz(0,1)=std::complex<double>(0);
	sz(1,1)=std::complex<double>(-1);
	
	return;
}	
inline void MOs::GellMann(matrix<std::complex<double> >& l1, matrix<std::complex<double> >& l2, matrix<std::complex<double> >& l3, matrix<std::complex<double> >& l4, matrix<std::complex<double> >& l5, matrix<std::complex<double> >& l6, matrix<std::complex<double> >& l7, matrix<std::complex<double> >& l8){
	MOs::Null(l1);
	MOs::Null(l2);
	MOs::Null(l3);
	MOs::Null(l4);
	MOs::Null(l5);
	MOs::Null(l6);
	MOs::Null(l7);
	MOs::Null(l8);

	l1(1,0)=std::complex<double>(1);
	l1(0,1)=std::complex<double>(1);
	
	l2(1,0)=std::complex<double>(0,1);
	l2(0,1)=std::complex<double>(0,-1);
	
	l3(0,0)=std::complex<double>(1);
	l3(1,1)=std::complex<double>(-1);
	
	l4(2,0)=std::complex<double>(1);
	l4(0,2)=std::complex<double>(1);
	
	l5(2,0)=std::complex<double>(0,1);
	l5(0,2)=std::complex<double>(0,-1);
	
	l6(2,1)=std::complex<double>(1);
	l6(1,2)=std::complex<double>(1);
	
	l7(2,1)=std::complex<double>(0,1);
	l7(1,2)=std::complex<double>(0,-1);
	
	l8(0,0)=std::complex<double>(1.0/std::sqrt(3.0));
	l8(1,1)=std::complex<double>(1.0/std::sqrt(3.0));
	l8(2,2)=std::complex<double>(-2.0/std::sqrt(3.0));

}
template <class T>
inline matrix<T> MOs::Transpose(const matrix<T>& A){
	//Transposes a Matrix
	size_t rows = A.GetRows(), cols= A.GetColumns();
	matrix<T> temp(cols,rows);
	for (size_t i=0; i < rows; i++){
		for(size_t j=0; j < cols; j++){
			temp(j,i) = A(i,j);
		}
	}
	return temp;
}
inline matrix<std::complex<double> > MOs::Dagger(const matrix<std::complex<double> >& A){
	//Take the Hermitian conjugate of a complex matrix
	size_t cols = A.GetColumns(), rows= A.GetRows();
	matrix<std::complex<double> > temp(cols,rows);
	for (size_t i=0; i < rows; i++){
		for(size_t j=0; j < cols; j++){
			temp(j,i) = conj(A(i,j));
		}
	}
	return temp;
}

inline matrix<std::complex<double> > MOs::Conjugate(const matrix<std::complex<double> >& A){
	//Take the Hermitian conjugate of a complex matrix
	size_t rows = A.GetRows(), cols= A.GetColumns();
	matrix<std::complex<double> > temp(rows,cols);
	for (size_t i=0; i < rows; i++){
		for(size_t j=0; j < cols; j++){
			temp(i,j) = conj(A(i,j));
		}
	}
	return temp;
}
template <class T>
inline T MOs::Trace(const matrix<T>& A){
	//Finds the trace of a matrix
	size_t rows = A.GetRows(), cols= A.GetColumns();
	if(rows !=  cols)
		UFs::MyError("Routine MOs::Trace: Matrix is not square");
	T temp=0.0;
	for (size_t i=0; i < rows; i++){
		temp = temp + A(i,i);
	}
	return temp;
}

template <class T>
inline matrix<T> MOs::TraceOutA(const matrix<T>& rho, size_t dimA){
	// Traces out first system (dimension dimA) of composite Hilbert space
	size_t rows = rho.GetRows(), cols= rho.GetColumns();
	if(rows !=  cols)
		UFs::MyError("Routine MOs::TraceOutA: Matrix is not square");
	if( rows % dimA != 0)
		UFs::MyError("Routine MOs::TraceOutA: dim(rho)/dim(system B) is not an integer");
	size_t dimB= rows/dimA;
	matrix<T> rhoB(dimB, dimB);	
	T temp=0.0;
	for (size_t i=0; i < dimB; i++){
		for (size_t j=0; j<dimB; j++){
			for (size_t k=0; k < dimA; k++){
				temp = temp + rho(i+dimB*k, j+dimB*k);
			}
			rhoB(i,j) = temp;
			temp = 0.0;
		}
	}
	return rhoB;
}
template <class T>
inline matrix<T> MOs::TraceOutB(const matrix<T>& rho, size_t dimB){
	// Traces out second system (dimension dimB) of composite Hilbert space
	size_t rows = rho.GetRows(), cols= rho.GetColumns();
	if(rows !=  cols)
		UFs::MyError("Routine MOs::TraceOutB: Matrix is not square");
	if( rows % dimB != 0)
		UFs::MyError("Routine MOs::TraceOutB: dim(rho)/dim(system B) is not an integer");
	size_t dimA= rows/dimB;
	matrix<T> rhoA(dimA, dimA);	
	T temp=0.0;
	for (size_t i=0; i < dimA; i++){
		size_t offsetX = i*dimB;
		for (size_t j=0; j<dimA; j++){
			size_t offsetY = j*dimB;
			for (size_t k=0; k < dimB; k++){
				temp = temp + rho(offsetX+k, offsetY+k);
			}
			rhoA(i,j) = temp;
			temp = 0.0;
		}
	}
	return rhoA;
}
template <class T>
inline matrix<T> MOs::TensorProduct(const matrix<T>& A, const matrix<T>& B){
	//Works out the TensorProduct of two matricies A tensor B
	//Note that if A is i x j and B is p x q then A \otimes B is an ip x jq rmatrix
	size_t rows1 = A.GetRows(), rows2 = B.GetRows(), cols1 = A.GetColumns(), cols2 = B.GetColumns();
	size_t rows_new = rows1*rows2, cols_new = cols1*cols2, n, m;
	matrix<T> temp(rows_new,cols_new);
	//a11 B, a12 B ... a1j B
	//ai1 B, ai2 B ... aij B
	for(size_t i=0; i<rows1; i++) {
		for(size_t j=0; j<cols1; j++){ 
			for(size_t p=0; p<rows2; p++){
				for(size_t q=0; q<cols2; q++){
					n = i*rows2+p;    
					m = j*cols2+q; //  0 (0 + 1)  + 1*dimb=2 + (0 + 1 )  (j*dimb+q) 		
					temp(n,m) = A(i,j)*B(p,q);
				}
			}
		}
	}
	return temp;
}

inline void  MOs::diagmult(matrix<std::complex<double> >& dest, std::complex<double>*diagelem,  const matrix<std::complex<double> >& source)
{
	size_t rows= dest.GetRows(), cols = dest.GetColumns();
	for (size_t i=0; i < rows; i++){
		for(size_t j=0; j < cols; j++){
			dest(i,j) = diagelem[i]*source(i,j);
		}
	}
}

inline void  MOs::multdiag(matrix<std::complex<double> >& dest,  const matrix<std::complex<double> >& source, std::complex<double> *diagelem)
{
	size_t rows= dest.GetRows(), cols = dest.GetColumns();
	for (size_t i=0; i < rows; i++){
		for(size_t j=0; j < cols; j++){
			dest(i,j) = diagelem[j]*source(i,j);
		}
	}
}


#endif /* MatrixOperations_h */
