/*
Name: Matrix.hpp
Author: Andreas Brunner

Dependences: Matrix.hpp, UtilityFunctions.hpp, BLAS

Brief Discription: This is the sparse Matrix extension. As in "matrix", this works only
with real/complex matrices and stores the entries in 	column-major-order.
In column-major storage the columns are stored one after the other. The linear offset p 
from the beginning of the array to any given element A(i,j) can then be computed as:
		p = j*Nrows+i
	where Nrows is the number of rows in the matrix. If i wanted to know the i and j 
	associtated with a given p then i woul duse 
		i=p %Nrows
		j= floor(p/Nrows)

Limitations:	Only support the obvious sparse representation saving non-zero elements' positions.
				
Version History
	v0: March 05th, 2009. Start
*/

#ifndef spMatrix_h
#define spMatrix_h
#include "Matrix.hpp"
#define SPTHRESHOLD 1e-10	// entries smaller than this will be declared zero

template <class U>
class spMatrix : public matrix<U>
{	

	// Multiplication (does not catch an error of a S1 = real and S2 being complex)
	template <class S1, class S2>
	friend spMatrix<S1> operator*(const S2& beta, const spMatrix<S1>& A);//multiplication by a scalar beta*A
	template <class S1, class S2>
	friend spMatrix<S1> operator*(const spMatrix<S1>& A, const S2& beta);//multiplication by a scalar A*beta
	template <class V>
	friend matrix<V> operator*(const spMatrix<V>& A, const matrix<V>& B);//sparse*complex matrix multplication A*B
	template <class V>
	friend matrix<V> operator*(const matrix<V>& A, const spMatrix<V>& B);//complex*sparse
/*	template <class V>	// TODO: when adding this I get a "ambiguous definition" error (conflict with spMat=spMat*spMat)
	friend matrix<V> operator*(const spMatrix<V>& A, const spMatrix<V>& B);//dense = sparse*sparse 
*/	template <class V>
	friend spMatrix<V> operator*(const spMatrix<V>& A, const spMatrix<V>& B);//sparse = sparse*sparse 

public:
	//Construtors
	spMatrix();
	spMatrix(size_t rows, size_t cols);
	spMatrix(size_t rows, size_t cols, size_t elements);
	spMatrix(const spMatrix<U>&);
	spMatrix(const matrix<U>&);

	//No need for spMatrix.initialize(), matrix.initialize is sufficient!

	//Destructor TODO: I get terrible memory/stack errors when implementing a destructor!
//	~spMatrix(); 



	//Assignment operator
/*	matrix<T>& operator=(const matrix<T>& m);
	template <class S>
	matrix<T>& operator=(const matrix<S>& m); //Still would like to have real assigend by complex -- take real part
*/	

/*	
	//overloading functions. 
	matrix<T> operator+(const matrix<T>& A);
	matrix<T> operator-(const matrix<T>& A);
	matrix<T> operator+(const matrix<T>& A) const;
	matrix<T> operator-(const matrix<T>& A) const;
	matrix<T>& operator+=(const matrix<T>& A);
	matrix<T>& operator-=(const matrix<T>& A);
*/
	void Update();
	size_t GetElements() const;
	size_t* GetPositions() const;
	double GiveRatio(); //		return sqrt( elements_/size_ )
protected:
	size_t* positions_;		// which positions are non-zero
	size_t elements_;			// how many of those
};


//------------ sparse Member funtions------------------

template <class U>
inline spMatrix<U>::spMatrix(){
	this->rows_=0;
	this->cols_=0;
	this->size_=0;
	this->LD_=0;
	this->outputstyle_ = Column;
	this->mat_ = new U [this->size_];
	elements_=0;
}
template <class U>
inline spMatrix<U>::spMatrix(size_t a, size_t b){
	this->rows_=a;
	this->cols_=b;
	this->size_=a*b;
	this->LD_=a;
	this->outputstyle_ = Column;
	this->mat_ = new U [this->size_];// dont care about memory, most elements are assigned
	elements_=0;
}
template <class U>
inline spMatrix<U>::spMatrix(size_t a, size_t b, size_t elements){
	this->rows_=a;
	this->cols_=b;
	this->size_=a*b;
	this->LD_=a;
	this->outputstyle_ = Column;
	this->mat_ = new U [this->size_];// dont care about memory, most elements are assigned
	elements_=elements;
	positions_ = new size_t [elements_];
}
template <class U>
inline spMatrix<U>::spMatrix(const matrix<U>& M) : matrix<U>::matrix(M) {
/*	this->rows_=M.GetRows();	// TODO: shouldnt I be able to access M.cols_ directly?
	this->cols_=M.GetColumns();
	this->size_=M.size();
	this->LD_=M.GetLD();
	this->outputstyle_ = Matrix; // how to GET the outputstyle of M?
	this->mat_ = new U [this->size_];
*/	elements_ = 0;
	for (size_t i=0; i<this->size_; i++) {
		if ( std::abs(this->mat_[i] = M[i]) >= (double)SPTHRESHOLD )
			elements_++;
	}
	positions_ = new size_t [elements_];
	int counter=0;
	for (size_t i=0; i<this->size_; i++){
		if ( std::abs(this->mat_[i]) >= (double)SPTHRESHOLD) {
			positions_[counter] = i;
			counter++;
		}	
	}
}
template <class U>
inline spMatrix<U>::spMatrix(const spMatrix<U>& M) : matrix<U>::matrix(M) {
/*	this->rows_=M.GetRows();	// TODO: shouldnt I be able to access M.cols_ directly?
	this->cols_=M.GetColumns();
	this->size_=M.size();
	this->LD_=M.GetLD();
	this->outputstyle_ = Matrix; // how to GET the outputstyle of M?
	this->mat_ = new U [this->size_];
	*/
	elements_ = M.GetElements();
	positions_ = new size_t [elements_];
	memcpy(positions_, M.GetPositions(), sizeof(positions_) );
	for (size_t k=0; k<M.elements_; k++){
		this->mat_[k]=M.mat_[k];
	}	
}
/*
template <class V>
inline spMatrix<V>::~spMatrix() {
	//destructor, deletes the matrix from memory when we leave the scope. mat_ is taken care of, ~matrix gets called afterwards!
	if (elements_ != 0)
		delete[] (positions_);  
}
*/
template <class V>
void spMatrix<V>::Update(){
	if (elements_ != 0) {
		delete[] positions_;	
		elements_ = 0;	
	}
	for (size_t i=0; i<this->size_; i++) {
		if ( abs(this->mat_[i]) >= (double)SPTHRESHOLD )
			elements_++; 
		/* TODO: should we set others to zero?
		else
			this->mat_[i]=0;*/
	}
	positions_ = new size_t [elements_];
	int counter=0;
	for (size_t i=0; i<this->size_; i++){
		if ( abs(this->mat_[i]) >= (double)SPTHRESHOLD) {
			positions_[counter] = i;
			counter++;
		}	
	}

}
// Multiplication (does not catch an error of a W = real and V being complex)
template <class V, class W>
inline spMatrix<V> operator*(const W& beta, const spMatrix<V>& A){//multiplication by a scalar beta*A
	spMatrix<V> result( A );
	size_t p;
	for (size_t k=0; k<A.elements_; k++){
		p=A.positions_[k];
		result[p]=beta*A[p];
	}
	return result;
}
template <class V, class W>
inline spMatrix<V> operator*(const spMatrix<V>& A, const W& beta){//multiplication by a scalar A*beta
	return beta*A; // 'zero' elements that are really tiny but !=0 will not be in positions_ even if beta lifts them over threshold
}
template <class V>
inline matrix<V> operator*(const spMatrix<V>& A, const matrix<V>& B){
	matrix<V> result( A.GetRows(), B.GetColumns() );
	size_t i, j, p, Arows=A.rows_;	
	for (size_t k=0; k<A.elements_; k++){
		p=A.positions_[k];
		i=p%Arows;
		j=p/Arows;
		for (size_t l=0; l<B.GetColumns(); l++)
			result(i,l)+=A(i,j)*B(j,l);	// TODO: use arrays instead of (i,j) operators to speed up?
	}	
	return result;
}
template <class V>
inline matrix<V> operator*(const matrix<V>& A, const spMatrix<V>& B){
	matrix<V> result( A.GetRows(), B.GetColumns() );
	size_t i, j, p, Brows=B.rows_;	
	for (size_t k=0; k<B.elements_; k++){
		p=B.positions_[k];
		i=p%Brows;
		j=p/Brows;
		for (size_t l=0; l<A.GetRows(); l++)
			result(l,j)+=A(l,i)*B(i,j); // TODO: use arrays instead of (i,j) operators to speed up?
	}	
	return result;
}
template <class V>
inline spMatrix<V> operator*(const spMatrix<V>& A, const spMatrix<V>& B){//sparse*sparse 
	spMatrix<V> result( A.GetRows(), B.GetColumns() );
	size_t i, j, I, J, p, Arows=A.rows_, Brows=B.rows_;	
	size_t resElemCount=0; // to fill result's position array properly and set result.elements later
	size_t * pos;
	pos = new size_t [ min(A.elements_*B.elements_, A.rows_*B.cols_ ) ];
	for (size_t k=0; k<A.elements_; k++){
		p=A.positions_[k];
		i=p%Arows;
		j=p/Arows;
		for (size_t l=0; l<B.elements_; l++){
			I=	B.positions_[l]%Brows;	
			if (j==I ){
				J = B.positions_[l]/Brows;
				result(i,J)+=A.mat_[k]*B.mat_[l]; // TODO: use arrays instead of (i,j) operators to speed up?
				pos[resElemCount]=i*J; //result.positions_[resElemCount]=i*J;
				resElemCount++;
			}
		}
	}
	result.elements_=resElemCount;
	result.positions_ = new size_t [resElemCount];
	memcpy(result.positions_, pos, sizeof(result.positions_));
	return result;
}
/*							guess i don't really need this in the end

template <class V>
inline matrix<V> operator*(const spMatrix<V>& A, const spMatrix<V>& B){//dense = sparse*sparse. doesnt make terribly much sense
	// same as spMatrix = spMatrix*spMatrix, just w/o setting position and elements arrays
	matrix<V> result( A.GetRows(), B.GetColumns() );
	size_t i, j, I, J, p, Arows=A.rows_, Brows=B.rows_;	
	for (size_t k=0; k<A.elements_; k++){
		p=A.positions_[k];
		i=p%Arows;
		j=p/Arows;
		for (size_t l=0; l<B.elements_; l++){
			I=	B.positions_[l]%Brows;	
			if (j==I){
				J = B.positions_[l]/Brows;
				result(i,J)+=A.mat_[k]*B.mat_[l];
			}
		}
	}
	return result;
}*/

/*
	//scan B.positions_ once and save information to avoid scanning for each of A's entries
	size_t* pos_cmo=B.positions_;
	size_t pos_rmo[B.elements_];
	size_t nrInCol[B.cols_], nrInRow[B.rows_];
	// fill up...
	size_t i,j, counter=0;
	while(counter<B.rows_){
		for (int k=0; k<B.elements_;k++){
			i=		// this doesnt make sense, does it?
	}
		
	return res; 
*/

template <class V>
size_t spMatrix<V>::GetElements() const{
	return elements_;
}
template <class V>
size_t* spMatrix<V>::GetPositions() const{ // well, it's not really const, is it?
	return positions_;
}

template <class V>
double spMatrix<V>::GiveRatio(){
	double s = sqrt(this->size_);
	double el = sqrt(elements_);
	return el/s;
}


#endif /* spMatrix_h */

