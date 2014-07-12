/*
Name: Matrix.hpp
Author: Jay Gambetta

Dependences: UtilityFunctions.hpp, BLAS
Brief Discription: This is my Matrix class. It works only with real/complex matrices and stores the entries in 
	column-major-order. In column-major storage the columns are stored one after the other. The linear 
	offset p from the beginning of the array to any given element A(i,j) can then be computed as:
		p = j*Nrows+i
	where Nrows is the number of rows in the matrix. If i wanted to know the i and j associtated with a given p then i woul duse 
		i=p %Nrows
		j= floor(p/Nrows)
		
	multiplication is done with the C wrapper of the fortran blas library.
	
Limitations: Sparse Matrices are still not included.

Version History
	v0: May 27th, 2007. Started the class file
	v1:	February 12th, 2008. Change the name from MatrixReal.h to RMc.h to stand for real matrix class
	v2: June 15th, 2008. Decided to use Column-major order as this is what lapack uses.
	v3: November 21st, 2008. Changed the name to Matrix.h and combined both real and complex matricies by templates
	v4: Janurary 3rd, 2009. Change matrix multiplication to use the c wrapper for the blas library under macs Accelerate platform.
	v5: Feburary 2nd, 2009 made it use a blas which can be implement
*/

#ifndef Matrix_h
#define Matrix_h
#include "UtilityFunctions.hpp"
#include "Headers.hpp"


enum OutputStyle {Column, List, Matrix};


template <class T>
class matrix 
{	
	//friend functions get to use the private varibles of the class as well as have different classes as inputs
	template <class S>
	friend std::ostream& operator << (std::ostream& output, const matrix<S>& A);//overloading << to output a matrix
	template <class S>
	friend std::istream& operator >> (std::istream& input, const matrix<S>& A);//overloading >> to read in a matrix
	
	// Multiplication (does not catch an error of a S1 = real and S2 being complex)
	template <class S1, class S2>
	friend matrix<S1> operator*(const S2& beta, const matrix<S1>& A);//multiplication by a scalar beta*A
	template <class S1, class S2>
	friend matrix<S1> operator*(const matrix<S1>& A, const S2& beta);//multiplication by a scalar A*beta
	//Real
	friend matrix<double> operator*(const matrix<double>& A, const matrix<double>&B); //real matrix multiplication A*B
	// //Complex
	friend matrix<std::complex<double> > operator*(const matrix<std::complex<double> >& A, const matrix<std::complex<double> >& B);//complex matrix multplication A*B
	friend matrix<std::complex<double> > operator*(const matrix<double>& A, const matrix<std::complex<double> >& B);//real-complex matrix multplication A*B
	friend matrix<std::complex<double> > operator*(const matrix<std::complex<double> >& A, const matrix<double>& B);//real-complex matrix multplication A*B
		
public:
	//Construtors
	matrix();
	matrix(size_t  rows, size_t cols);	
	matrix(size_t  size);	//Makes a square matrix and rows = sqrt(size) columns = sqrt(dims)
	matrix(const matrix<T>& m);
	matrix(const matrix<T>& m, const char uplo); 
	
	//Initialize an empty matrix() to matrix(size_t  rows, size_t cols)
	void initialize(size_t  rows, size_t cols);
	//Clear used memory
	void clear();
			
	//Destructor
 	virtual ~matrix(); 
		
	//Assignment operator
	matrix<T>& operator=(const matrix<T>& m);
	template <class S>
	matrix<T>& operator=(const matrix<S>& m); //Still would like to have real assigend by complex -- take real part
		
	//Addressing elements by vector representation
	T& operator[](size_t element);
	T operator[](size_t element) const;
	//Addressing elements by matrix representation
	T& operator()(size_t row, size_t col);
	T operator()(size_t row, size_t col) const;
	
	//overloading functions. 
	matrix<T> operator+(const matrix<T>& A);
	matrix<T> operator-(const matrix<T>& A);
	matrix<T> operator+(const matrix<T>& A) const;
	matrix<T> operator-(const matrix<T>& A) const;
	matrix<T>& operator+=(const matrix<T>& A);
	matrix<T>& operator-=(const matrix<T>& A);
		
	//Member Functions
	size_t GetColumns() const;//gives the number of columns
	size_t GetRows() const;//gives the number of rows
	size_t GetLD() const;//gives the leading dimension -- number of rows
	size_t size() const;//gives the size of the underlying vector
	void resize(size_t row, size_t col);//sets the size of the underlying vector	
	matrix<T>& SetOutputStyle(enum OutputStyle outputstyle);//sets the style the matrix is display by <<
	T* GetMat() const;//gives you the address of element 0 then *(c+i) gives you the ith element
	
protected:
	size_t rows_, cols_, size_, LD_;
	//rows_ and cols_ are the rows and columns of the matrix
	//size_ = rows*colums dimensions of the vector representation   
	//LD is the leading dimeonsion and for Column major order is in general eqaul to rows
	enum OutputStyle outputstyle_;
	//outputstyle_ is the output style used by << 
	T * mat_;
	// the ptr to the vector containing the matrix
};
//------------Member funtions------------------
template <class T>
inline matrix<T>::matrix() : rows_(0), cols_(0), size_(0), LD_(0), outputstyle_(Column){ mat_=0;}
// constructs an empty matrix using the ....
template <class T>
inline matrix<T>::matrix(size_t rows, size_t cols) : rows_(rows), cols_(cols), size_(rows*cols), LD_(rows), outputstyle_(Column), mat_(new T [size_]){}
//constructs an empty matrix of size rows, cols and sets outputstyle to zero
template <class T>
inline matrix<T>::matrix(size_t dim2) : size_(dim2), outputstyle_(Column), mat_(new T [size_]){
	//constructs a square matrix of dims sqrt(size)*sqrt(size)
	if (size_ == 0)
		rows_ = size_;
	else if (size_ == 1)
		rows_ = size_;
	else if (size_ == 2){
		UFs::MyError("matrix constructor matrix(dim): the number you enterd does not have a interger sqrt");
	}
	else
		for(rows_ = 2; rows_ < size_; ++rows_)
		{
			if (size_ == rows_*rows_)
				break;
			if (rows_*rows_ > size_){
				UFs::MyError("matrix constructor matrix(dim): the number you enterd does not have a interger sqrt");
			}
		}
	
	cols_=rows_;
	LD_=rows_;	
	
}
//constructs an empty matrix of size rows, cols and sets outputstyle to zero
template <class T>
inline matrix<T>::matrix(const matrix<T>& rhs) : rows_(rhs.rows_), cols_(rhs.cols_), size_(rhs.size_), LD_(rows_), outputstyle_(rhs.outputstyle_), mat_(new T [size_]){
	 // Copy constructor, copies the matric to another matrix
	for(size_t p = 0; p < size_; p++){
		mat_[p] = rhs.mat_[p];
	}
}
template <class T>
inline matrix<T>::matrix(const matrix<T>& rhs, const char uplo) : rows_(rhs.rows_), cols_(rhs.cols_), size_(rhs.size_), LD_(rows_), outputstyle_(rhs.outputstyle_), mat_(new T [size_]){
	 // Copy constructor, copies the matric to another matrix but only the upper or lower triangle
	if(uplo=='U')
	{
		for (size_t i=0; i< rows_; i++)
			for(size_t j=i; j<cols_; j++)
				mat_[j*rows_+i]=rhs.mat_[j*rows_+i];
	}
	else if(uplo=='L')
	{
		for (size_t i=0; i< rows_; i++)
			for(size_t j=0; j<=i; j++)
				mat_[j*rows_+i]=rhs.mat_[j*rows_+i];
	}
	else{
		UFs::MyError("matrix copy constructor did not have a valid option -- these are Lower, Upper");
	}		
}
template <class T>
inline void matrix<T>::initialize(size_t  rows, size_t cols){
	if (rows_ != rows || cols_ != cols ) {//if the rows are different size delete re-construct the matrix
		if (mat_ != 0) delete [] (mat_);
		rows_=rows;
		cols_=cols;
		size_=rows_*cols_;
		LD_=rows;
		mat_= new T[size_];
		// std::cout << "Assignement (reassigment), size " << size_ << " rows_ " << rows_ << " cols " << cols_ << " LD "<< LD_ <<  std::endl;
	}
}
template <class T>
inline void matrix<T>::clear(){
	if (!mat_ || !size_) return;
	rows_=cols_=size_=0;
	delete[] (mat_);
	mat_ = 0;
}

template <class T>
inline void matrix<T>::resize(size_t rows, size_t cols){
	if (rows_ == rows && cols_ == cols ) return;
	T * tempmat= new T[size_=rows*cols];
	
	for (size_t i=0;  i<rows; i++)
		for (size_t j=0;  j<cols; j++)
			if(i<rows_ && j<cols_) tempmat[j*rows +i] = mat_[j*rows_ +i];
			else tempmat[j*rows +i]=0.0;
	LD_=rows_=rows;
	cols_=cols;
	delete[] (mat_);
	mat_ = tempmat;	
	
	
}

template <class T>
inline matrix<T>::~matrix(){
	//destructor, deletes the matrix from memory when we leave the scope
	if (mat_ != 0)
		delete[] (mat_);
}
template <class T>
inline matrix<T>& matrix<T>::operator=(const matrix<T>& rhs){
	// overloading the assignement operator
	// postcondition: normal assignment via copying has been performed;
	//		if matrix and rhs were different sizes, matrix
	//		has been resized to match the size of rhs
	//if (&rhs != this){// only do this if they are different.
//	std::cout << "??are I this??" << std::endl;
//	std::cout << rows_ << std::endl;
//	std::cout << cols_ << std::endl;
//	std::cout << rhs << std::endl;
//	std::cout << rhs.GetRows() << std::endl;
//	std::cout << rhs.GetColumns() << std::endl;
		if (rows_ != rhs.rows_ || cols_ != rhs.cols_ ) {//if the rows are different size delete re-construct the matrix
//			std::cout << "aif1" << std::endl;
			if (mat_ != 0) delete [] (mat_);
			rows_=rhs.rows_;
			cols_=rhs.cols_;
			size_=rows_*cols_;
			LD_=rhs.LD_;
			mat_= new T[size_];
//			std::cout << "Assignement (reassigment), size " << size_ << " rows_ " << rows_ << " cols " << cols_ << " LD "<< LD_ <<  std::endl;
		}
//		std::cout << "are past" << std::endl;
		for (size_t p=0; p < size_; p++){// the copying of the matrix
			//std::cout << p << std::endl;
			mat_[p] = T(rhs.mat_[p]);
		}
		//LD_=rhs.LD_;//added as maybe we need to use LD more the just a wrapper for rows in blas/lapack routines;
	//}
	return *this;
}
template <class T> template <class S>
inline matrix<T>& matrix<T>::operator=(const matrix<S>& rhs){
	// overloading the assignement operator
	// postcondition: normal assignment via copying has been performed;
	//		if matrix and rhs were different sizes, matrix
	//		has been resized to match the size of rhs
	//if (&rhs != this){// only do this if they are different.
//	std::cout << "am I this??" << std::endl;
		if (rows_ != rhs.GetRows() || cols_ != rhs.GetColumns() ) {//if the rows are different size delete re-construct the matrix
			if (mat_ != 0) delete [] (mat_);
			rows_=rhs.GetRows();
			cols_=rhs.GetColumns();
			size_=rows_*cols_;
			LD_=rhs.GetLD();
			mat_= new T[size_];
		}
		for (size_t p=0; p < size_; p++){// the copying of the matrix
			mat_[p] = T(rhs[p]);
		}
//		S* c = rhs.GetMat();
//		for (size_t p=0; p < size_; p++){// the copying of the matrix
//			mat_[p] = T(*(c+p));
//		}
		//LD_=rhs.LD_;//added as maybe we need to use LD more the just a wrapper for rows in blas/lapack routines;
	//}
	return *this;
}
template <class T>
inline T& matrix<T>::operator[](size_t p){
	// returns the pth element of the vector representing the matrix,  remember [] is the operator acting on the object matrix (read backwards).
	if (p >=size_){
		UFs::MyError("matrix class operator []: Matrix subscript out of bounds");
	}
	return mat_[p]; 
}
template <class T>
inline T matrix<T>::operator[](size_t p) const{
	// returns the pth element of the vector representing the matrix if the operator is acting on a const
	if (p >=size_){
		UFs::MyError("matrix class operator [] const: Matrix subscript out of bounds");
	}
	return mat_[p]; 
}
template <class T>
inline T& matrix<T>::operator()(size_t i, size_t j){
	//return the data at position i,j, remember () is the operator acting on the object matrix (read backwards).
	if (i >= rows_ || j >= cols_){
		UFs::MyError("matrix class operator (): Matrices subscript out of bounds");
	}
	return mat_[j*rows_+i];
}
template <class T>
inline T matrix<T>::operator()(size_t i, size_t j) const{
	// return the data at i,j if the operator is acting on a const 
	if (i >= rows_ || j >= cols_){
		UFs::MyError("matrix class operator () const: Matrices subscript out of bounds");
	}
	return mat_[j*rows_+i];
}
template <class T>
inline size_t matrix<T>::GetRows() const{
	//returns the rows of the matrix
	return rows_;
}
template <class T>
inline size_t matrix<T>::GetColumns() const{
	//returns the colums of the matrix
	return cols_;
}
template <class T>
inline size_t matrix<T>::GetLD() const{
	//returns the leading dimension 
	return LD_;
}
template <class T>
inline size_t matrix<T>::size() const{
	//returns the size of the underlying vector 
	return size_;
}
template <class T>
inline T* matrix<T>::GetMat() const{
	//returns the ptr for the matrix data
	return mat_;
}
template <class T>
inline matrix<T>& matrix<T>::SetOutputStyle(enum OutputStyle outputstyle){
	//sets the outputstyle 
	outputstyle_=outputstyle;
	return *this;
}
template <class T>
inline matrix<T> matrix<T>::operator+(const matrix<T>& A){
	//overloads the + for matrix addition, can this be more efficient
	if(rows_ != A.rows_ || cols_ != A.cols_)
		UFs::MyError("matrix class operator +: Matrices are not the same size");
	matrix<T> temp(rows_,cols_);
	for (unsigned int p=0; p < size_; p++){			
		temp.mat_[p] = mat_[p] + A.mat_[p];
	}
	return temp;
}
template <class T>
inline matrix<T> matrix<T>::operator-(const matrix<T>& A){
	//overloads the - for matrix substraction, can this be more efficient
	if(rows_ != A.rows_ || cols_ != A.cols_)
		UFs::MyError("matrix class operator -: Matrices are not the same size");
	matrix<T> temp(rows_,cols_);
	for (unsigned int p=0; p < size_; p++){			
		temp.mat_[p] = mat_[p] - A.mat_[p];
	}
	return temp;
}
template <class T>
inline matrix<T> matrix<T>::operator+(const matrix<T>& A) const{
	//overloads the + for matrix addition if it is a const matrix, can this be more efficient
	if(rows_ != A.rows_ || cols_ != A.cols_)
		UFs::MyError("matrix class operator + const: Matrices are not the same size");
	matrix<T> temp(rows_,cols_);
	for (unsigned int p=0; p < size_; p++){			
		temp.mat_[p] = mat_[p] + A.mat_[p];
	}
	return temp;
}
template <class T>
inline matrix<T> matrix<T>::operator-(const matrix<T>& A) const{
	//overloads the - for matrix substraction, can this be more efficient
	if(rows_ != A.rows_ || cols_ != A.cols_)
		UFs::MyError("matrix class operator - const: Matrices are not the same size");
	matrix<T> temp(rows_,cols_);
	for (unsigned int p=0; p < size_; p++){			
		temp.mat_[p] = mat_[p] - A.mat_[p];
	}
	return temp;
}
template <class T>
inline matrix<T>& matrix<T>::operator+=(const matrix<T>& A){
	//overloads the += for matrix addition and assignment, can this be more efficient
	if(rows_ != A.rows_ || cols_ != A.cols_)
		UFs::MyError("matrix class operator +=: Matrices are not the same size");
	for (size_t p=0; p < size_; p++){
		mat_[p] += A.mat_[p];
	}
	return *this;
}
template <class T>
inline matrix<T>& matrix<T>::operator-=(const matrix<T>& A){
	//overloads the -= for matrix subtraction and assignement, can this be more efficient
	if(rows_ != A.rows_ || cols_ != A.cols_)
		UFs::MyError("matrix class operator -=: Matrices are not the same size");
	for (size_t p=0; p < size_; p++){
		mat_[p] -= A.mat_[p];
	}
	return *this;
}
//============Friend Functions ======================================
template <class T>
std::ostream& operator << (std::ostream& output, const matrix<T>& A){
	// overloads << and has three outputs styles
	// Column is in a coulmn vector
	// List is a row of the actual column vector
	// Matrix is the matrix
	if(A.outputstyle_ == List)
	{
		for (size_t p=0; p < A.size_; p++){
			output << A.mat_[p]  << '\t';
		}
	}
	else if(A.outputstyle_ == Matrix)
	{ 
		for (size_t i=0; i<A.rows_; i++){
			for (size_t j=0; j<A.cols_; j++){
				output << A.mat_[j*A.rows_ +i] << "\t";
			}
		output << std::endl;
		}
	}
	
	else if(A.outputstyle_ == Column)
	{
		for (size_t p=0; p < A.size_; p++){
			output << A.mat_[p]  << '\n';
		}
	}
	else{
		UFs::MyError("matrix operator << is not assigned with a valid option -- these are Column, List, Matrix");
	}
	
   	return output;
}
template <class T>
std::istream& operator >> (std::istream& input, const matrix<T>& A){
	// overloads the >> to read in a row into column format 
	for (size_t j=0; j<A.cols_; j++){
		for (size_t i=0; i<A.rows_; i++){
			input >> A.mat_[j*A.rows_ +i];
		}
	}
	return input;
}
template <class S1, class S2>
matrix<S1> operator*(const matrix<S1>& A, const S2& beta){
	//overlads A*beta
	size_t rows= A.rows_, cols = A.cols_;
	matrix<S1> temp(rows,cols);
	for (size_t i=0; i < rows; i++){
		for(size_t j=0; j < cols; j++){
			temp(i,j) = beta*A(i,j);
		}
	}
	return temp;
}
template <class S1, class S2>
matrix<S1> operator*(const S2& beta, const matrix<S1>& A){
	//overloads beta*A
	size_t rows= A.rows_, cols = A.cols_;
	matrix<S1> temp(rows,cols);
	for (size_t i=0; i < rows; i++){
		for(size_t j=0; j < cols; j++){
			temp(i,j) = beta*A(i,j);
		}
	}
	return temp;
}
matrix<double> operator*(const matrix<double>& A, const matrix<double>& B){
	//overloads A*B for real matricies and uses the blas dgemm routine
	//cblas_dgemm(CblasXMajor,op,op,N,M,K,alpha,A,LDA,B,LDB,beta,C,LDC)
	// C-> alpha*op(A)*op(B) +beta C
	matrix<double> C(A.rows_,B.cols_);
	double alpha =1.0, beta =0.0;
	dgemm_(&Trans[0], &Trans[0], &A.rows_, &B.cols_, &A.cols_, &alpha, A.mat_, &A.LD_, B.mat_, &B.LD_, &beta, C.mat_, &C.LD_);
	//cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, A.rows_, B.cols_, A.cols_, 1.0, A.mat_, A.LD_, B.mat_, B.LD_, 0.0, C.mat_, C.LD_);
	std::cout << C << std::endl;
	return C;
}
matrix<std::complex<double> > operator*(const matrix<std::complex<double> >& A, const matrix<std::complex<double> >& B){
	//overloads A*B for complex matricies and uses the blas zgemm routine
	//cblas_zgemm(CblasXMajor,op,op,N,M,K,alpha,A,LDA,B,LDB,beta,C,LDC)
	// C-> alpha*op(A)*op(B) +beta C
	static matrix<std::complex<double> > C;
	C.initialize(A.rows_,B.cols_);
	std::complex<double> alpha =1.0, beta =0.0;
	zgemm_(&Trans[0], &Trans[0], &A.rows_, &B.cols_, &A.cols_, &alpha, A.mat_, &A.LD_, B.mat_, &B.LD_, &beta, C.mat_, &C.LD_);
	// cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, A.rows_, B.cols_, A.cols_, &alpha, A.mat_, A.LD_, B.mat_, B.LD_, &beta, C.mat_, C.LD_);
	return C;
}
matrix<std::complex<double> > operator*(const matrix<double>& A, const matrix<std::complex<double> >& B){
	//overloads A*B for complex matricies and uses the blas zgemm routine
	//cblas_zgemm(CblasXMajor,op,op,N,M,K,alpha,A,LDA,B,LDB,beta,C,LDC)
	// C-> alpha*op(A)*op(B) +beta C
	matrix<std::complex<double> > C(A.rows_,B.cols_), Ac(A.rows_,A.cols_);
	Ac=A;
	std::complex<double> alpha =1.0, beta =0.0;
	zgemm_(&Trans[0], &Trans[0], &Ac.rows_, &B.cols_, &Ac.cols_, &alpha, Ac.mat_, &Ac.LD_, B.mat_, &B.LD_, &beta, C.mat_, &C.LD_);
	// cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, Ac.rows_, B.cols_, Ac.cols_, &alpha, Ac.mat_, Ac.LD_, B.mat_, B.LD_, &beta, C.mat_, C.LD_);
	return C;
}
matrix<std::complex<double> > operator*(const matrix<std::complex<double> >& A, const matrix<double>& B){
	//overloads A*B for complex matricies and uses the blas zgemm routine
	//cblas_zgemm(CblasXMajor,op,op,N,M,K,alpha,A,LDA,B,LDB,beta,C,LDC)
	// C-> alpha*op(A)*op(B) +beta C
	matrix<std::complex<double> > C(A.rows_,B.cols_), Bc(B.rows_,B.cols_);
	Bc=B;
	std::complex<double> alpha =1.0, beta =0.0;
	zgemm_(&Trans[0], &Trans[0], &A.rows_, &Bc.cols_, &A.cols_, &alpha, A.mat_, &A.LD_, Bc.mat_, &Bc.LD_, &beta, C.mat_, &C.LD_);
	// cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, A.rows_, Bc.cols_, A.cols_, &alpha, A.mat_, A.LD_, Bc.mat_, Bc.LD_, &beta, C.mat_, C.LD_);
	return C;
}


#endif /* Matrix_h */

