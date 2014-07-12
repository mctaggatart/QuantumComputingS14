/*
Name: TestBlas.cpp
Author: Jay Gambetta

Dependences: None  
Brief Discription: XXXX
Limitations: None

Version History
	v0: January 12, 2009.

*/
#include <iostream> 
#include <Accelerate/Accelerate.h>
using namespace std;

int main (void)
{
	const size_t rows=3, size2=rows*rows;
	
	complex<double> A[size2], B[size2], C[size2];
	
	A[0*rows+0]=3.1;  A[1*rows+0]=1.3;  A[2*rows+0]=-5.7;	/* matrix A */
	A[0*rows+1]=1.0;  A[1*rows+1]=-6.9; A[2*rows+1]=5.8;	
	A[0*rows+2]=3.4;  A[1*rows+2]=7.2;  A[2*rows+2]=-8.8;
	
	B[0*rows+0]=3.1;  B[1*rows+0]=1.3;  B[2*rows+0]=-5.7;	/* matrix B */
	B[0*rows+1]=1.0;  B[1*rows+1]=-6.9; B[2*rows+1]=5.8;	
	B[0*rows+2]=3.4;  B[1*rows+2]=7.2;  B[2*rows+2]=-8.8;
	
	complex<double> alpha =1.0, beta =0.0;
	
	complex<double>* Apt= A;
	complex<double>* Bpt= B;
	complex<double>* Cpt= C;	
		
	size_t LDA = rows;
	size_t LDB = rows;
	size_t LDC = rows;
	cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, rows, rows, rows, &alpha, Apt, LDA, Bpt, LDB, &beta, Cpt, LDC);
	for(size_t i = 0; i < rows; ++i)
	{
		for(size_t j = 0; j < rows; ++j)
		{
			cout << C[j*rows+i] << "\t";
		}
		cout << endl;
	}
	
	return 0;

}
// 
// 
// #include<iostream>
// #include<complex>
// 
// 
// // define complex type
// typedef std::complex<double> zomplex;
// 
// const zomplex Im = zomplex(0.0, 1.0); 
// 
// /* ZGEMM */
// extern "C" void zgemm_(const char*, const char*, const int*, const int*, const int*,
// const double*, void*, const int*, void*, const int*, const double*, void*, const int*);
// 
// 
// const int m = 300; // the dimension of the matrices in the test
// const int n = m*m; // the number of elements in a matrix
// const int ntimes = 100; // the number of times we do the computation C = AxB;
// 
// 
// int main(void)
// {
// zomplex* A = new zomplex[m*m];
// zomplex* B = new zomplex[m*m];
// zomplex* C = new zomplex[m*m];

// 
// int i;
// for(i=0; i<n; i++){
//    A[i]=zomplex((double) (rand()-1)/RAND_MAX, (double) (rand()-1)/RAND_MAX);
//    B[i]=zomplex((double) (rand()-1)/RAND_MAX, (double) (rand()-1)/RAND_MAX);
//    C[i]=zomplex((double) (rand()-1)/RAND_MAX, (double) (rand()-1)/RAND_MAX);
// }
// 
// const char trans[] = {'N'};
// const double alpha[] = {1.4, 0.5};
// const double beta[] = {2.3, 0.2};
// 
// 
// for (i = 0; i < ntimes; ++i)
// {
// zgemm_(trans, trans, &m, &m, &m, alpha, A, &m, B, &m, beta, C, &m);
// }
// 
// 
// delete[] A;
// delete[] B;
// delete[] C;
// }

// CXX = g++ 
// CC = g++ 
// CXXFLAGS = -Wall -W 
// #LOADLIBES= -lf77blas -latlas -lg2c 
// LOADLIBES= -lgotoblas -lg2c -lpthread 
// #LOADLIBES= -lblas -lg2c 
// #LOADLIBES= -lacml -lg2c 
// 
// TARGETS = zgemm_test 
// 
// all: $(TARGETS) 
// 
// 
// .PHONY: clean 
// clean: 
// rm -rf *.o *~ *orig *exe $(TARGETS) 
