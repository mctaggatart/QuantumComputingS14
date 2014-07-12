/*
Name: TestLAPACK1.cpp
Author: Jay Gambetta

Dependences: Lapack, Matrix.hpp
Brief Discription: Solves the equation A x = b using the lapack routines 
Limitations: None

Version History
	v0: June 10, 2008.
	v1: January 3rd, 2009 - changed it to support my new matrix template -- removed hermitian and symmetric matrix
	v2: January 13th, 2009 - romoved dependence on my matrix class so it can be a test for just lapack
	

*/

#include <iostream> 
using namespace std;

extern "C"{
extern void dgesv_(size_t*,size_t*, double*,size_t*,size_t*, double*, size_t*,size_t*);}

int main ()
{	
	
	const size_t size = 3; 	/* dimension of matrix */
	const size_t size2 = size*size;
	size_t c1, c2, pivot[size], ok;
	double b[size];	/* single precision! */
	double A[size2];



	A[0*size+0]=3.1;  A[1*size+0]=1.3;  A[2*size+0]=-5.7;	/* matrix A */
	A[0*size+1]=1.0;  A[1*size+1]=-6.9; A[2*size+1]=5.8;	
	A[0*size+2]=3.4;  A[1*size+2]=7.2;  A[2*size+2]=-8.8;	

	b[0]=-1.3;			/* if you define b as a matrix then you */
	b[1]=-0.1;			/* can solve multiple equations with */
	b[2]=1.8;			/* the same A but different b */ 	
	double* Apt= A;

	c1=size;			/* and put all numbers we want to pass */
	c2=1;    			/* to the routine in variables */

	/* find solution using LAPACK routine SGESV, all the arguments have to */
	/* be pointers and you have to add an underscore to the routine name */
	dgesv_(&c1, &c2, Apt, &c1, pivot, b, &c1, &ok);      

	/*
	 parameters in the order as they appear in the function call
	    order of matrix A, number of right hand sides (b), matrix A,
	    leading dimension of A, array that records pivoting, 
	    result vector b on entry, x on exit, leading dimension of b
	    return value */ 

	for (size_t j=0; j<size; j++) cout << b[j] << endl;	/* print vector x */
	return 0;
}
