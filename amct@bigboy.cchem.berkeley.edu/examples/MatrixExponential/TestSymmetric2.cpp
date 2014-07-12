/*
Name: TestHermitian.cpp
Author: Jay Gambetta

Dependences: MatrixExponential.hpp  
Brief Discription: An Example showing how to take the expondential of a Hermitian Matrix
Limitations: None

Version History
	v0: October  9, 2008.
	v1: Janurary 4th, 2008 -- changes to support new matrix class

*/
#include <MatrixExponential.hpp>
using namespace std;

int main ()
{
	size_t dim=4;
	
	matrix<double> A(dim,dim), B(dim,dim);
	A.SetOutputStyle(Matrix);
	B.SetOutputStyle(Matrix);

	
	for(size_t i = 0; i < dim; i++){
		A(i,i)=i;
		for (size_t j=0; j < i; j++){
			A(i,j)=i+j;
			A(j,i) = A(i,j);
		}
	}
	
	matrix<double> C(A);
	cout << "The matrix you entered is "<< endl;
	cout << A << endl;
	
	cout << "The Exp of A is " << endl;
	B = ExpM::EigenMethod(A);
	cout <<  B << endl;
	cout << "The matrix you entered is "<< endl;
	cout << C << endl;
	cout << "The Exp of 2A is " << endl;
	B = ExpM::EigenMethod(C,2.0);
	cout <<  B << endl;

	return 0;
}

