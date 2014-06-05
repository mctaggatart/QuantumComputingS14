/*
Name: TestCavity.cpp
Author: Jay Gambetta

Dependences: MatrixOperations.hpp, QuantumOperations.hpp
Brief Discription: This program test the superoperators
Limitations: None

Version History
	v0: January 4th, 2009

*/
#include <MatrixOperations.hpp>
#include <QuantumOperations.hpp>
using namespace std;
int main(){
	
	unsigned int dim = 4;//Hilbert Space 
	matrix<complex<double> > A(dim,dim), rho(dim,dim);
	
	for(size_t i = 0; i < dim; i++){
		A(i,i)=i;
		for (size_t j=0; j < i; j++){
			A(i,j)=complex<double>(i,j);
			A(j,i) = conj(A(i,j));
		}
	}
	
	for(size_t i = 0; i < dim; i++){
		rho(i,i)=i;
		for (size_t j=0; j < i; j++){
			rho(i,j)=complex<double>(i,j);
			rho(j,i) = conj(A(i,j));
		}
	}
	A.SetOutputStyle(Matrix);
	rho.SetOutputStyle(Matrix);
	cout << "A is" << endl;
	cout << A << endl;
	cout << "rho is " << endl;
	cout << rho << endl;
	
	matrix<complex<double> > c(dim,dim);
	c.SetOutputStyle(Matrix);
	cout << "Testing Jump" << endl; 
	c=QOs::SuperJump(A,rho);
	cout << c << endl;
	cout << "Testing Homo" << endl; 
	c=QOs::SuperHomo(A,rho);
	cout << c << endl;
	cout << "Testing Back" << endl; 
	c=QOs::SuperHomo(A+MOs::Dagger(A),rho);
	cout << c << endl;
	cout << "Testing Damping" << endl; 
	c=QOs::SuperDamp(A,rho);
	cout << c << endl;
	cout << "Testing Anti" << endl; 
	c=QOs::SuperAnti(A,rho);
	cout << c << endl;
	cout << "Testing Ham" << endl; 
	c=QOs::SuperHami(A,rho);
	cout << c << endl;
	
	
	
	return 0;
}