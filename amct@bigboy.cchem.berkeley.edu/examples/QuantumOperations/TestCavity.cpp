/*
Name: TestCavity.cpp
Author: Jay Gambetta

Dependences: MatrixOperations.hpp, QuantumOperations.hpp
Brief Discription: This program test the cavity quantum states and expectation value
Limitations: None

Version History
	v0: January 4th, 2009

*/
#include <MatrixOperations.hpp>
#include <QuantumOperations.hpp>
using namespace std;



int main(){
	
	matrix<complex<double> > a(10,10), ad(10,10), rho(10,10), n(10,10);
	MOs::Destroy(a);
	ad =MOs::Dagger(a);
	n = ad*a;
	
	n.SetOutputStyle(Matrix);
	a.SetOutputStyle(Matrix);
	ad.SetOutputStyle(Matrix);
	rho.SetOutputStyle(Matrix);


	cout << "n is: " << endl;
	cout << n << endl;
	cout << "a is: " << endl;
	cout << a << endl;
	cout << "ad is: " << endl;
	cout << ad << endl;
	
	size_t m =3;
	QOs::FockState(rho,m);
	cout << "rho is: " << endl;
	cout << rho << endl;
	cout << "Testing Expectations" << endl;
	complex<double> c;
	c = QOs::Expectation(a,rho);
	cout << c << endl;
	c = QOs::Expectation(n,rho);
	cout << c << endl;
	
	
	return 0;
}