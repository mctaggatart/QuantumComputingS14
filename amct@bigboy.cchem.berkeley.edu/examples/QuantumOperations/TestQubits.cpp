/*
Name: TestQubits.cpp
Author: Jay Gambetta

Dependences: MatrixOperations.hpp, QuantumOperations.hpp
Brief Discription: This program test the cavity quantum states and expectation value
Limitations: None

Version History
	v0: January 4th, 2009

*/
#include <MatrixOperations.hpp>
#include <QuantumOperations.hpp>
#include <UsefullFunctions.hpp>
using namespace std;



int main(){
	
	matrix<complex<double> > sx(2,2), sy(2,2), sz(2,2), si(2,2);
	MOs::Pauli(sx, sy, sz);
	MOs::Identity(si);
	
	matrix<complex<double> > sii(MOs::TensorProduct(si,si));
	matrix<complex<double> > sxi(MOs::TensorProduct(sx,si));
	matrix<complex<double> > syi(MOs::TensorProduct(sy,si));
	matrix<complex<double> > szi(MOs::TensorProduct(sz,si));
	matrix<complex<double> > six(MOs::TensorProduct(si,sx));
	matrix<complex<double> > sxx(MOs::TensorProduct(sx,sx));
	matrix<complex<double> > syx(MOs::TensorProduct(sy,sx));
	matrix<complex<double> > szx(MOs::TensorProduct(sz,sx));
	matrix<complex<double> > siy(MOs::TensorProduct(si,sy));
	matrix<complex<double> > sxy(MOs::TensorProduct(sx,sy));
	matrix<complex<double> > syy(MOs::TensorProduct(sy,sy));	
	matrix<complex<double> > szy(MOs::TensorProduct(sz,sy));
	matrix<complex<double> > siz(MOs::TensorProduct(si,sz));
	matrix<complex<double> > sxz(MOs::TensorProduct(sx,sz));
	matrix<complex<double> > syz(MOs::TensorProduct(sy,sz));
	matrix<complex<double> > szz(MOs::TensorProduct(sz,sz));
	
	sxi.SetOutputStyle(Matrix);
	syi.SetOutputStyle(Matrix);
	szi.SetOutputStyle(Matrix);
	six.SetOutputStyle(Matrix);
	siy.SetOutputStyle(Matrix);
	siz.SetOutputStyle(Matrix);
	
	
	cout << "sxi is: " << endl;
	cout << sxi << endl;
	cout << "syi is: " << endl;
	cout << syi << endl;
	cout << "szi is: " << endl;
	cout << szi << endl;
	cout << "six is: " << endl;
	cout << six << endl;
	cout << "siy is: " << endl;
	cout << siy << endl;
	cout << "siz is: " << endl;
	cout << siz << endl;
	
	matrix<complex<double> > rho1(2,2), rho2(2,2), Rho(4,4);
	
	rho1.SetOutputStyle(Matrix);
	rho2.SetOutputStyle(Matrix);
	Rho.SetOutputStyle(Matrix);
	
	QOs::QubitState(rho1,0,1,0);
	cout << "rho 1 is: " << endl;
	cout << rho1 << endl;
	cout << "sy " << QOs::Expectation(sy,rho1) << endl;
	
	QOs::QubitState(rho2,1,0,0);
	cout << "rho 2 is: " << endl;
	cout << rho2 << endl;
	cout << "sx " << QOs::Expectation(sx,rho2) << endl;
	
	cout << "rho total is: " << endl;
	Rho=MOs::TensorProduct(rho1,rho2);
	cout << Rho << endl;
	
	cout << "-----------------------------" << endl;
	cout << "sxi " << QOs::Expectation(sxi,Rho) << endl;
	cout << "syi " << QOs::Expectation(syi,Rho) << endl;
	cout << "szi " << QOs::Expectation(szi,Rho) << endl;
	cout << "six " << QOs::Expectation(six,Rho) << endl;
	cout << "siy " << QOs::Expectation(siy,Rho) << endl;
	cout << "siz " << QOs::Expectation(siz,Rho) << endl;
	
	cout << QOs::Concurrence(Rho) << endl;
	cout << QOs::Purity(Rho) << endl;
	
	MOs::Null(Rho);
	Rho(0,0)=0.5; 
	Rho(3,3)=0.5;
	Rho(0,3)=0.5;
	Rho(3,0)=0.5;
	
	cout << "Rho Bell state: " << endl;
	cout << Rho << endl;
	cout << QOs::Concurrence(Rho) << endl;
	cout << QOs::Purity(Rho) << endl;
	
	
	return 0;
}
