/*
Name: PauliGellMann.cpp
Author: Jay Gambetta

Dependences: MatrixOperations.hpp
Brief Discription: Test the Puali and the GellMan calls from matrix operations
Limitations: None

Version History
	v0: October 10, 2008.

*/
#include <MatrixOperations.hpp>
using namespace std;

int main ()
{
	matrix<complex<double> > s1(2,2), s2(2,2), s3(2,2);
	matrix<complex<double> > l1(3,3), l2(3,3), l3(3,3), l4(3,3), l5(3,3), l6(3,3), l7(3,3), l8(3,3);
	s1.SetOutputStyle(Matrix);
	s2.SetOutputStyle(Matrix);
	s3.SetOutputStyle(Matrix);
	l1.SetOutputStyle(Matrix);
	l2.SetOutputStyle(Matrix);
	l3.SetOutputStyle(Matrix);
	l4.SetOutputStyle(Matrix);
	l5.SetOutputStyle(Matrix);
	l6.SetOutputStyle(Matrix);
	l7.SetOutputStyle(Matrix);
	l8.SetOutputStyle(Matrix);
	
	MOs::Pauli(s1, s2, s3);
	
	cout << "Pauli operators "<< endl;
	cout << "s1: " << endl;
	cout << s1 << endl;
	cout << "s2: " << endl;
	cout << s2 << endl;
	cout << "s3: " << endl;
	cout << s3 << endl;
	
	MOs::GellMann(l1, l2, l3, l4, l5, l6, l7, l8);
	
	cout << "Gell Mann operators "<< endl;
	cout << "l1: " << endl;
	cout << l1 << endl;
	cout << "l2: " << endl;
	cout << l2 << endl;
	cout << "l3: " << endl;
	cout << l3 << endl;
	cout << "l4: " << endl;
	cout << l4 << endl;
	cout << "l5: " << endl;
	cout << l5 << endl;
	cout << "l6: " << endl;
	cout << l6 << endl;
	cout << "l7: " << endl;
	cout << l7 << endl;
	cout << "l8: " << endl;
	cout << l8 << endl;
	
	
	
	
	return 0;
}
