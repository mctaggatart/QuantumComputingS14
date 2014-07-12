/*
Name: CavityOperators.cpp
Author: Jay Gambetta

Dependences: MatrixOperations.hpp
Brief Discription: Test the Destroy and the identity calls from matrix operations
Limitations: None

Version History
	v0: October 10, 2008.

*/
#include <MatrixOperations.hpp>
using namespace std;

int main ()
{
	matrix<complex<double> > a(10,10), ad(10,10), id(10,10), n(10,10);
	MOs::Destroy(a);
	ad =MOs::Dagger(a);
	n = ad*a;
	MOs::Identity(id);
	
	n.SetOutputStyle(Matrix);
	a.SetOutputStyle(Matrix);
	ad.SetOutputStyle(Matrix);
	id.SetOutputStyle(Matrix);
	
	cout << "n is: " << endl;
	cout << n << endl;
	cout << "a is: " << endl;
	cout << a << endl;
	cout << "ad is: " << endl;
	cout << ad << endl;
	cout << "id is: " << endl;
	cout << id << endl;

	return 0;
}
