/*
Name: TestReal.cpp
Author: Jay Gambetta

Dependences: MatrixOperations.hpp
Brief Discription: This program test my matrix operation class on Real Matricies
Limitations: None

Version History
	v0: February 12, 2008.
	v1: January 3rd, 2009. Updated to support my new matrix class

*/
#include <MatrixOperations.hpp>
using namespace std;

int main(){
	matrix<double> a(2,2), b(2,2), c(2,2), d(4,4); // we can also declare matrices of type int, float, double etc.
	a.SetOutputStyle(Matrix);
	b.SetOutputStyle(Matrix);
	c.SetOutputStyle(Matrix);
	d.SetOutputStyle(Matrix);
	
	cout << "Enter real matrix a: " << endl;
	cin >> a;
	cout << "a is: " << endl;
	cout << a << endl;
	cout << "Enter real matrix b: " << endl;
	cin >> b;	
	cout << "b is: " << endl;
	cout << b << endl;
	
	cout << "A is square " << MOs::IsSquare(a) << endl;
	cout << "A is diagonal " << MOs::IsDiagonal(a) << endl;
	cout << "A is scalar " << MOs::IsScalar(a) << endl;
	cout << "A is Null " << MOs::IsNull(a) << endl;
	cout << "A is Symmetric " << MOs::IsSymmetric(a) << endl;
	cout << "A is skew symmetric " << MOs::IsSkewSymmetric(a) << endl;
	cout << "A is uppertriangular " << MOs::IsUpperTriangular(a) << endl;
	cout << "A is lower triangular " << MOs::IsLowerTriangular(a) << endl;
	
	
	double alpha;
	alpha = MOs::Trace(a);
    cout << "The trace of a is: "<< alpha << endl; 
	
	c = MOs::Transpose(a);
	cout << endl << "Result of Transpose(a):" << endl;
	cout << c << endl;
		
	d=MOs::TensorProduct(a,b);
	cout << endl << " a tensor b" << endl;
	cout << d << endl;
	
	d=MOs::TensorProduct(b,a);
	cout << endl << " b tensor a" << endl;
	cout << d << endl;
	
	return 0;
}