/*
Name: TestAssignement.cpp
Author: Jay Gambetta

Dependences: None  
Brief Discription: Tests the matrix class and its friends operators for the template double
Limitations: None

Version History
	v0: January  2, 2009.

*/
using namespace std;
#include "Matrix.hpp"

int main ()
{
	int rows=3;
	
	//Testing Assignement and input output
	matrix<double> A(rows,rows), C(rows*rows), B;
	cout << "rows of A " << A.GetRows() << " columns of A " << A.GetColumns() << endl; 
	cout << "rows of B " << B.GetRows() << " columns of B " << B.GetColumns() << endl; 
	cout << "rows of C " << C.GetRows() << " columns of C " << C.GetColumns() << endl; 
	B=C;
	cout << "rows of B " << B.GetRows() << " columns of B " << B.GetColumns() << endl; 	

	return 0;
}
