/*
Name: TestReal.cpp
Author: Jay Gambetta

Dependences: None  
Brief Discription: Tests the matrix class and its friends operators for the template double
Limitations: None

Version History
	v0: January  2, 2009.

*/
using namespace std;
#include "Matrix.hpp"

int main (int argc, char const *argv[])
{
	cout << "Running program " << argv[0] << endl;
	int rows=3;
	
	//Testing Assignement and input output
	matrix<complex<double> > A(rows,rows), C(rows*rows), D(rows,rows), E(rows,rows);
	cout << "rows of C " << C.GetRows() << " columns of C " << C.GetColumns() << endl; 
	
	cout << A << endl;
	cout << "Enter Matrix A" << endl;
	cin >> A;
	cout << A << endl;
	cout << "Printing as a Matrix" << endl;
	A.SetOutputStyle(Matrix);
	C.SetOutputStyle(Matrix);
	D.SetOutputStyle(Matrix);
	E.SetOutputStyle(Matrix);
	cout << A << endl;
	cout << "lets copy it to B" << endl;
	matrix<complex<double> > B(A);
	cout << B << endl;
	B.SetOutputStyle(Matrix);
	
	//assinging elements using matrix representation 
	cout << "lets change B_10" << endl;
	cout << "enter a value for B_10:" << endl;
	cin >> B(1,0); 
	cout << B << endl;
	// assigning elements using vector representation
	cout << "lets change B_10" << endl;
	cout << "enter a value for B_10:" << endl;
	cin >> B[1]; 
	cout << B << endl;
	cout << "To show that A was not changed" << endl;
	cout << A << endl;
	
	cout << endl << "The number of rows and columns of A are: " << A.GetRows() << " x "<< A.GetColumns() << endl;
	
	cout << "Re Enter Matrix B" << endl;
	cin >> B;
	cout << B << endl;
	
	
	//Testing addition / substraction
	cout << " ------------------------------------ Testing operations now --------------------------" << endl;

	
	C = A + B;
	cout << endl << "Result of A + B:" << endl << C << endl;
	
	C = A - B;
	cout << endl << "Result of A - b:" << endl << C << endl;
	
	C = A * B;
	cout << endl << "Result of A * B (which is now C) is:" << endl << C << endl;
	 
	cout << endl << "Result of A + (B*C) is:" << endl << (A+B*C) << endl;
	
	C = A*5.0;
	cout << endl << "Result of A*5 is:" << endl << C << endl;

	for(int i=0; i<rows; i++){
		for(int j=0; j<rows; j++){
			D(i,j)=complex<double>(5,4);
			E(i,j)=complex<double>(8,3);
		}
	}

	D += A;
	cout << endl << "Result of D + A using +=  is:" << endl << D<< endl;
	
	E-= A;
	cout << endl << "Result of E - A using -=  is:" << endl << E << endl;

	matrix<double>  A2(rows,rows), B2(rows,rows);
	
	for(int i=0; i<rows; i++){
		for(int j=0; j<rows; j++){
			A2(i,j)=2;
			B2(i,j)=3;
		}
	}
	cout << "recap A is" << endl;
	cout << A << endl;
	C = A * B2;
	cout << endl << "Result of A * B2 is:" << endl << C << endl;
	
	C = A2 * A;
	cout << endl << "Result of A2 * A is:" << endl << C << endl;
	
	return 0;
}
