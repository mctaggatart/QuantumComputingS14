/*
Name: TestLAPACK2.cpp
Author: Jay Gambetta

Dependences: Lapack, Matrix.hpp
Brief Discription: Finds the right, left eigenvectors and eigenvalues of a real matrix
Limitations: None

Version History
	v0: June 10, 2008.
	v1: January 3rd, 2009 - changed it to support my new matrix template -- removed hermitian and symmetric matrix

*/
#include <Matrix.hpp>
using namespace std;

/* finding the eigenvalues of a real matrix */

int main(){

	int dim;
	cout << "enter the size of the matrix "<< endl;
	cin >> dim;
	matrix<double> A(dim,dim), VL(dim,dim), VR(dim,dim);
	A.SetOutputStyle(Matrix);
	VL.SetOutputStyle(Matrix);
	VR.SetOutputStyle(Matrix);

	cout << "enter the matrix" << endl;
	cin >> A;

	cout << "The matrix you entered is "<< endl;
	cout << A << endl;


	//SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
	//options
	size_t N = dim, LDA=dim, LDVL=dim, LDVR=dim;
	int LWORK=4*dim, ok;
	char JOBVL = 'V';
	char JOBVR = 'V';
	double* Apt = A.GetMat(); // assigns Apt to point to the matrix array contained in A;
	double* VLpt = VL.GetMat(); 
	double* VRpt = VR.GetMat();
	double *WR, *WI, *WORK;

	WR =new double[LDA];
	WI =new double[LDA];
	WORK = new double[LWORK];

	dgeev_(&JOBVL,&JOBVR,&N,Apt,&LDA,WR,WI,VLpt,&LDVL,VRpt,&LDVR,WORK, &LWORK, &ok);

	if (ok==0){
		cout << "The eigenvalues are " << endl;
		for (size_t i=0; i<N; i++){
			cout << " real "<< WR[i] << " imag " << WI[i] << endl;
   		}
		cout << "The right eigenvector matrix is "<< endl;
		cout << VR << endl;
		cout << "The left eigenvector matrix is " << endl;
		cout << VL << endl;
	}
	else cout << "An error occured" << endl;

	delete [] WI;
	delete [] WR;
	delete [] WORK;
	
	return 0;
}