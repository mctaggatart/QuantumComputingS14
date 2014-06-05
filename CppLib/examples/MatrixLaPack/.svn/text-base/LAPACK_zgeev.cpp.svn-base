/*
Name: TestLAPACK4.cpp
Author: Jay Gambetta

Dependences: Lapack, Matrix.hpp
Brief Discription: Finds the right, left eigenvectors and eigenvalues of a complex matrix
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
	matrix<complex<double> > A(dim,dim), VL(dim,dim), VR(dim,dim);
	A.SetOutputStyle(Matrix);
	VL.SetOutputStyle(Matrix);
	VR.SetOutputStyle(Matrix);

	cout << "enter the complex matrix" << endl;
	cin >> A;

	cout << "the matrix you entered is " << endl;
	cout << A << endl;


	//  SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
	//options
	size_t N = dim, LDA=dim, LDVL=dim, LDVR=dim;
	int LWORK=2*dim, ok;
	char JOBVL = 'V';
	char JOBVR = 'V';
	complex<double>* Apt = A.GetMat(); // assigns Apt to point to the matrix array contained in A;
	complex<double>* VLpt = VL.GetMat(); 
	complex<double>* VRpt = VR.GetMat();
	complex<double> *W, *WORK;
	double *RWORK;

	W =new complex<double>[LDA];
	WORK = new complex<double>[LWORK];
	RWORK = new double[LWORK]; 

	zgeev_(&JOBVL,&JOBVR,&N,Apt,&LDA,W,VLpt,&LDVL,VRpt,&LDVR,WORK,&LWORK,RWORK,&ok);

	if (ok==0){
		cout << "The eigenvalues are " << endl;
		for (size_t i=0; i<N; i++){
			cout << W[i] << endl;
	   	}
		cout << "The right eigenvector matrix is "<< endl;
		cout << VR << endl;
		cout << "The left eigenvector matrix is " << endl;
		cout << VL << endl;
	
	}
	else cout << "An error occured" << endl;

	delete [] W;
	delete [] WORK;
	delete [] RWORK;
	
	return 0;
}

