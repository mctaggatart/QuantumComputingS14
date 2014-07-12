/*
Name: MatrixExponential.hpp
Author: Jay Gambetta

Dependences: MatrixOperations.hpp, Lapack   
Brief Discription: A header file which uses my matrix class to calculate the expondential of a matrix. 

Limitations: currently I only use the diagonlization method (under LAPACK) but i want extend this to use other methods

Version History
	v0: October  9, 2008.
	v1: January 3rd, 2009 - changed it to support my new matrix template -- removed hermitian and symmetric matrix

*/
#ifndef MatrixExponential_h
#define MatrixExponential_h
#include "MatrixOperations.hpp"

namespace ExpM { //All my matrix operators
	
	matrix<double> EigenMethod(const matrix<double>& A);
	matrix<double> EigenMethod(const matrix<double>& A, const double& x);
	matrix<std::complex<double> > EigenMethod(const matrix<std::complex<double> >& A);									//Z eigenvectors of A, W eigenvalues of A
	matrix<std::complex<double> > EigenMethod(const matrix<std::complex<double> >& Ain, const std::complex<double>& x, matrix<std::complex<double> >* Z=NULL, double *W=NULL );	
	
}
//------------Member funtions------------------
inline matrix<double > ExpM::EigenMethod(const matrix<double >& A){
	
	if (MOs::IsSymmetric(A)==0) UFs::MyError("Routine ExpM::EigenMethod: The input matrix is not Symmetric.");
	size_t N = A.GetRows(), LDA=A.GetLD(), IL, IU, M, LDZ=N;
	double VL, VU, ABSTOL=std::numeric_limits<double>::min(); 
	double *WORK;
 	int *IWORK;
	//finding the optimal work sizes
	int LWORK=-1, LIWORK=-1, ok=0;
	WORK = new double[1];
 	IWORK = new int[1];
	dsyevr_(&Jobz[0],&Range[0],&UpLo[1],&N,NULL,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,NULL,NULL,&LDZ,NULL,WORK,&LWORK,IWORK,&LIWORK,&ok);
	
	//Setting the optimal workspace
	LWORK = (int)WORK[0]  ;
 	LIWORK = IWORK[0];
	
	delete WORK;
	delete IWORK;

	//Running the algorithim 
	matrix<double> Z(LDZ,LDZ), temp(A);
	double *W;
	W =new double[LDA];
	int *ISUPPZ;
	ISUPPZ = new int[2*N]; 
	WORK = new double[LWORK];
	IWORK = new int[LIWORK];
	dsyevr_(&Jobz[0],&Range[0],&UpLo[1],&N,temp.GetMat(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z.GetMat(),&LDZ,ISUPPZ,WORK,&LWORK,IWORK,&LIWORK,&ok);	
	
	// size_t N = A.GetRows(), LDA=A.GetLD();
	// matrix<double> temp(N,N);
	// 
	// int LWORK=3*N-1, ok;
	// 
	// double *W, *WORK;
	// W =new double[N];
	// WORK = new double[LWORK];
	// 
	// 
	// dsyev_(&Jobz[0],&UpLo[1],&N,A.GetMat(),&LDA,W,WORK,&LWORK,&ok);
	if (ok!=0)
		UFs::MyError("Routine ExpM::EigenMethod: An error occured in geting the diagonalization of A");

	//Getting the Dig*Dagger(Z)
	for (size_t i=0; i < N; i++){
		for(size_t j=0; j < N; j++){
			temp(i,j) = std::exp(W[i])*Z(j,i);
		}
	}

	temp = Z*temp;
	delete [] W; 
	delete [] ISUPPZ;
	delete [] WORK;
	delete [] IWORK;
	
	return temp;
}
inline matrix<double > ExpM::EigenMethod(const matrix<double >& A, const double& x){
 
	if (MOs::IsSymmetric(A)==0) UFs::MyError("Routine ExpM::EigenMethod: The input matrix is not Symmetric.");
	size_t N = A.GetRows(), LDA=A.GetLD(), IL, IU, M, LDZ=N;
	double VL, VU, ABSTOL=std::numeric_limits<double>::min(); 
	double *WORK;
 	int *IWORK;
	//finding the optimal work sizes
	int LWORK=-1, LIWORK=-1, ok=0;
	WORK = new double[1];
 	IWORK = new int[1];
	dsyevr_(&Jobz[0],&Range[0],&UpLo[1],&N,NULL,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,NULL,NULL,&LDZ,NULL,WORK,&LWORK,IWORK,&LIWORK,&ok);
	
	//Setting the optimal workspace
	LWORK = (int)WORK[0]  ;
 	LIWORK = IWORK[0];
	
	delete WORK;
	delete IWORK;

	//Running the algorithim 
	matrix<double> *Z, temp(A);
	Z = new matrix<double>(LDZ,LDZ);
	double *W;
	W =new double[LDA];
	int *ISUPPZ;
	ISUPPZ = new int[2*N]; 
	WORK = new double[LWORK];
	IWORK = new int[LIWORK];
	dsyevr_(&Jobz[0],&Range[0],&UpLo[1],&N,temp.GetMat(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z->GetMat(),&LDZ,ISUPPZ,WORK,&LWORK,IWORK,&LIWORK,&ok);	
	
	// size_t N = A.GetRows(), LDA=A.GetLD();
	// matrix<double> temp(N,N);
	// 
	// int LWORK=3*N-1, ok;
	// 
	// double *W, *WORK;
	// W =new double[N];
	// WORK = new double[LWORK];
	// 
	// 
	// dsyev_(&Jobz[0],&UpLo[1],&N,A.GetMat(),&LDA,W,WORK,&LWORK,&ok);
	if (ok!=0)
		UFs::MyError("Routine ExpM::EigenMethod: An error occured in geting the diagonalization of A");

	//Getting the Dig*Dagger(Z)
	for (size_t i=0; i < N; i++){
		for(size_t j=0; j < N; j++){
			temp(i,j) = std::exp(x*W[i])*(*Z)(j,i);
		}
	}

	temp = (*Z)*temp;
	delete [] W; 
	delete [] ISUPPZ;
	delete [] WORK;
	delete [] IWORK;
	delete Z;
	
	return temp;
}

inline matrix<std::complex<double> > ExpM::EigenMethod(const matrix<std::complex<double> >& A){
	
	if (MOs::IsHermitian(A)==0) UFs::MyError("Routine ExpM::EigenMethod: The input matrix is not Hermitian.");
	size_t N = A.GetRows(), LDA=A.GetLD(), IL, IU, M, LDZ=N;
	double VL, VU, ABSTOL=std::numeric_limits<double>::min(); 
	std::complex<double> *WORK;
	double *RWORK;
	int *IWORK;
	//finding the optimal work sizes
	int LWORK=-1, LRWORK=-1, LIWORK=-1, ok=0;
	WORK = new std::complex<double>[1];
	RWORK = new double[1];
	IWORK = new int[1];
	zheevr_(&Jobz[0],&Range[0],&UpLo[1],&N,NULL,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,NULL,NULL,&LDZ,NULL,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&ok);
	
	//Setting the optimal workspace
	LWORK = (int)real(WORK[0]);
	LRWORK = (int)RWORK[0];
	LIWORK = IWORK[0];
	
	delete WORK;	
	delete RWORK;
	delete IWORK;

	//Running the algorithim 
	matrix<std::complex<double> >* Z;
	static matrix<std::complex<double> > temp;
	temp = A;
	Z = new  matrix<std::complex<double> >(LDZ,LDZ);
	double *W;
	W =new double[LDA];
	int *ISUPPZ;
	ISUPPZ = new int[2*N]; 
	WORK = new std::complex<double>[LWORK];
	RWORK = new double[LRWORK];
	IWORK = new int[LIWORK];
	zheevr_(&Jobz[0],&Range[0],&UpLo[1],&N,temp.GetMat(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z->GetMat(),&LDZ,ISUPPZ,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&ok);	
	
	// size_t N = A.GetRows(), LDA=A.GetLD();
	// matrix<std::complex<double> > temp(A);
	// 
	// int LWORK=2*N-1, ok;
	// 
	// double *W, *RWORK; 
	// std::complex<double> *WORK;
	// W =new double[N];
	// WORK = new std::complex<double>[LWORK];
	// RWORK = new double[3*N-2];
	//
	//zheev_(&Jobz[0],&UpLo[1],&N,A.GetMat(),&LDA,W,WORK,&LWORK,RWORK,&ok);
	// std::cout << ok << std::endl;
	if (ok!=0)
		UFs::MyError("Routine ExpM::EigenMethod: An error occured in geting the diagonalization of A");

	//Getting the Dig*Dagger(Z)
	for (size_t i=0; i < N; i++){
		for(size_t j=0; j < N; j++){
			temp(i,j) = std::exp(W[i])*conj((*Z)(j,i));
		}
	}


	temp = (*Z)*temp;
	delete [] W; 
	delete [] ISUPPZ;
	delete [] WORK;
	delete [] RWORK;
	delete [] IWORK;
	delete Z;
	
	return temp;
}

inline matrix<std::complex<double> > ExpM::EigenMethod(const matrix<std::complex<double> >& A, const std::complex<double>& x, matrix<std::complex<double> >* Z, double *W){	

	bool returnZ(Z!=NULL), returnW(W!=NULL);

	//if (MOs::IsHermitian(A)==0) UFs::MyError("Routine ExpM::EigenMethod: The input matrix is not Hermitian.");
	size_t N = A.GetRows(), LDA=A.GetLD(), IL, IU, M, LDZ=N;
	double VL, VU, ABSTOL=std::numeric_limits<double>::min(); 
	std::complex<double> *WORK;
	double *RWORK;
	int *IWORK;
	//finding the optimal work sizes
	int LWORK=-1, LRWORK=-1, LIWORK=-1, ok=0;
	WORK = new std::complex<double>[1];
	RWORK = new double[1];
	IWORK = new int[1];
	zheevr_(&Jobz[0],&Range[0],&UpLo[1],&N,NULL,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,NULL,NULL,&LDZ,NULL,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&ok);

	//Setting the optimal workspace
	LWORK = (int)real(WORK[0]);
	LRWORK = (int)RWORK[0];
	LIWORK = IWORK[0];

	delete [] WORK;	
	delete [] RWORK;
	delete [] IWORK;

	//Running the algorithim 
	if (!returnZ) Z = new matrix<std::complex<double> >(LDZ,LDZ);
	if (!returnW) W = new double[LDA];
	
	matrix<std::complex<double> > temp(A);
	int *ISUPPZ;
	ISUPPZ = new int[2*N]; 
	WORK = new std::complex<double>[LWORK];	RWORK = new double[LRWORK];
	IWORK = new int[LIWORK];
	
	//std::cout << " A " << temp << std::endl;
		zheevr_(&Jobz[0],&Range[0],&UpLo[1],&N,temp.GetMat(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z->GetMat(),&LDZ,ISUPPZ,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&ok);	
	
	// size_t N = A.GetRows(), LDA=A.GetLD();
	// matrix<std::complex<double> > temp(A);
	// 
	// int LWORK=2*N-1, ok;
	// 
	// double *W, *RWORK; 
	// std::complex<double> *WORK;
	// W =new double[N];
	// WORK = new std::complex<double>[LWORK];
	// RWORK = new double[3*N-2];
	//
	//zheev_(&Jobz[0],&UpLo[1],&N,A.GetMat(),&LDA,W,WORK,&LWORK,RWORK,&ok);
	if (ok!=0)
		UFs::MyError("Routine ExpM::EigenMethod: An error occured in geting the diagonalization of A");

	//Getting the Dig*Dagger(Z)
	for (size_t i=0; i < N; i++){
		for(size_t j=0; j < N; j++){
			temp(i,j) = std::exp(x*W[i])*conj((*Z)(j,i));
		}
	//	std::cout << i << " " << std::exp(x*W[i]) << " ";
	}
	//std::cout << '\n'; 
	temp.SetOutputStyle(Matrix);
	Z->SetOutputStyle(Matrix);
	//std::cout << "Z " <<  *Z << std::endl;
	
	temp = (*Z)*temp;
	
	
	//std::cout << "temp " <<  temp << std::endl;
	
	if(!returnW) delete [] W; 
	if(!returnZ) delete Z; 
	
	delete [] ISUPPZ;
	delete [] WORK;
	delete [] RWORK;
	delete [] IWORK;
	
	return temp;
}
#endif /* MatrixExponential_h */

