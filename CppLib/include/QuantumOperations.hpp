/*
Copyright 2010 Jay M. Gambetta

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

	http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/
/*

Dependences: Matrix.hpp
Brief Discription: 	A header Containing My quantum operators, states, superoperators etc
Limitations: None

Version History
	v0: Dec 10, 2007
	v1:	June 27, 2008. Changed to support my new storage for hermitian matrices
	v2: January 3rd, 2009 - changed it to support my new matrix template -- removed hermitian and symmetric matrix
	v3: May 20, 2009 - Fixed some bugs, added OutputMatrix
	v4: September 12, 2009 -  Move OutputMatrix to UsefullFunctions

*/
#ifndef QuantumOperations_h
#define QuantumOperations_h
#include "UsefullFunctions.hpp"
#include "MatrixOperations.hpp"
#include "MatrixExponential.hpp"

namespace QOs { //All my quantum operators
	
	
	
	//States
//	template <class T>
	template <class T>
	void FockState(matrix<T>& rho, const size_t&  m);
	//void CoherentState(matrix<std::complex<double> >& rho, const std::complex<double>& alpha);
	template <class T>
	void ThermalState(matrix<T>& rho, const double& nth);
	//void SqueezedState(matrix<std::complex<double> >& rho, const std::complex<double>& alpha, const std::complex<double>& zeta);
	void QubitState(matrix<std::complex<double> >& rho, const double& x, const double& y, const double& z);
	void TwoQubitState(matrix<std::complex<double> >& Rho, const std::vector<double>& pauli);
	
	//State representations 
	//void Qfunction(const matrix<std::complex<double> >& rho, matrix<std::complex<double> >& Qxy);
	//void Wignerfunction(const matrix<std::complex<double> >& rho, matrix<std::complex<double> >& Wxy);  

	
	//Expectation Values
	template <class T>
	std::complex<double> Expectation(const matrix<T>& A, const matrix<std::complex<double> >& rho);
	
	//Dense Superoperators gets only the lower triangle part
	//-iH rho + i rho H  
	matrix<std::complex<double> > SuperHami(const matrix<std::complex<double> >& H, const matrix<std::complex<double> >& rho);
	//2a rho ad - adg a rho - rho a adg
	matrix<std::complex<double> > SuperDamp(const matrix<std::complex<double> >& a, const matrix<std::complex<double> >& rho);
	//2 a rho ad
	matrix<std::complex<double> > SuperJump(const matrix<std::complex<double> >& a, const matrix<std::complex<double> >& rho);
	//a rho + rho adg
	matrix<std::complex<double> > SuperAnti(const matrix<std::complex<double> >& A, const matrix<std::complex<double> >& rho);
	// //a rho + rho adg - <a+adg> rho
	matrix<std::complex<double> > SuperHomo(const matrix<std::complex<double> >& A, const matrix<std::complex<double> >& rho);
	// exp(alpha ad - alpha* a)
	matrix<std::complex<double> > Displace(const size_t dim, const std::complex<double> alpha);
	
	//quantities
	double Purity(const matrix<std::complex<double> >& rho);
	void EigenEnergies(const matrix<std::complex<double> >& A, const size_t num, double* vals);
	void EigenSystem(const matrix<std::complex<double> >& A, const size_t num, double* vals, matrix<std::complex<double> >& Z);
	void MatrixElements(const matrix<std::complex<double> >& A, const size_t num, double* vals, matrix<std::complex<double> >& B);
	
	//2 qubit quanties
	double Concurrence(const matrix<std::complex<double> >& rho);


}
//------------Member funtions------------------
inline void QOs::QubitState(matrix<std::complex<double> >& rho, const double& x, const double& y, const double& z){
	//Initial Qubit state
	// z is one puts the qubit in the state 1,0 and z equal -1 puts the qubit in the state 0,1 
	const double TINY = std::numeric_limits<double>::epsilon();
	if ((std::pow(x,2)+std::pow(y,2)+std::pow(z,2) ) > (1+TINY))
		UFs::MyError("Routine QOs::QubitState: The state you enterend is not on the bloch sphere, its purity is greater then 1");
	if (rho.GetRows() != 2)
		UFs::MyError("Routine QOs::QubitState: The dimensions of the qubit space is not 2");
	rho(0,0) =0.5*(1.0+z);
	rho(1,0) =0.5*std::complex<double>(x,y);
	rho(0,1) =std::conj(rho(1,0));
	rho(1,1) =0.5*(1.0-z);
}
inline void QOs::TwoQubitState(matrix<std::complex<double> >& Rho, const std::vector<double>&pauli){
	if (Rho.GetRows() != 4)
		UFs::MyError("Routine QOs::TwoQubitState: The dimensions of the 2 qubit space is not 4"); 
	std::complex<double> I(0,1);
	
	//sii, sxi, syi, szi, six, sxx, syx, szx, siy, sxy, syy, szy, siz, sxz, syz, szz
	// 0    1    2    3    4    5    6    7    8    9    10   11   12   13   14   15
	 
	Rho(0,0)=0.25*(1 + pauli[12] + pauli[3] + pauli[15]);
	Rho(1,0)=0.25*(pauli[4] + I*pauli[8] + pauli[7] + I*pauli[11]);
	Rho(2,0)=0.25*(pauli[1] + pauli[13] + I*pauli[2] + I*pauli[14]);
	Rho(3,0)=0.25*(pauli[5] + I*pauli[9] + I*pauli[6] - pauli[10]);
	Rho(1,1)=0.25*(1 - pauli[12] + pauli[3] - pauli[15]);
	Rho(2,1)=0.25*(pauli[5] - I*pauli[9] + I*pauli[6] + pauli[10]);
	Rho(3,1)=0.25*(pauli[1] - pauli[13] + I*pauli[2] - I*pauli[14]);
	Rho(2,2)=0.25*(1 + pauli[12] - pauli[3] - pauli[15]);
	Rho(3,2)=0.25*(pauli[4] + I*pauli[8] - pauli[7] - I*pauli[11]);
	Rho(3,3)=0.25*(1 - pauli[12] - pauli[3] + pauli[15]);
	
	for (size_t i=0; i< 4; i++)
		for(size_t j=i+1; j<4; j++)
			Rho(i,j)=std::conj(Rho(j,i));
}
template <class T>
inline void QOs::FockState(matrix<T>& rho, const size_t& n){
	//Initial fock state, rho=|n><n| 
	if (n > rho.GetRows())
		UFs::MyError("Routine QOs::FockState: the fock state chosen exceeds the hilbert space dimensions");
	MOs::Null(rho);
	rho(n,n)=T(1.0);
}
// inline void QOs::CoherentState(matrix<std::complex<double> >& rho, const std::complex<double>& alpha){
// 	//Coherent state rho=exp(-|alpha|^2) sum_(n,m) alpha^n (alpha^*)^m/sqrt(n!m!) ket{n}bra{m}
// 	
// 	// exp (log(rho_nm)) = exp(-|a|^2 + n*log(alpha) + m*log(alpha)- 0.5*SFs::LogFactorial(n) -0.5*SFs::LogFactorial(m));
// 	std::complex<double> norm = exp(-alpha*conj(alpha));
// 	for (size_t n = 0; n < rho.GetRows(); ++n){
// 		for(size_t m = 0; m < rho.GetColumns(); ++m){
			// if n < 10 && m <10
			//rho(n,m)=norm*pow(alpha)
			//else
// 			rho(n,m)=exp(-conj(alpha)*alpha + n*log(alpha) + m*log(conj(alpha))- 0.5*SFs::LogFactorial(n) -0.5*SFs::LogFactorial(m));
// 		}
// 	}
// }
template <class T>
inline void QOs::ThermalState(matrix<T>& rho, const double& nth){
	//Thermal state rho=sum( p_n |n><n| ), where p_n = 1/(1+nth) (nth/(1+nth))^n
	MOs::Null(rho);
	double pn;
	for (size_t n=0; n<rho.GetRows(); ++n){
		pn = 1/(1+nth)*std::pow( nth/(1+nth), (int)n);
		rho(n,n)=T(pn);
	}
}
template <class T>
inline std::complex<double> QOs::Expectation(const matrix<T>& A, const matrix<std::complex<double> >& rho ){
	//Calculates the expectation value of the operator A
   	size_t dim = rho.GetRows();
	if (dim != A.GetRows() || dim != A.GetColumns())
			UFs::MyError("Routine QOs::Expectation: matrix dimensions and rho dont agree");
 	std::complex<double> temp=0;
	for(size_t i=0; i< dim; i++){
		for(size_t p=0; p<dim; p++){
			temp += A(i,p)*rho(p,i);
		}
	}
	return temp;
}
inline matrix<std::complex<double> > QOs::SuperHami(const matrix<std::complex<double> >& H, const matrix<std::complex<double> >& rho){
	//Gets the superoperator rhonnew = -iH rho + i rho H 
	//cblas_zher2k(CblasXMajor,uplo,op,N,K,alpha,A,LDA,B,LDB,beta,C,LDC)
	// C-> alpha*A*B^dagger + alpha^* B A^dagger +beta C
	size_t dim = H.GetRows();		
	matrix<std::complex<double> > C(dim,dim);
	std::complex<double> alpha(0.0,-1.0); 
	double beta = 0.0;
	zher2k_(&UpLo[1], &Trans[0], &dim, &dim, &alpha, H.GetMat(), &dim, rho.GetMat(), &dim, &beta, C.GetMat(), &dim);
	// cblas_zher2k(CblasColMajor, CblasLower, CblasNoTrans, dim, H.GetColumns(), &alpha, H.GetMat(), H.GetLD(), rho.GetMat(), rho.GetLD(), 0.0, C.GetMat(), C.GetLD());

	return C;
}
inline matrix<std::complex<double> > QOs::SuperDamp(const matrix<std::complex<double> >& a, const matrix<std::complex<double> >& rho){
	//Gets the superoperator rhonnew = 2A*rho*A^dag -A^dag*A*rho -rho*A^dag*A 
	//cblas_zgemm(CblasXMajor,op,op,N,M,K,alpha,A,LDA,B,LDB,beta,C,LDC)
	// C-> alpha*op(A)*op(B) +beta C
	//cblas_zher2k(CblasXMajor,uplo,op,N,K,alpha,A,LDA,B,LDB,beta,C,LDC)
	// C-> alpha*A*B^dagger + alpha^* B A^dagger +beta C
	size_t dim = rho.GetRows();
	if (dim != a.GetRows() || dim != a.GetColumns())
			UFs::MyError("Routine QOs::SuperDamp: matrix dimensions and rho dont agree");
	matrix<std::complex<double> > ap(dim,dim), C(dim,dim);
	std::complex<double> alpha(1.0,0.0), beta(0.0,0.0); 
	double gamma = 0.0, delta = -1.0;
	zgemm_(&Trans[0], &Trans[0], &dim, &dim, &dim, &alpha, a.GetMat(), &dim, rho.GetMat(), &dim, &beta, ap.GetMat(), &dim);
	zher2k_(&UpLo[1], &Trans[2], &dim, &dim, &alpha, a.GetMat(), &dim, ap.GetMat(), &dim, &gamma, C.GetMat(), &dim);
	zher2k_(&UpLo[1], &Trans[0], &dim, &dim, &alpha, a.GetMat(), &dim, ap.GetMat(), &dim, &delta, C.GetMat(), &dim);
	// cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &alpha, a.GetMat(), dim, rho.GetMat(), dim, &beta, ap.GetMat(), dim);
	// 	cblas_zher2k(CblasColMajor, CblasLower, CblasConjTrans, dim, dim, &alpha, a.GetMat(), dim, ap.GetMat(), dim, 0.0, C.GetMat(), dim);
	// 	cblas_zher2k(CblasColMajor, CblasLower, CblasNoTrans, dim, dim, &alpha, a.GetMat(), dim, ap.GetMat(), dim, -1.0, C.GetMat(), dim);
		
	return C;
}
inline matrix<std::complex<double> > QOs::SuperJump(const matrix<std::complex<double> >& a, const matrix<std::complex<double> >& rho){
	//Gets the superoperator rhonnew = 2A*rho*A^dag 
	//cblas_zgemm(CblasXMajor,op,op,N,M,K,alpha,A,LDA,B,LDB,beta,C,LDC)
	// C-> alpha*op(A)*op(B) +beta C
	//cblas_zher2k(CblasXMajor,uplo,op,N,K,alpha,A,LDA,B,LDB,beta,C,LDC)
	// C-> alpha*A*B^dagger + alpha^* B A^dagger +beta C
	size_t dim = rho.GetRows();
	if (dim != a.GetRows() || dim != a.GetColumns())
			UFs::MyError("Routine QOs::SuperDamp: matrix dimensions and rho dont agree");
	matrix<std::complex<double> > ap(dim,dim), C(dim,dim);;
	std::complex<double> alpha(1.0,0.0), beta(0.0,0.0);
	double gamma = 0.0; 
	zgemm_(&Trans[0], &Trans[0], &dim, &dim, &dim, &alpha, a.GetMat(), &dim, rho.GetMat(), &dim, &beta, ap.GetMat(), &dim);
	zher2k_(&UpLo[1], &Trans[0], &dim, &dim, &alpha, a.GetMat(), &dim, ap.GetMat(), &dim, &gamma, C.GetMat(), &dim);
	// cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, dim, dim, dim, &alpha, a.GetMat(), dim, rho.GetMat(), dim, &beta, ap.GetMat(), dim);
	// cblas_zher2k(CblasColMajor, CblasLower, CblasNoTrans, dim, dim, &alpha, a.GetMat(), dim, ap.GetMat(), dim, 0.0, C.GetMat(), dim);
	
	return C;
}
inline matrix<std::complex<double> > QOs::SuperAnti(const matrix<std::complex<double> >& A, const matrix<std::complex<double> >& rho){
	//Gets the superoperator rhonnew = A rho + rho A 
	//cblas_zher2k(CblasXMajor,uplo,op,N,K,alpha,A,LDA,B,LDB,beta,C,LDC)
	// C-> alpha*A*B^dagger + alpha^* B A^dagger +beta C
	size_t dim = A.GetRows();		
	matrix<std::complex<double> > C(dim,dim);
	std::complex<double> alpha(1.0,0.0); 
	double beta = 0.0;
	zher2k_(&UpLo[1], &Trans[0], &dim, &dim, &alpha, A.GetMat(), &dim, rho.GetMat(), &dim, &beta, C.GetMat(), &dim);
	// cblas_zher2k(CblasColMajor, CblasLower, CblasNoTrans, dim, A.GetColumns(), &alpha, A.GetMat(), A.GetLD(), rho.GetMat(), rho.GetLD(), 0.0, C.GetMat(), C.GetLD());
	return C;
}
inline matrix<std::complex<double> > QOs::SuperHomo(const matrix<std::complex<double> >& A, const matrix<std::complex<double> >& rho){
	//Gets the superoperator rhonnew = a rho + rho adg - <a+adg> rho
	//cblas_zher2k(CblasXMajor,uplo,op,N,K,alpha,A,LDA,B,LDB,beta,C,LDC)
	// C-> alpha*A*B^dagger + alpha^* B A^dagger +beta C
	size_t dim = A.GetRows();		
	matrix<std::complex<double> > C(rho,'L');
	// std::cout << C << std::endl;
	std::complex<double> alpha(1.0,0.0), ave = -QOs::Expectation(A, rho);
	ave+=std::conj(ave);
	double beta = std::real(ave);
	// std::cout << ave << std::endl;
	// std::cout << beta << std::endl;
	zher2k_(&UpLo[1], &Trans[0], &dim, &dim, &alpha, A.GetMat(), &dim, rho.GetMat(), &dim, &beta, C.GetMat(), &dim);
	// cblas_zher2k(CblasColMajor, CblasLower, CblasNoTrans, dim, A.GetColumns(), &alpha, A.GetMat(), A.GetLD(), rho.GetMat(), rho.GetLD(), beta, C.GetMat(), C.GetLD());
	return C;
}
inline matrix<std::complex<double> > QOs::Displace(const size_t dim, const std::complex<double> alpha){
	matrix<std::complex<double> > tmp(dim, dim);
	MOs::Destroy(tmp);
	tmp = alpha*MOs::Dagger(tmp) + conj(alpha)*tmp;
	return tmp=ExpM::EigenMethod(tmp);
}
inline double QOs::Concurrence(const matrix<std::complex<double> >& rho){
	if (4 != rho.GetRows())
			UFs::MyError("Routine QOs::Concurrence: Rho is not a 2 qubit state");
			
		
	matrix<std::complex<double> > Syy(4,4), rho_con(4,4);
	MOs::Null(Syy);
	Syy(0,3)=-1.0;
	Syy(1,2)=1.0;
	Syy(2,1)=1.0;
	Syy(3,0)=-1.0;
		
	Syy = rho*Syy*MOs::Conjugate(rho)*Syy;

	//  SUBROUTINE ZGEEV( JOBVL, JOBVR, N, A, LDA, W, VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
	//options
	size_t N = 4, LDA=4, LDVL=4, LDVR=4;
	int ok, LWORK=2*4;
	std::complex<double>* VLpt =Syy.GetMat(); 
	std::complex<double>* VRpt=Syy.GetMat();
	std::complex<double> *W, *WORK;
	double *RWORK;

	W =new std::complex<double>[LDA];
	WORK = new std::complex<double>[LWORK];
	RWORK = new double[LWORK]; 
	
	zgeev_(&Jobz[1],&Jobz[1], &N, Syy.GetMat(), &LDA, W, VLpt, &LDVL, VRpt, &LDVR, WORK, &LWORK, RWORK, &ok);
	// zgeev_(&JOBVL,&JOBVR, &N, (__CLPK_doublecomplex *)Syy.GetMat(), &LDA, (__CLPK_doublecomplex *)W, (__CLPK_doublecomplex *)VLpt, &LDVL, (__CLPK_doublecomplex *)VRpt, &LDVR, (__CLPK_doublecomplex *)WORK, &LWORK, RWORK, &ok);
	// 
	// std::cout << W[0] <<  std::endl;
	// std::cout << W[1] <<  std::endl;
	// std::cout << W[2] <<  std::endl;
	// std::cout << W[3] <<  std::endl;
	
	double temp = 0.0;
	for (size_t i=0; i < 4; i++){
		if( std::real(W[i])> temp )
			temp=std::real(W[i]);
	}
	
	temp = 2.0*std::sqrt(temp)-std::real(std::sqrt(W[0])+std::sqrt(W[1])+std::sqrt(W[2])+std::sqrt(W[3]));
	// std::cout << imag(sqrt(W[3])-sqrt(W[2])-sqrt(W[1])-sqrt(W[0])) << std::endl;
	delete [] W;
	delete [] WORK;
	delete [] RWORK;

	return USs::Max<double>(0,temp);
	
}
inline double QOs::Purity(const matrix<std::complex<double> >& rho){
	size_t dim = rho.GetRows();
 	std::complex<double> temp=0;
	for(size_t i=0; i< dim; i++){
		for(size_t p=0; p<dim; p++){
			temp += rho(i,p)*rho(p,i);
		}
	}
	return std::real(temp);
}
inline void QOs::EigenEnergies(const matrix<std::complex<double> >& A, const size_t num, double* vals){
	
	if (MOs::IsHermitian(A)==0) UFs::MyError("Routine QOs::EigenEnergies: The input matrix is not Hermitian.");
	
	
	size_t N = A.GetRows(), LDA=A.GetLD(), IL=1, IU=num, M=num, LDZ=N;
	double VL, VU, ABSTOL=std::numeric_limits<double>::min(); 
	std::complex<double> *WORK;
	double *RWORK;
	int *IWORK;
	//finding the optimal work sizes
	int LWORK=-1, LRWORK=-1, LIWORK=-1, ok=0;
	WORK = new std::complex<double>[1];
	RWORK = new double[1];
	IWORK = new int[1];
	zheevr_(&Jobz[1],&Range[2],&UpLo[1],&N,NULL,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,NULL,NULL,&LDZ,NULL,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&ok);
	//std::cout << "here" << M;
	//Setting the optimal workspace
	LWORK = (int)real(WORK[0]);
	LRWORK = (int)RWORK[0];
	LIWORK = IWORK[0];
	
	//Running the algorithim 
	matrix<std::complex<double> > Z(LDZ,M), temp(A);
	double *W;
	W =new double[LDA];
	int *ISUPPZ;
	ISUPPZ = new int[2*M]; 
	delete [] WORK;
	delete [] RWORK;
	delete [] IWORK;
	WORK = new std::complex<double>[LWORK];
	RWORK = new double[LRWORK];
	IWORK = new int[LIWORK];
	zheevr_(&Jobz[1],&Range[2],&UpLo[1],&N,temp.GetMat(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z.GetMat(),&LDZ,ISUPPZ,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&ok);	
	if (ok==0){
		for (size_t i=0; i<num; i++){
			vals[i] = W[i];
	   	}
	}
	
	else
		UFs::MyError("Routine QOs::EigenEnergies: An error occured in geting the diagonalization of A");



	delete [] W; 
	delete [] ISUPPZ;
	delete [] WORK;
	delete [] RWORK;
	delete [] IWORK;

}
inline void QOs::EigenSystem(const matrix<std::complex<double> >& A, const size_t num, double* vals, matrix<std::complex<double> >& Z){
// returns the eigen values and the unitary with the eigenvectors on the diagonal

	
	if (MOs::IsHermitian(A)==0) UFs::MyError("Routine QOs::EigenSystem: The input matrix is not Hermitian.");
	
	size_t N = A.GetRows(), LDA=A.GetLD(), IL=1, IU=num, M=num, LDZ=N;
	if (Z.GetRows()!=1 && Z.GetColumns()!=M) UFs::MyError("Routine QOs::EigenSystem: The input matrix is not of the zize dim X num .");
	double VL, VU, ABSTOL=std::numeric_limits<double>::min(); 
	std::complex<double> *WORK;
	double *RWORK;
	int *IWORK;
	//finding the optimal work sizes
	int LWORK=-1, LRWORK=-1, LIWORK=-1, ok=0;
	WORK = new std::complex<double>[1];
	RWORK = new double[1];
	IWORK = new int[1];
	zheevr_(&Jobz[0],&Range[2],&UpLo[1],&N,NULL,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,NULL,NULL,&LDZ,NULL,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&ok);
	//std::cout << "here" << M;
	//Setting the optimal workspace
	LWORK = (int)real(WORK[0]);
	LRWORK = (int)RWORK[0];
	LIWORK = IWORK[0];
	
	//Running the algorithim 
	matrix<std::complex<double> > temp(A);
	double *W;
	W =new double[LDA];
	int *ISUPPZ;
	ISUPPZ = new int[2*M]; 
	delete [] WORK;
	delete [] RWORK;
	delete [] IWORK;
	WORK = new std::complex<double>[LWORK];
	RWORK = new double[LRWORK];
	IWORK = new int[LIWORK];
	zheevr_(&Jobz[0],&Range[2],&UpLo[1],&N,temp.GetMat(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z.GetMat(),&LDZ,ISUPPZ,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&ok);	
	if (ok==0){
		for (size_t i=0; i<num; i++){
			vals[i] = W[i];
	   	}
	}
	else
		UFs::MyError("Routine QOs::EigenSystem: An error occured in geting the diagonalization of A");

	
	//std::cout << Z << std::endl;

	delete [] W; 
	delete [] ISUPPZ;
	delete [] WORK;
	delete [] RWORK;
	delete [] IWORK;

}
inline void QOs::MatrixElements(const matrix<std::complex<double> >& A, const size_t num, double* vals, matrix<std::complex<double> >& B){
	
	if (MOs::IsHermitian(A)==0) UFs::MyError("Routine QOs::MatrixElements: The input matrix is not Hermitian.");
	
	
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
	
	//Running the algorithim 
	matrix<std::complex<double> > Z(LDZ,LDZ), temp(A);
	double *W;
	W =new double[LDA];
	int *ISUPPZ;
	ISUPPZ = new int[2*N]; 
	delete [] WORK;
	delete [] RWORK;
	delete [] IWORK;
	WORK = new std::complex<double>[LWORK];
	RWORK = new double[LRWORK];
	IWORK = new int[LIWORK];
	zheevr_(&Jobz[0],&Range[0],&UpLo[1],&N,temp.GetMat(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z.GetMat(),&LDZ,ISUPPZ,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&ok);	

	if (ok==0){
		for (size_t i=0; i<num; i++){
			vals[i] = W[i];
	   	}
	   	
	   	B = MOs::Dagger(Z)*B*Z;
	   	//B.SetOutputStyle(Matrix);
	   	//std::cout << B << std::endl;

	}
	
	else
		UFs::MyError("Routine QOs::EigenEnergies: An error occured in geting the diagonalization of A");


	delete [] W; 
	delete [] ISUPPZ;
	delete [] WORK;
	delete [] RWORK;
	delete [] IWORK;

}
#endif /* QuantumOperations_h */

