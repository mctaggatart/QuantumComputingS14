/*
Name: Headers.hpp
Author: Jay Gambetta

Dependences: BLAS
Brief Discription: 

Limitations: Sparse Matrices are still not included

Version History
	v0: Feb 2nd, 2007. Started the file
*/

#ifndef Headers_h
#define Headers_h
// #include <Accelerate/Accelerate.h>

const char Trans[] = {'N','T','C'};
/*  Trans (input) CHARACTER*1.
		On entry, TRANSA specifies the form of op( A ) to be used in the matrix multiplication as follows:
 			= 'N' no transpose;
			= 'T' transpose of A;
 			= 'C' hermitian conjugate of A.
*/
const char UpLo[] = {'U','L'};
/*  UpLo    (input) CHARACTER*1
 			= 'U':  Upper triangle of A is stored;
			= 'L':  Lower triangle of A is stored.
*/
const char Jobz[] = {'V','N'};
/*	Jobz    (input) CHARACTER*1
 			= 'N':  Compute eigenvalues only;
			= 'V':  Compute eigenvalues and eigenvectors.
*/
const char Range[] = {'A','V','I'};
/*  Range   (input) CHARACTER*1
				= 'A': all eigenvalues will be found.
				= 'V': all eigenvalues in the half-open interval (VL,VU] will be found.
				= 'I': the IL-th through IU-th eigenvalues will be found.
*/

#ifdef __cplusplus
extern "C" {
#endif

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

void dgemm_(const char* TransA, const char* TransB, const size_t* M, const size_t* N, const size_t* K, const double* alpha, const double* A, const size_t* lda, const double* B, const size_t* lba, const double* beta, double* C, size_t* ldc);
void zgemm_(const char* TransA, const char* TransB, const size_t* M, const size_t* N, const size_t* K,
const std::complex<double>* alpha, const std::complex<double>* A, const size_t* lda, const std::complex<double>* B, const size_t* ldb, const std::complex<double>* beta, std::complex<double>* C, size_t* ldc);
void zher2k_(const char* UpLo, const char* Trans, const size_t* N, const size_t* K, const std::complex<double>* alpha, const std::complex<double>* A, const size_t* lda, const std::complex<double>* B, const size_t* ldb, const double* beta, std::complex<double>* C, const size_t* ldc);

/*
 * ===========================================================================
 * Prototypes for lapack
 * ===========================================================================
 */

int dgeev_(const char* jobvl, const char* jobvr, const size_t* n, double* a, size_t* lda, double* wr, double* wi, double* vl, size_t* ldvl, double* vr, size_t* ldvr, double* work, int* lwork, int* info);
int zgeev_(const char* jobvl, const char* jobvr, const size_t* n, std::complex<double>* a, size_t* lda, std::complex<double>* w, std::complex<double>* vl, size_t* ldvl, std::complex<double>* vr, size_t* ldvr, std::complex<double>* work, int* lwork, double* rwork, int* info);
int dsyev_(const char* jobz, const char* uplo, const size_t* n, double* a, size_t* lda, double* w, double* work, int* lwork, int* info);
int zheev_(const char* jobz, const char* uplo, const size_t* n, std::complex<double>* a, size_t* lda, double* w, std::complex<double>* work, int* lwork, double* rwork, int* info);	
int dsyevr_(const char* jobz, const char* range, const char* uplo, const size_t* n, double* a, size_t* lda, double* vl, double* vu, size_t* il, size_t* iu, double* abstol, size_t* m, double* w, double* z, size_t* ldz, int* isuppz, double* work, int* lwork, int* iwork, int* liwork, int* info);
int zheevr_(const char* jobz, const char* range, const char* uplo, const size_t* n, std::complex<double>* a, size_t* lda, double* vl, double* vu, size_t* il, size_t* iu, double* abstol, size_t* m, double* w, std::complex<double>* z, size_t* ldz, int* isuppz, std::complex<double>* work, int* lwork, double* rwork, int* lrwork, int* iwork, int* liwork, int* info);


#ifdef __cplusplus
}
#endif

#endif /* Headers_h */

