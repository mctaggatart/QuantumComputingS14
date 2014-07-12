/*
Name: TestLAPACK3.cpp
Author: Jay Gambetta

Dependences: Lapack, CMc.h  
Brief Discription: Finds the eigenvectors and eigenvalues of a hermitian matrix
Limitations: None

Version History
	v0: June 10, 2008.
	v1: January 3rd, 2009 - changed it to support my new matrix template -- removed hermitian and symmetric matrix

*/
#include <MatrixOperations.hpp>
// #include <Accelerate/Accelerate.h>
using namespace std;

int main(){

	size_t dim;
	// cout << "enter the size of the matrix "<< endl;
	// 	cin >> dim;
	dim=3;
	matrix<complex<double> > A(dim,dim), D(dim,dim), B(dim,dim);
	A.SetOutputStyle(Matrix);
	D.SetOutputStyle(Matrix);
	B.SetOutputStyle(Matrix);
	

	// lower triagonal
	// cout << "enter the hermitian matrix" << endl;
	// cin >> A;
	// cout << "The matrix you entered is "<< endl;
	A[0*dim+0]=3.1;  A[1*dim+0]=1.0;  A[2*dim+0]=-5.7;	/* matrix A */
	A[0*dim+1]=1.0;  A[1*dim+1]=-6.9; A[2*dim+1]=5.8;	
	A[0*dim+2]=-5.7;  A[1*dim+2]=5.8;  A[2*dim+2]=-8.8;
	cout << A << endl;


	/*SUBROUTINE ZHEEVR( 1-JOBZ, 2-RANGE, 3-UPLO, 4-N, 5-A, 6-LDA, 7-VL, 8-VU, 9-IL, 10-IU,11-ABSTOL, 12-M, 13-W, 14-Z, 15-LDZ, 16-ISUPPZ, 17-WORK, 18-LWORK, 19-RWORK, 20-LRWORK, 21-IWORK, 22-LIWORK, 23-INFO )
	options
	*  Arguments
	*  =========
	*
	*  JOBZ    (input) CHARACTER*1
	*          = 'N':  Compute eigenvalues only;
	*          = 'V':  Compute eigenvalues and eigenvectors.
	*
	*  RANGE   (input) CHARACTER*1
	*          = 'A': all eigenvalues will be found.
	*          = 'V': all eigenvalues in the half-open interval (VL,VU]
	*                 will be found.
	*          = 'I': the IL-th through IU-th eigenvalues will be found.
	********** For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and
	********** ZSTEIN are called
	*
	*  UPLO    (input) CHARACTER*1
	*          = 'U':  Upper triangle of A is stored;
	*          = 'L':  Lower triangle of A is stored.
	*
	*  N       (input) INTEGER
	*          The order of the matrix A.  N >= 0.
	*
	*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
	*          On entry, the Hermitian matrix A.  If UPLO = 'U', the
	*          leading N-by-N upper triangular part of A contains the
	*          upper triangular part of the matrix A.  If UPLO = 'L',
	*          the leading N-by-N lower triangular part of A contains
	*          the lower triangular part of the matrix A.
	*          On exit, the lower triangle (if UPLO='L') or the upper
	*          triangle (if UPLO='U') of A, including the diagonal, is
	*          destroyed.
	*
	*  LDA     (input) INTEGER
	*          The leading dimension of the array A.  LDA >= max(1,N).
	*
	*  VL      (input) DOUBLE PRECISION
	*  VU      (input) DOUBLE PRECISION
	*          If RANGE='V', the lower and upper bounds of the interval to
	*          be searched for eigenvalues. VL < VU.
	*          Not referenced if RANGE = 'A' or 'I'.
	*
	*  IL      (input) INTEGER
	*  IU      (input) INTEGER
	*          If RANGE='I', the indices (in ascending order) of the
	*          smallest and largest eigenvalues to be returned.
	*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
	*          Not referenced if RANGE = 'A' or 'V'.
	*
	*  ABSTOL  (input) DOUBLE PRECISION
	*          The absolute error tolerance for the eigenvalues.
	*          An approximate eigenvalue is accepted as converged
	*          when it is determined to lie in an interval [a,b]
	*          of width less than or equal to
	*
	*                  ABSTOL + EPS *   max( |a|,|b| ) ,
	*
	*          where EPS is the machine precision.  If ABSTOL is less than
	*          or equal to zero, then  EPS*|T|  will be used in its place,
	*          where |T| is the 1-norm of the tridiagonal matrix obtained
	*          by reducing A to tridiagonal form.
	*
	*          See "Computing Small Singular Values of Bidiagonal Matrices
	*          with Guaranteed High Relative Accuracy," by Demmel and
	*          Kahan, LAPACK Working Note #3.
	*
	*          If high relative accuracy is important, set ABSTOL to
	*          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that
	*          eigenvalues are computed to high relative accuracy when
	*          possible in future releases.  The current code does not
	*          make any guarantees about high relative accuracy, but
	*          furutre releases will. See J. Barlow and J. Demmel,
	*          "Computing Accurate Eigensystems of Scaled Diagonally
	*          Dominant Matrices", LAPACK Working Note #7, for a discussion
	*          of which matrices define their eigenvalues to high relative
	*          accuracy.
	*
	*  M       (output) INTEGER
	*          The total number of eigenvalues found.  0 <= M <= N.
	*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
	*
	*  W       (output) DOUBLE PRECISION array, dimension (N)
	*          The first M elements contain the selected eigenvalues in
	*          ascending order.
	*
	*  Z       (output) COMPLEX*16 array, dimension (LDZ, max(1,M))
	*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
	*          contain the orthonormal eigenvectors of the matrix A
	*          corresponding to the selected eigenvalues, with the i-th
	*          column of Z holding the eigenvector associated with W(i).
	*          If JOBZ = 'N', then Z is not referenced.
	*          Note: the user must ensure that at least max(1,M) columns are
	*          supplied in the array Z; if RANGE = 'V', the exact value of M
	*          is not known in advance and an upper bound must be used.
	*
	*  LDZ     (input) INTEGER
	*          The leading dimension of the array Z.  LDZ >= 1, and if
	*          JOBZ = 'V', LDZ >= max(1,N).
	*
	*  ISUPPZ  (output) INTEGER array, dimension ( 2*max(1,M) )
	*          The support of the eigenvectors in Z, i.e., the indices
	*          indicating the nonzero elements in Z. The i-th eigenvector
	*          is nonzero only in elements ISUPPZ( 2*i-1 ) through
	*          ISUPPZ( 2*i ).
	********** Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
	*
	*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
	*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
	*
	*  LWORK   (input) INTEGER
	*          The length of the array WORK.  LWORK >= max(1,2*N).
	*          For optimal efficiency, LWORK >= (NB+1)*N,
	*          where NB is the max of the blocksize for ZHETRD and for
	*          ZUNMTR as returned by ILAENV.
	*
	*          If LWORK = -1, then a workspace query is assumed; the routine
	*          only calculates the optimal sizes of the WORK, RWORK and
	*          IWORK arrays, returns these values as the first entries of
	*          the WORK, RWORK and IWORK arrays, and no error message
	*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
	*
	*  RWORK   (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LRWORK))
	*          On exit, if INFO = 0, RWORK(1) returns the optimal
	*          (and minimal) LRWORK.
	*
	* LRWORK   (input) INTEGER
	*          The length of the array RWORK.  LRWORK >= max(1,24*N).
	*
	*          If LRWORK = -1, then a workspace query is assumed; the
	*          routine only calculates the optimal sizes of the WORK, RWORK
	*          and IWORK arrays, returns these values as the first entries
	*          of the WORK, RWORK and IWORK arrays, and no error message
	*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
	*
	*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
	*          On exit, if INFO = 0, IWORK(1) returns the optimal
	*          (and minimal) LIWORK.
	*
	* LIWORK   (input) INTEGER
	*          The dimension of the array IWORK.  LIWORK >= max(1,10*N).
	*
	*          If LIWORK = -1, then a workspace query is assumed; the
	*          routine only calculates the optimal sizes of the WORK, RWORK
	*          and IWORK arrays, returns these values as the first entries
	*          of the WORK, RWORK and IWORK arrays, and no error message
	*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
	*
	*  INFO    (output) INTEGER
	*          = 0:  successful exit
	*          < 0:  if INFO = -i, the i-th argument had an illegal value
	*          > 0:  Internal error
	*/
	
	size_t N = dim, LDA=dim, IL, IU, M, LDZ=dim; 
	double VL, VU, ABSTOL=numeric_limits<double>::min(); 
	complex<double> *WORK;
	double *RWORK;
	int *IWORK;
	//finding the optimal work sizes
	int LWORK=-1, LRWORK=-1, LIWORK=-1, ok=0;
	WORK = new complex<double>[1];
	RWORK = new double[1];
	IWORK = new int[1];
	zheevr_(&Jobz[0],&Range[0],&UpLo[1],&N,NULL,&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,NULL,NULL,&LDZ,NULL,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&ok);
	
	//Setting the optimal workspace
	LWORK = (int)real(WORK[0]);
	LRWORK = (int)RWORK[0];
	LIWORK = IWORK[0];
	
	//Running the algorithim 
	matrix<complex<double> > Z(LDZ,LDZ);
	Z=A;
	Z.SetOutputStyle(Matrix);
	double *W;
	W =new double[LDA];
	int *ISUPPZ;
	ISUPPZ = new int[2*N]; 
	WORK = new complex<double>[LWORK];
	RWORK = new double[LRWORK];
	IWORK = new int[LIWORK];
	zheevr_(&Jobz[0],&Range[0],&UpLo[1],&N,A.GetMat(),&LDA,&VL,&VU,&IL,&IU,&ABSTOL,&M,W,Z.GetMat(),&LDZ,ISUPPZ,WORK,&LWORK,RWORK,&LRWORK,IWORK,&LIWORK,&ok);	
	
	//outputing the optimal values
	cout << LWORK <<endl;
	cout << LRWORK <<endl;
	cout << LIWORK <<endl;
	cout << M << endl;
	// 
	if (ok==0){
		cout << "The eigenvalues are " << endl;
		for (size_t i=0; i<N; i++){
			D(i,i)=W[i];
	   	}
		cout << D << endl;
		cout << "The eigenvector matrix is " << endl;
		cout << Z << endl;
		cout << A << endl;
	}
	else cout << "An error occured" << endl;


	cout << "The orginal matrix is UDUdg " << endl;
	B = Z*D*MOs::Dagger(Z);
	cout << B << endl;
	cout << "The diagonal matrix is UDUdg " << endl;
	B=(MOs::Dagger(Z)*B*Z);
	cout << B << endl;
	
	for(size_t j=0; j<2*N; j++){
		cout << IWORK[j] << endl;
	}

	delete [] W; 
	delete [] WORK;
	delete [] RWORK;
	delete [] IWORK;
	delete [] ISUPPZ;
	
	return 0;
}
