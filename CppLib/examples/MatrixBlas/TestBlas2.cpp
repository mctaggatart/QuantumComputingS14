#include <iostream>
#include <complex>
#include <Headers.hpp>

using namespace std;

int main (void)
{
	const size_t rows=3, size2=rows*rows;
	
	complex<double> A[size2], B[size2], C[size2];
	
	A[0*rows+0]=3.1;  A[1*rows+0]=1.3;  A[2*rows+0]=-5.7;	/* matrix A */
	A[0*rows+1]=1.0;  A[1*rows+1]=-6.9; A[2*rows+1]=5.8;	
	A[0*rows+2]=3.4;  A[1*rows+2]=7.2;  A[2*rows+2]=-8.8;
	
	B[0*rows+0]=3.1;  B[1*rows+0]=1.3;  B[2*rows+0]=-5.7;	/* matrix B */
	B[0*rows+1]=1.0;  B[1*rows+1]=-6.9; B[2*rows+1]=5.8;	
	B[0*rows+2]=3.4;  B[1*rows+2]=7.2;  B[2*rows+2]=-8.8;
	
	complex<double> alpha =1.0, beta =0.0;
	
	complex<double>* Apt= A;
	complex<double>* Bpt= B;
	complex<double>* Cpt= C;	
		
	size_t LDA = rows;
	size_t LDB = rows;
	size_t LDC = rows;
	
	zgemm_(&Trans[0], &Trans[0], &rows, &rows, &rows, &alpha, Apt, &LDA, Bpt, &LDB, &beta, Cpt, &LDC);
		for(size_t i = 0; i < rows; ++i)
	{
		for(size_t j = 0; j < rows; ++j)
		{
			cout << C[j*rows+i] << "\t";
		}
		cout << endl;
	}
	
	return 0;

}


