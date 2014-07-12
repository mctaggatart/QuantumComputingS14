/*
Name: TestSuperDamp.cpp
Author: Jay Gambetta

Dependences: MatrixOperations.hpp, QuantumOperations.hpp
Brief Discription: This program benchmarks different wast of getting the SuperOperator for -i[H,rho].
Limitations: None

Version History
	v0: January 4th, 2009

run in terminal
	i="0"; while [ $i -lt 8 ]; do i=$[$i+1]; echo $i; ./TestSuperDamp.e $i TestSuperDamp.dat; done

*/
#include <MatrixOperations.hpp>
#include <QuantumOperations.hpp>
using namespace std;

inline matrix<complex<double> > SuperDamp(const matrix<complex<double> >& a, const matrix<complex<double> >& rho){
	size_t dim= rho.GetRows(), rows1 = a.GetRows(), cols1 = a.GetColumns();
	if (dim != rows1 || dim != cols1)
		UFs::MyError("dimensions of operators and rho dont agree");

	complex<double> tempc;
	matrix<complex<double> > temprho(dim,dim);
	for (unsigned int i=0; i < dim; i++){
		for(unsigned int p=0; p < dim; p++){
			for(unsigned int q=0; q < dim; q++){
				//diagonal part
				temprho[i*dim+i] += 2.0*a(i,p)*rho(p,q)*conj(a(i,q)) - conj(a(p,i))*a(p,q)*rho(q,i) - rho(i,p)*conj(a(q,p))*a(q,i);
				//cout << imag(tempdia) << endl;
				//lowerTriangle
				for(unsigned int j=0; j < i; j++){
					temprho(i,j)+=2.0*a(i,p)*rho(p,q)*conj(a(j,q)) - conj(a(p,i))*a(p,q)*rho(q,j) - rho(i,p)*conj(a(q,p))*a(q,j);
				}
			}
		}
	}
	return temprho;
}


int main (int argc, char const *argv[]){
	
	cout << "Running program " << argv[0] << endl;
	//cout << CLOCKS_PER_SEC << endl;
	double to_ms =1000.0/CLOCKS_PER_SEC;
	cout << to_ms << endl;
	size_t dim=int(pow(2.0,atof(argv[1])));
	string const outfile = argv[2];
	ofstream dataout;
	UFs::AppendFile(outfile,dataout, 16);
	
	/* initialize random seed: */
	srand ( time(NULL) );
	
	matrix<complex<double> > A(dim,dim), rho(dim, dim); 
	
	for(size_t i = 0; i < dim; i++){
		rho(i,i)=(rand()%100)/100.0;
		for (size_t j=0; j < i; j++){
			rho(i,j)=complex<double>((rand()%100)/100.0,(rand()%100)/100.0);
			rho(j,i) = conj(rho(i,j));
		}
	}
	
	for(size_t i = 0; i < dim; i++){
		A(i,i)=(rand()%100)/100.0;
		for (size_t j=0; j < i; j++){
			A(i,j)=complex<double>((rand()%100)/100.0,(rand()%100)/100.0);
			A(j,i) = conj(A(i,j));
		}
	}
	
	matrix<complex<double> > c(dim,dim);
	c.SetOutputStyle(Matrix);
	//A.SetOutputStyle(Matrix);
	//rho.SetOutputStyle(Matrix);
	dataout << argv[1];
	
	cout << "here" << endl;
	//Method 1
	int clo = clock();
	c=SuperDamp(A,rho);
	dataout << '\t' << to_ms*(clock() - clo);
	//cout << A << endl;
	//cout << rho << endl;
	// cout << c << endl;

	//Method 2
	clo = clock();
	matrix<complex<double> > ap(dim,dim);
	complex<double> alpha(1.0,0.0),beta(0.0,0.0); 
	size_t M=A.GetRows();
	size_t K=A.GetColumns();
	size_t N=rho.GetColumns();
	double gamma = 0.0, delta = -1.0;
	
	size_t LDA = N;
	size_t LDB = N;
	size_t LDC = N;
	zgemm_(&Trans[0], &Trans[0], &M, &N, &K, &alpha, A.GetMat(), &LDA, rho.GetMat(), &LDB, &beta, ap.GetMat(), &LDC);
	zher2k_(&UpLo[1], &Trans[2], &N, &K, &alpha, A.GetMat(), &LDA, ap.GetMat(), &LDC, &gamma, c.GetMat(), &LDC);
	zher2k_(&UpLo[1], &Trans[0], &N, &K, &alpha, A.GetMat(), &LDA, ap.GetMat(), &LDC, &delta, c.GetMat(), &LDC);
	// cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, M,N, K, &alpha, A.GetMat(), LDA, rho.GetMat(), LDB, &beta, ap.GetMat(), LDC);
	// cblas_zher2k(CblasColMajor, CblasLower, CblasConjTrans, N, K, &alpha, A.GetMat(), LDA, ap.GetMat(), LDC, 0.0, c.GetMat(), LDC);
	// cblas_zher2k(CblasColMajor, CblasLower, CblasNoTrans, N, K, &alpha, A.GetMat(), LDA, ap.GetMat(), LDC, -1.0, c.GetMat(), LDC);
	dataout << '\t' << to_ms*(clock() - clo);
	// for (size_t i=0; i< K; i++)
	// 	for(size_t j=i+1; j<K; j++)
	// 		n(i,j)=conj(n(j,i));
	// cout << c << endl;
	
	// //Method 3
	clo = clock();
	c=QOs::SuperDamp(A,rho);
	dataout << '\t' << to_ms*(clock() - clo);
	// //cout << c << endl;
	
	//Method 4
	clo = clock();
	matrix<complex<double> > nn(dim,dim);
	complex<double> I(0.0,-1.0);
	nn=MOs::Dagger(A)*A;
	c=2.0*A*rho*MOs::Dagger(A)-nn*rho-rho*nn;
	dataout << '\t' << to_ms*(clock() - clo) << endl;
	
	
	return 0;
}