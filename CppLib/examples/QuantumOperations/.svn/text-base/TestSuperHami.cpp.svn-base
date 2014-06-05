/*
Name: TestSuperHami.cpp
Author: Jay Gambetta

Dependences: MatrixOperations.hpp, QuantumOperations.hpp
Brief Discription: This program benchmarks different wast of getting the SuperOperator for -i[H,rho].
Limitations: None

Version History
	v0: January 4th, 2009

run in terminal
	i="0"; while [ $i -lt 10 ]; do i=$[$i+1]; echo $i; ./TestSuperHami.e $i TestSuperHami.dat; done

*/
#include <MatrixOperations.hpp>
#include <QuantumOperations.hpp>
using namespace std;

inline matrix<complex<double> > SuperHami(const matrix<complex<double> >& H, const matrix<complex<double> >& rho){
	size_t dim= rho.GetRows(), rows1 = H.GetRows(), cols1 = H.GetColumns();
	if (dim != rows1 || dim != cols1)
		UFs::MyError("dimensions of operators and rho dont agree");
		
	complex<double> I(0.0,1.0);
	matrix<complex<double> > temprho(dim,dim);
	// MOs::Null(temprho);
	for (size_t i=0; i < dim; i++){
		for(size_t p=0; p < dim; p++){
			//diagonal part
			temprho(i,i) += -I*H(i,p)*rho(p,i) + I*rho(i,p)*H(p,i) ;
			//cout << imag(tempdia) << endl;
			//lowerTriangle
			for(size_t j=0; j < i; j++){
				temprho(i,j)+=-I*H(i,p)*rho(p,j) + I*rho(i,p)*H(p,j);;
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
	//c.SetOutputStyle(Matrix);
	//A.SetOutputStyle(Matrix);
	//rho.SetOutputStyle(Matrix);
	dataout << argv[1];
	
	cout << "here" << endl;
	//Method 1
	int clo = clock();
	c=SuperHami(A,rho);
	dataout << '\t' << to_ms*(clock() - clo);
	//cout << A << endl;
	//cout << rho << endl;
	//cout << c << endl;

	//Method 2
	clo = clock();
	complex<double> alpha(0.0,-1.0); 
	double beta =0.0;
	size_t N=A.GetRows();
	size_t K=A.GetColumns();
	
	size_t LDA = N;
	size_t LDB = N;
	size_t LDC = N;
	zher2k_(&UpLo[1], &Trans[0], &N, &K, &alpha, A.GetMat(), &LDA, rho.GetMat(), &LDB, &beta, c.GetMat(), &LDC);
	//cblas_zher2k(CblasColMajor, CblasLower, CblasNoTrans, N, K, &alpha, A.GetMat(), LDA, rho.GetMat(), LDB, beta, c.GetMat(), LDC);
	dataout << '\t' << to_ms*(clock() - clo);
	//cout << c << endl;
	
	//Method 3
	clo = clock();
	c=QOs::SuperHami(A,rho);
	dataout << '\t' << to_ms*(clock() - clo);
	//cout << c << endl;
	
	//Method 4
	clo = clock();
	complex<double> I(0.0,-1.0);
	c=-I*A*rho+I*rho*A;
	dataout << '\t' << to_ms*(clock() - clo) << endl;
	
	
	return 0;
}