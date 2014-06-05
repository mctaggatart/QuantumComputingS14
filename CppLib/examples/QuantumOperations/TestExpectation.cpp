/*
Name: TestExpectation.cpp
Author: Jay Gambetta

Dependences: MatrixOperations.hpp, QuantumOperations.hpp
Brief Discription: This program benchmarks different wast of getting the expectation values.
Limitations: None

Version History
	v0: January 4th, 2009

run in terminal
	./TestExpectation.e 10 TestExpectation.dat
	i="0"; while [ $i -lt 12 ]; do i=$[$i+1]; echo $i; ./TestExpectation.e $i TestExpectation.dat; done

*/
#include <MatrixOperations.hpp>
#include <QuantumOperations.hpp>
using namespace std;

inline complex<double> Expectation(const matrix<complex<double> >& A, const matrix<complex<double> >& rho ){
   	size_t dim = rho.GetRows();
	if (dim != A.GetRows() || dim != A.GetColumns())
			UFs::MyError("matrix dimensions and rho dont agree");
 	complex<double> temp=0;
	for(size_t i=0; i< dim; i++){
		for(size_t p=0; p<dim; p++){
			temp += A(i,p)*rho(p,i);
		}
	}
	return temp;
}

int main (int argc, char const *argv[]){
	
	cout << "Running program " << argv[0] << endl;
	size_t dim=int(pow(2.0,atof(argv[1])));
	double to_ms =1000.0/CLOCKS_PER_SEC;
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
	
	complex<double> c;
	dataout << argv[1];
	
	//Method 1
	int clo = clock();
	c=MOs::Trace(A*rho);
	dataout << '\t' << to_ms*(clock() - clo);
	cout << c << endl;

	//Method 2
	clo = clock();
	c=Expectation(A,rho);
	dataout << '\t' << to_ms*(clock() - clo);
	cout << c << endl;
	
	//Method 3
	clo = clock();
	c=QOs::Expectation(A,rho);
	dataout << '\t' << to_ms*(clock() - clo) << endl;;
	cout << c << endl;
	
	return 0;
}