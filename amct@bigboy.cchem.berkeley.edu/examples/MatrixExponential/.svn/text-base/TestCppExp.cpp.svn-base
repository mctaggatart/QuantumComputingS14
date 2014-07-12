/*
Name: TestCppExp.cpp
Author: Jay Gambetta

Dependences: MatrixExponential.hpp
Brief Discription: This program benchmarks my c++ matrix expontential
Limitations: None

Version History
	v0: Feburary 4th, 2009

run in terminal
	i="0"; while [ $i -lt 12 ]; do i=$[$i+1]; echo $i; ./TestCppExp.e $i TestCppExp.dat; done

*/
#include <MatrixExponential.hpp>
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
	complex<double> I(0.0,1.0);
	
	/* initialize random seed: */
	srand ( time(NULL) );
	
	matrix<complex<double> > A(dim,dim);
	
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
	
	int clo = clock();
	c=ExpM::EigenMethod(A,-I);
	dataout << '\t' << to_ms*(clock() - clo) << endl;
	
	
	return 0;
}