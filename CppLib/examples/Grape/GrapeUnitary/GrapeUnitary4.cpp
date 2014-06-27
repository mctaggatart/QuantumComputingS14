/*
Name: Grape.cpp
Author: Jay Gambetta

Dependences: GrapeUnitary.hpp  
Brief Discription: This demonstrates grape for unitary for PHI_4 on a three level system like the transmon

Version History
	v0: Feb  11th, 2009.

time ./GrapeUnitary4.e GrapeUnitary4.dat


n = 10000
real	0m8.287s
user	0m7.917s
sys	0m0.134s

real	0m0.027s
user	0m0.023s
sys	0m0.003s


*/
#include <OptimizeEvolution.hpp>
using namespace std;

int main (int argc, char const *argv[]){

	verbose=no;
	cout << "Running program " << argv[0] << endl;

	//Grape inputs
	size_t num_time=10, dim = 3, num_controls =2;
	size_t max_iter=10000;
	double tolerance=std::numeric_limits<double>::min(), fidelity, base_a=2.0, epsilon=1, sigma=1, tgate=4, dt;
	dt=tgate/double(num_time);
	
	
	OptimizeEvolution sys(dim, num_controls, num_time, dt, "Unitary4");
	sys.SetNumericalParameters(fidelity=0.9999900, base_a, epsilon, tolerance, max_iter);
	
	matrix<complex<double> > a(dim,dim), ad(dim,dim), n(dim,dim), Hdrift(dim,dim), Hcontrol(dim,dim);
	double delta=0.0, Delta=-2.02319;
	MOs::Destroy(a);
	ad = MOs::Dagger(a);
	n=ad*a;
	Hdrift = (delta-0.5*Delta)*n+ 0.5*Delta*n*n;
	sys.SetHdrift(Hdrift);
	
	AnalyticControl u0((a+ad)*0.5, dt, num_time, 1, &sys, "u4_control0");
	u0.ShiftedGaussian(M_PI, tgate/sigma, NULL);
	
	Control u1((-complex<double>(0.0,1.0)*a+complex<double>(0.0,1.0)*ad)*0.5, dt, num_time, 1, &sys, "u4_control1");
	
	matrix<complex<double> > U_desired(dim,dim);
	U_desired.SetOutputStyle(Matrix);
	MOs::Null(U_desired);
	U_desired(1,0)=std::complex<double>(1,0);
	U_desired(0,1)=std::complex<double>(1,0);
	// U_desired(2,2)=std::complex<double>(1,0);
	sys.SetUDesired(U_desired);	
	sys.SetTrueRhoDesired(U_desired);
	sys.Phi = &OptimizeEvolution::Phi4Sub2;
	sys.gradPhi = &OptimizeEvolution::GradPhi4Sub2;	
	
	//run grape
	double to_ms =1000.0/CLOCKS_PER_SEC;
	int clo ;
	clo = clock();
	sys.UnitaryTransfer();
	clo = to_ms*(clock() - clo);
	cout << "time " << clo << endl;

	return 0;
}
