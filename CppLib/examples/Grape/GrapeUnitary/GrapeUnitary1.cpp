/*
Name: Grape.cpp
Author: Jay Gambetta

Dependences: GrapeUnitary.hpp  
Brief Discription: This demonstrates grape for unitary using Phi3 and doing a pi pulse on a qubit
Limitations: None

Version History
	v0: Feb  11th, 2009.

time ./GrapeUnitary1.e GrapeUnitary1.dat

*/
#include <OptimizeEvolution.hpp>
using namespace std;

int main (int argc, char const *argv[]){
	
	verbose=no;
	cout << "Running program " << argv[0] << endl;
	
	//Grape inputs
	size_t num_time=10000, dim = 2, num_controls =1, max_iter=1000;//10000
	double tolerance=std::numeric_limits<double>::min(), fidelity, base_a=1.5, epsilon=5000, tgate=4, dt, alpha;
	dt=tgate/double(num_time);
	
	
	OptimizeEvolution sys(dim, num_controls, num_time, dt, "Unitary1");
	sys.SetNumericalParameters(fidelity=0.9999900, base_a, epsilon, tolerance, max_iter);
	
	matrix<complex<double> > a(dim,dim), ad(dim,dim), n(dim,dim), Hdrift(dim,dim), Hcontrol(dim,dim);
	double delta=0;
	MOs::Destroy(a); ad = MOs::Dagger(a);
	n=ad*a;
	Hdrift = delta*n;
	Hcontrol=0.5*(a+ad);
	sys.SetHdrift(Hdrift);
 	
	Control u0(Hcontrol, dt, num_time, 1, &sys, "u1_control");
	u0.ShiftedGaussian(M_PI, alpha=2.5, NULL);
	u0.Normalize(M_PI/10);
	//Initial condition
	matrix<complex<double> > U_desired(dim,dim);
	U_desired.SetOutputStyle(Matrix);
	MOs::Null(U_desired);
		U_desired(1,0)=std::complex<double>(0,-1);
		U_desired(0,1)=std::complex<double>(0,-1);
	//	U_desired(0,1)=std::complex<double>(1,0);
	//U_desired(1,0)=std::complex<double>(1,0);
		//		sys.SetOmega(0.09);
		//		sys.SetOmega(.01); //good one
		sys.SetOmega(50.0);
		sys.SetUDesired(U_desired);
	sys.SetTrueRhoDesired(U_desired);
	//run grape	
	sys.UnitaryTransfer();
	
	
		
	return 0;
}
