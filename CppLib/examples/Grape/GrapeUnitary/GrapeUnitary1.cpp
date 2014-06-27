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
	size_t num_time=1000, dim = 2, num_controls =1, max_iter=100000;
		double tolerance= 10*numeric_limits<double>::min();
	  double fidelity, base_a=1.5, epsilon=5000, tgate=4, dt, alpha;
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
 	//changed to reg control from analytic
	Control u0(Hcontrol, dt, num_time, 1, &sys, "u1_control");
		u0.ShiftedGaussian(M_PI, alpha=2.5, NULL);
	//edited random control to not actually be random
	//	u0.RandomControl(0,1);
	//u0.Sine(0,0.1);
	//Initial condition
       	matrix<complex<double> > U_desired(dim,dim);
	matrix<complex<double> > Rho_des(dim,dim);
	Rho_des.SetOutputStyle(Matrix);
	U_desired.SetOutputStyle(Matrix);
	MOs::Null(U_desired);
	MOs::Null(Rho_des);
	Rho_des(1,1)=1;
	//	U_desired(1,1)=1;
	
	U_desired(1,0)=std::complex<double>(0,-1);
	U_desired(0,1)=std::complex<double>(0,-1);
	sys.SetUDesired(U_desired);
	sys.SetTrueRhoDesired(U_desired);
	//run grape	
	sys.UnitaryTransfer();
	
	
		
	return 0;
}
