/*
Name: TestComplexRK45.cpp
Author: Jay Gambetta

Dependences: OrdinaryDiffEqs.hpp  
Brief Discription: Demonstrates the use of my ODE solver using the RK45 method with an adative stepper to solve a complex vector
	This program test the Euler method for the differential equation
		d_t a = - kappa a/2 -i Delta a -i epsilon
	The full solution is 
		a = a_s + exp(-kappa t/2 -i Delta t)(a_0-a_s)
	wher
		a_s = -i epsilon/(kappa/2+i Delta )

Limitations: None

Version History
	v0: January 4th, 2009.

*/
#include <OrdinaryDiffEqs.hpp>
using namespace std;

class ODE_equations_child: public ODE_equations<vector<complex<double> > >{
  public:
	double kappa;
	double epsilon;
	double Delta;
	complex<double> I;
	ODE_equations_child(const double kappain, const double epsilonin, const double Deltain);
	~ODE_equations_child(){ cout << "~ODE_equations_child\n";};
	void derivs(double t, const vector<complex<double> >& a, vector<complex<double> >& da);
	void SaveSubset(std::vector<double>& data, const vector<complex<double> >& a);
};
inline ODE_equations_child::ODE_equations_child(const double kappain, const double epsilonin, const double Deltain): kappa(kappain), epsilon(epsilonin), Delta(Deltain) {
	 I=complex<double>(0.0,1.0);
}
inline void ODE_equations_child::derivs(double t, const vector<complex<double> >& a, vector<complex<double> >& da){
	for (size_t j=0; j <a.size(); j++){
		da[j] = -0.5*kappa*a[j] - I*epsilon -I*Delta*a[j];
	}
}
inline void ODE_equations_child::SaveSubset(std::vector<double>& data, const vector<complex<double> >& a){
	//none defined
}


int main (int argc, char const *argv[])
{	
	// verbose=no;
	size_t dim = 1;
	vector<complex<double> > a(dim);
	
	//Initial condition
	double tinitial=0.0, tend = 10.0, kappa=1.0, epsilon=1.0, Delta=1.0;
	complex<double> a0=0.0;
	a[0]=a0;
	
	//Constructing the system
	ODE_equations_child sys(kappa,epsilon,Delta);
	ODE_equations<vector<complex<double> > >* pointersys = &sys;
	//sys.gamma = 
	
	//Files 
	ofstream data_out;
	string filename = "TestComplexRK45.dat";
	UFs::OpenFile(filename, data_out, 16); 	// opens a file for saving the data
	ifstream datain;
	string filename_temp = "Temp.dat";
	
	//Numerical parameters  
	size_t savepoints = 100;
	double h;
	double eps_rel, eps_abs;
	double scale_y = 1.0, scale_dydt = 1.0;
	
	// derived units for error testing
	double t, error;
	complex<double> numerical_a, exact_a, asteady, I(0.0, 1.0);
	asteady=-I*epsilon/(kappa*0.5+I*Delta);

	for (size_t l=0; l < 8; l++){
		//Numerical parameters    	
		eps_abs = pow(10,(-double(l)));
		eps_rel = eps_abs; 
		h = eps_abs;
			
		//Solver declaring and solving
		ODE_solver<vector<complex<double> > > solver(tinitial, tend, a, pointersys, savepoints, filename_temp);
		solver.SetNumericalParameters(Adaptive, h, eps_abs, eps_rel, scale_y, scale_dydt);   
		solver.SetDataSaveType(Testing);
		solver.RK45AdapStep();
	
		// Reading in the solution
		datain.open(filename_temp.c_str(), ios::in);
		UFs::FileCheck(datain,filename_temp.c_str());
	
		//Testing agains known soultion
		for (size_t i=0; i<savepoints;i++){
			datain >> t;
			datain >> h;
			datain >> numerical_a;
			//cout << t <<  endl;
			exact_a= asteady + exp(-kappa*t*0.5 -I*Delta*t)*(a0-asteady);
			error = std::abs(exact_a  - numerical_a);
			// cout << error << endl; 
			data_out << eps_abs << " " << h << " " <<  t << " " <<  real(numerical_a) << " " << imag(numerical_a)  << " " << real(exact_a) << " " << imag(exact_a) << " " << error << endl;
		}
		datain.close();
		data_out << endl;  
	}	
	data_out.close();
	return 0;
}