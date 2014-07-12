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
#include <MatrixOperations.hpp>
#include <QuantumOperations.hpp>
#include <OrdinaryDiffEqs.hpp>
using namespace std;

class ODE_equations_child: public ODE_equations<matrix<complex<double> > >{
  public:
	double kappa;
	double epsilon;
	double Delta;
	complex<double> I;
	matrix<complex<double> > a, ad, H;
	ODE_equations_child(size_t dim);
	~ODE_equations_child(){ cout << "~ODE_equations_child\n";};
	void derivs(double t, const matrix<complex<double> >& rho, matrix<complex<double> >& drho);
	void SaveSubset(std::vector<double>& data, const matrix<complex<double> >& rho);
};
inline ODE_equations_child::ODE_equations_child(size_t dim): a(dim,dim), H(dim,dim) {
	I=complex<double>(0.0,1.0);
	MOs::Destroy(a);
	ad = MOs::Dagger(a);
	
}
inline void ODE_equations_child::derivs(double t, const matrix<complex<double> >& rho, matrix<complex<double> >& drho){
	//returns the lower triangle of the computation
	drho = QOs::SuperHami(H,rho)+0.5*kappa*QOs::SuperDamp(a,rho);
	//Fill in the upper triangle of rho
	for (size_t i=0; i< rho.GetRows(); i++)
		for(size_t j=i+1; j<rho.GetColumns(); j++)
			drho(i,j)=conj(drho(j,i));
}
inline void ODE_equations_child::SaveSubset(std::vector<double>& data, const matrix<complex<double> >& rho){
	//none defined
}

int main (int argc, char const *argv[])
{	
	// verbose=no;
	size_t dim = 10;
	matrix<complex<double> > rho(dim,dim), rho_temp(dim,dim);
	rho_temp.SetOutputStyle(Matrix);
	
	//Initial condition
	double tinitial=0.0, tend = 10.0;
	complex<double> a0=0.0;
	size_t m =0;
	QOs::FockState(rho,m);
	
	
	//Constructing the system
	ODE_equations_child sys(dim);
	ODE_equations<matrix<complex<double> > >* pointersys = &sys;
	sys.kappa = 1.0; 
	sys.epsilon = 1.0;
	sys.Delta = 1.0;
	sys.H=sys.epsilon*(sys.a+sys.ad)+ sys.Delta*sys.ad*sys.a;
	
	//Files 
	ofstream data_out;
	string filename = "TestMasterRKCK.dat";
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
	asteady=-I*sys.epsilon/(sys.kappa*0.5+I*sys.Delta);

	for (size_t l=0; l < 8; l++){
		//Numerical parameters    	
		eps_abs = pow(10,(-double(l)));
		eps_rel = eps_abs; 
		h = eps_abs;
			
		//Solver declaring and solving
		ODE_solver<matrix<complex<double> > > solver(tinitial, tend, rho, pointersys, savepoints, filename_temp);
		solver.SetNumericalParameters(Adaptive, h, eps_abs, eps_rel, scale_y, scale_dydt);   
		solver.SetDataSaveType(Testing);
		solver.RKCKAdapStep();
	
		// Reading in the solution
		datain.open(filename_temp.c_str(), ios::in);
		UFs::FileCheck(datain,filename_temp.c_str());
	
		//Testing agains known soultion
		for (size_t k=0; k<savepoints;k++){
			datain >> t;
			datain >> h;
			// for (size_t i=0; i< rho.GetRows(); i++)
			// 	for(size_t j=0; j<rho.GetColumns(); j++)
			// 		datain >> rho_temp(i,j);
			datain >> rho_temp;
			numerical_a = MOs::Trace(sys.a*rho_temp);
			//cout << t <<  endl;
			exact_a= asteady + exp(-sys.kappa*t*0.5 -I*sys.Delta*t)*(a0-asteady);
			error = std::abs(exact_a  - numerical_a);
			// cout << error << endl; 
			data_out << eps_abs << " " << h << " " <<  t << " " <<  real(numerical_a) << " " << imag(numerical_a)  << " " << real(exact_a) << " " << imag(exact_a) << " " << error << endl;
		}
		cout << "Rho End is hermitian " << MOs::IsHermitian(rho_temp) << endl;
		datain.close();
		data_out << endl;  
	}	
	data_out.close();
	return 0;
}