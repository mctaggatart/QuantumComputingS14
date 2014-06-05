/*
Name: TestRK4.cpp
Author: Jay Gambetta

Dependences: OrdinaryDiffEqs.hpp  
Brief Discription: Demonstrates the use of my ODE solver using the RK4 method with a constant stepper
	This program test the Euler method for the differential equation
		d_t x = gamma x + epsilon
	The full solution is 
		x = -epsilon/gamma + exp(gamma*t) (x_0 + epsilon/gamma)

Version History
	v0: January 4th, 2009.

*/
#include <OrdinaryDiffEqs.hpp>
using namespace std;

class ODE_equations_child: public ODE_equations<vector<double> > {
  public:
	double gamma;
	double epsilon;
	ODE_equations_child(const double, const double);
	~ODE_equations_child(){ cout << "~ODE_equations_child\n";};
	void derivs(double t, const vector<double>& y, vector<double>& dy);
	void SaveSubset(std::vector<double>& data, const vector<double>& y);
};
inline ODE_equations_child::ODE_equations_child(const double gammain, const double epsilonin): gamma(gammain), epsilon(epsilonin) {}
inline void ODE_equations_child::derivs(double t, const vector<double>& y, vector<double>& dy){
	for (size_t j=0; j <y.size(); j++){
		dy[j] = gamma*y[j] + epsilon;
	}
}
inline void ODE_equations_child::SaveSubset(std::vector<double>& data, const vector<double>& y){
	//none defined
}

int main (int argc, char const *argv[])
{	
	// verbose=no;
	size_t dim = 1;
	vector<double> y(dim);
	
	//Initial condition
	double y0=0.0, tinitial=0.0, tend = 10.0, gamma=-1.0, epsilon=1.0;
	y[0]=y0;
	
	//Constructing the system
	ODE_equations_child sys(gamma,epsilon);
	ODE_equations<vector<double> >* pointersys = &sys;
	//sys.gamma = 
	
	//Files 
	ofstream data_out;
	string filename = "TestRK4.dat";
	UFs::OpenFile(filename, data_out, 16); 	// opens a file for saving the data
	ifstream datain;
	string filename_temp = "Temp.dat";
	
	//Numerical parameters
	size_t savepoints = 100;
	double h;
	
	// derived units for error testing
	double epsilon_on_gamma = epsilon/gamma;
	double A = (y0 + epsilon_on_gamma);
	double t, ny, ey, error;

	for (size_t l=0; l < 11; l++){
		//Numerical parameters    	
		h = pow(10,(-double(l)));
		
		//Solver declaring and solving
		ODE_solver<vector<double> > solver(tinitial, tend, y, pointersys, savepoints, filename_temp);
		solver.SetNumericalParameters(Constant, h) ;   
		solver.SetDataSaveType(Complete);
		solver.RK4ConstStep();
	
		// Reading in the solution
		datain.open(filename_temp.c_str(), ios::in);
		UFs::FileCheck(datain,filename_temp.c_str());
	
		//Testing agains known soultion
		for (size_t i=0; i<savepoints;i++){
			datain >> t;
			datain >> ny;
			//cout << t <<  endl;
			ey= -epsilon_on_gamma + exp(gamma*t)*A;
			error = abs(ey  - ny);
			// cout << error << endl; 
			data_out << h << " " <<  t << " " <<  ny  << " " << ey << " " << error << endl;
		}
		datain.close();
		data_out << endl;  
	}	
	data_out.close();
	return 0;
}