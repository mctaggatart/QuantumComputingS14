/*
Name: TestMasterCavity.cpp
Author: Jay Gambetta

Dependences: OrdinaryDiffEqs.hpp  
Brief Discription: Demonstrates the use of my ODE solver for simulating a master eqaution of a cavity

Limitations: None

Version History
	v0: January 4th, 2009.

run in terminal
	i="0"; while [ $i -lt 9 ]; do i=$[$i+1]; echo $i; ./TestMasterCavity.e $i TestMasterCavitySize.dat; done



*/
#include <MatrixOperations.hpp>
#include <QuantumOperations.hpp>
#include <OrdinaryDiffEqs.hpp>
using namespace std;

class ODE_equations_child: public ODE_equations<matrix<complex<double> > >{
  // defines the class for the equations
  public:
	double kappa, epsilon, Delta;
	complex<double> I;
	matrix<complex<double> > a, ad, H;
	ODE_equations_child(size_t dim);
	~ODE_equations_child(){ cout << "~ODE_equations_child\n";};
	void derivs(double t, const matrix<complex<double> >& rho, matrix<complex<double> >& drho);
	void SaveSubset(std::vector<double>& data, const matrix<complex<double> >& rho);
};
inline ODE_equations_child::ODE_equations_child(size_t dim): a(dim,dim), H(dim,dim) {
	//constructor for the class
	I=complex<double>(0.0,1.0);
	MOs::Destroy(a);
	ad = MOs::Dagger(a);
	
}
inline void ODE_equations_child::derivs(double t, const matrix<complex<double> >& rho, matrix<complex<double> >& drho){
	//the master equation
	//returns the lower triangle of the computation
	drho = QOs::SuperHami(H,rho)+0.5*kappa*QOs::SuperDamp(a,rho);
	
	//Fill in the upper triangle of rho
	for (size_t i=0; i< rho.GetRows(); i++)
		for(size_t j=i+1; j<rho.GetColumns(); j++)
			drho(i,j)=conj(drho(j,i));
}
inline void ODE_equations_child::SaveSubset(std::vector<double>& data, const matrix<complex<double> >& rho){
	//none defined
	complex<double> temp = QOs::Expectation(a,rho);
	data.push_back(real(temp));
	data.push_back(imag(temp));
}

int main (int argc, char const *argv[])
{	
	verbose=no;
	cout << "Running program " << argv[0] << endl;
	size_t dim=size_t(pow(2.0,atof(argv[1])));
	string const outfile = argv[2];
	ofstream dataout;
	UFs::AppendFile(outfile,dataout, 16);
	dataout << argv[1];
	
	//Initial condition
	matrix<complex<double> > rho(dim,dim);
	double tinitial=0.0, tend = 10.0;
	size_t m =0;
	QOs::FockState(rho,m);
	
	int clo = clock();
	//Constructing the system
	ODE_equations_child sys(dim);
	ODE_equations<matrix<complex<double> > >* pointersys = &sys;
	sys.kappa = 1.0; 
	sys.epsilon = 1.0;
	sys.Delta = 1.0;
	sys.H=sys.epsilon*(sys.a+sys.ad)+ sys.Delta*sys.ad*sys.a;
	
	//Files 
	string filename = "TestMasterCavity.dat";
	
	//Numerical parameters  
	size_t savepoints = 100;
	double h;
	double eps_rel, eps_abs;
	double scale_y = 1.0, scale_dydt = 1.0;
	eps_abs = 1e-4;;
	eps_rel = eps_abs; 
	h = eps_abs;
	
	//Solver declaring and solving
	ODE_solver<matrix<complex<double> > > solver(tinitial, tend, rho, pointersys, savepoints, filename);
	solver.SetNumericalParameters(Adaptive, h, eps_abs, eps_rel, scale_y, scale_dydt);   
	solver.SetDataSaveType(Subset);
	solver.RKCKAdapStep();
	dataout << '\t' << (clock() - clo) << endl;
	
	return 0;
}