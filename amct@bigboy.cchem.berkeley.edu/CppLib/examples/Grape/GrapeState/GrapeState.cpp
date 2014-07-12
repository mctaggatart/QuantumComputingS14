/*
Name: Grape.cpp
Author: Jay Gambetta

Dependences: GrapeState.hpp  
Brief Discription: This demonstrates grape for state transformation for a nonlinear oscilator 
Limitations: None

Version History
	v0: Feb  11th, 2009.

time ./GrapeState.e GrapeState.dat

*/
#include <Grape.hpp>
using namespace std;

int main (int argc, char const *argv[]){
	
	verbose=no;
	cout << "Running program " << argv[0] << endl;
	
	//Grape inputs
	size_t num_time=1000, dim = 3, num_controls =1;
	Grape sys(dim, num_controls, num_time);
	size_t max_iter=100000;
	double tolerance=std::numeric_limits<double>::min(), fidelity=0.999, base_a=1.3, epsilon=5000, tgate=10*2*M_PI, dt;
	dt=tgate/double(num_time);
	size_t numsubpix[2];
	numsubpix[0] =1;
	numsubpix[1] =1;

	// cout << dt << endl;
	sys.SetNumericalParameters(fidelity, dt, base_a, epsilon, tolerance, max_iter);
	
	matrix<complex<double> > a(dim,dim), ad(dim,dim), n(dim,dim), Hdrift(dim,dim), Hcontrol(dim,dim);
	double omega=1, delta=10.0;
	MOs::Destroy(a);
	ad = MOs::Dagger(a);
	n=ad*a;
	Hdrift = (omega-0.5*delta)*n+ 0.5*delta*n*n;
	Hcontrol=a+ad;
	sys.SetHdrift(Hdrift);
 	sys.SetHcontrol(Hcontrol,0);
	vector<double> ucontrol(num_time), ucontrolfinal;
	for(size_t j = 0; j < num_time; ++j)
	{
		ucontrol[j]=M_PI*cos(omega*j*dt)/tgate*0.99;
	}
	sys.Setucontrol(ucontrol,0);
	
	//Initial condition
	matrix<complex<double> > rho_initial(dim,dim), rho_desired(dim,dim);
	QOs::FockState(rho_initial,0);
	QOs::FockState(rho_desired,1);
	sys.SetRhoDesired(rho_desired);
	sys.SetRhoInitial(rho_initial);
	
	//run grape
	sys.StateTransfer(&Grape::Phi0, &Grape::GradPhi0, NULL, NULL);
	
	
	// The finial control
	ucontrolfinal.assign(sys.Getucontrol(0),sys.Getucontrol(0) + num_time );
	string const outfile = argv[1];
	ofstream dataout;
	UFs::OpenFile(outfile,dataout, 16);
	for(size_t j =0; j < num_time; j++){
			dataout << j*dt << '\t' << ucontrolfinal[j] << endl;
	}
	dataout.close();
	return 0;
}