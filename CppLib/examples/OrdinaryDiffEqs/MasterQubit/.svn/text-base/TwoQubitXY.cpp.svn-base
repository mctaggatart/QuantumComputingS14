/*
Name: TwoQubitXY.cpp
Author: Jay Gambetta

Dependences: None  
Brief Discription: This numerically simulates the master equation for XY coupling of two qubits, it saves in the pauli basis as well as the purity, concurance and norm
Limitations: None

Version History
	v0: January  8, 2009.

*/
#include <QuantumOperations.hpp>
#include <OrdinaryDiffEqs.hpp>
using namespace std;

double	tanhpulse(double t, double sigma, double td, double ton, double amp){
	double one_on_sigma=1.0/sigma;
	return 0.5*amp*(tanh(one_on_sigma*(t-ton)) + tanh(-one_on_sigma*(t-ton-td)) );
}


class ODE_equations_child: public ODE_equations<matrix<complex<double> > >{
  public:
	double J, gamma1, gamma2, gamma_phi1, gamma_phi2, omega1, omega2, amp,tg;
	matrix<complex<double> > sii, sxi, syi, szi, six, sxx, syx, szx, siy, sxy, syy, szy, siz, sxz, syz, szz, spi, smi, sip, sim, H0, H;
	ODE_equations_child();
	~ODE_equations_child(){ cout << "~ODE_equations_child\n";};
	void derivs(double t, const matrix<complex<double> >& rho, matrix<complex<double> >& drho);
	void SaveSubset(std::vector<double>& data, const matrix<complex<double> >& rho);
};
inline ODE_equations_child::ODE_equations_child(){
	
}
inline void ODE_equations_child::derivs(double t, const matrix<complex<double> >& rho, matrix<complex<double> >& drho){
	//returns the lower triangle of the computation
	H = H0 + (omega1 -tanhpulse(t,0.001,tg,0.01,amp))*szi*0.5+ omega2*siz*0.5;
	drho = QOs::SuperHami(H,rho)+0.5*gamma1*QOs::SuperDamp(smi,rho)+0.5*gamma2*QOs::SuperDamp(sim,rho)+0.25*gamma_phi2*QOs::SuperDamp(szi,rho)+0.25*gamma_phi1*QOs::SuperDamp(siz,rho);
	//Fill in the upper triangle of rho
	for (size_t i=0; i< rho.GetRows(); i++)
		for(size_t j=i+1; j<rho.GetColumns(); j++)
			drho(i,j)=conj(drho(j,i));
}
inline void ODE_equations_child::SaveSubset(std::vector<double>& data, const matrix<complex<double> >& rho){
	//none defined
	// complex<double> temp;
	double temp2;
	    
	// temp = QOs::Expectation(sii,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(sxi,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(syi,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(szi,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(six,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(sxx,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(syx,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(szx,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(siy,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(sxy,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(syy,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(szy,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(siz,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(sxz,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(syz,rho);
	// data.push_back(real(temp));
	// temp = QOs::Expectation(szz,rho);
	// data.push_back(real(temp));
	
	temp2 = real(MOs::Trace(rho));
	data.push_back(temp2);
	
	temp2 = QOs::Purity(rho);
	data.push_back(temp2);
	
	temp2 = QOs::Concurrence(rho);
	data.push_back(temp2);
		
}

int main (int argc, char const *argv[])
{
	matrix<complex<double> > sx(2,2), sy(2,2), sz(2,2), si(2,2), sp(2,2), sm(2,2);
	MOs::Pauli(sx, sy, sz);
	MOs::Identity(si);
	MOs::Destroy(sp);
	sm =MOs::Dagger(sp);
	

	
	//Constructing the system
	ODE_equations_child sys;
	ODE_equations<matrix<complex<double> > >* pointersys = &sys;
	sys.J = 0.06*2.0*M_PI;
	sys.gamma1 =1.0/600.0;
	sys.gamma2 = 1.0/600.0;
	sys.gamma_phi1 = 0.0;
	sys.gamma_phi2 = 0.0;
	sys.omega1 = 9*2.0*M_PI-8*2.0*M_PI;
	sys.omega2 = 8*2.0*M_PI-8*2.0*M_PI;
	sys.amp=2.0*M_PI;
	sys.tg = 500.0;
	
	//All the qubit-cavity operators
	sys.sii=MOs::TensorProduct(si,si);
	sys.sxi=MOs::TensorProduct(sx,si);
	sys.syi=MOs::TensorProduct(sy,si);
	sys.szi=MOs::TensorProduct(sz,si);
	sys.six=MOs::TensorProduct(si,sx);
	sys.sxx=MOs::TensorProduct(sx,sx);
	sys.syx=MOs::TensorProduct(sy,sx);
	sys.szx=MOs::TensorProduct(sz,sx);
	sys.siy=MOs::TensorProduct(si,sy);
	sys.sxy=MOs::TensorProduct(sx,sy);
	sys.syy=MOs::TensorProduct(sy,sy);	
	sys.szy=MOs::TensorProduct(sz,sy);
	sys.siz=MOs::TensorProduct(si,sz);
	sys.sxz=MOs::TensorProduct(sx,sz);
	sys.syz=MOs::TensorProduct(sy,sz);
	sys.szz=MOs::TensorProduct(sz,sz);
	
	sys.spi=MOs::TensorProduct(sp,si);
	sys.smi=MOs::TensorProduct(sm,si);
	sys.sip=MOs::TensorProduct(si,sp);
	sys.sim=MOs::TensorProduct(si,sm);
	
	sys.H0 = sys.J*(sys.sim*sys.spi + sys.sip*sys.smi);
	
	matrix<complex<double> > rho1(2,2), rho2(2,2), rho(4,4);
	//Initial condition
	double tinitial=0.0, tend = 50.0;
	QOs::QubitState(rho1,0,0,1);
	QOs::QubitState(rho2,0,0,-1);
	rho=MOs::TensorProduct(rho1,rho2);
	
	
	//Files 
	string filename = "TestQubitXY.dat";
	
	//Numerical parameters  
	size_t savepoints = 1000;
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
	
	
		
	
	return 0;
}
