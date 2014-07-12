/*
Name: GrapeCabity.cpp
Author: fmotzoi

Dependences: GrapeUnitary.hpp  
Brief Discription: This demonstrates grape for lab frame vs rotating frame with large pixels

Version History
	v0: Feb  11th, 2009.
	v1: Oct 8, 2010
*/
#include <OptimizeEvolution.hpp>
using namespace std;

int main (int argc, char const *argv[]){

	verbose=no;
	cout << "Running program " << argv[0] << endl;

	//Grape inputs
	double wa=0*2*M_PI;  
	double wr = -1*2*M_PI; 		
	double wr2 = 1*2*M_PI; 		
	double g = 0.05*2*M_PI;
	size_t max_iter=100000;
	complex<double> ii(0.0,1.0);
	double tolerance=0.0000000000000000000001, fidelity, base_a=2.0, epsilon=1, sigma=1, tgate=30, dt;
	size_t rwasub=10, mult=1;
	size_t num_time=tgate*rwasub*mult, dimQ=2, dimCav=3, dimCav2=2, dim= dimQ*dimCav*dimCav2, num_controls =4;
	dt=tgate/double(num_time);
	cout << dim << endl;
	OptimizeEvolution sys(dim, num_controls, num_time, dt, "TwoCavity");
	sys.SetNumericalParameters(fidelity=0.9999900, base_a, epsilon, tolerance, max_iter);
	
	matrix<complex<double> > HcontrolZ(dim,dim), HcontrolZ2(dim,dim), HcontrolX(dim,dim), HcontrolY(dim,dim);
	matrix<complex<double> > U_desired(dim,dim),  U_rot(dim,dim);
	U_desired.SetOutputStyle(Matrix);
	MOs::Identity(U_desired);
	
	
	matrix<complex<double> > Hdrift(dim, dim), Hcouple(dim, dim), HdriftQ(dimQ, dimQ),  HdriftCav(dimCav, dimCav);  
	matrix<complex<double> > IdentQ(dimQ, dimQ), IdentCav(dimCav, dimCav), IdentCav2(dimCav2, dimCav2);
	MOs::Identity(IdentQ); MOs::Identity(IdentCav);MOs::Identity(IdentCav2);
	
	matrix<complex<double> > a(dimCav,dimCav), ad(dimCav,dimCav), n(dimCav,dimCav);
	matrix<complex<double> > a2(dimCav2,dimCav2), ad2(dimCav2,dimCav2), n2(dimCav2,dimCav2);
	matrix<complex<double> > o(dimQ,dimQ), od(dimQ,dimQ), nq(dimQ,dimQ);
	
	Hdrift.SetOutputStyle(Matrix); HdriftCav.SetOutputStyle(Matrix); HdriftQ.SetOutputStyle(Matrix);
	HcontrolX.SetOutputStyle(Matrix);HcontrolY.SetOutputStyle(Matrix);HcontrolZ.SetOutputStyle(Matrix);HcontrolZ2.SetOutputStyle(Matrix);
	
	MOs::Destroy(a); ad = MOs::Dagger(a); n=ad*a;
	MOs::Destroy(a2); ad2 = MOs::Dagger(a2); n2=ad2*a2;
	MOs::Destroy(o); od = MOs::Dagger(o); nq=od*o;
	U_desired= MOs::TensorProduct(IdentCav2, MOs::TensorProduct(IdentCav, 0.5*(o + od)));
		
	HdriftCav = wr*n;
	double Delta = -0.330 * 2 * M_PI;
	HdriftQ = wa*nq + Delta/2*( nq*nq - 1.0*nq);
	Hcouple = MOs::TensorProduct(IdentCav2,g * ( MOs::TensorProduct(a, od) + MOs::TensorProduct(ad, o) ));
	Hdrift =  MOs::TensorProduct(IdentCav2, MOs::TensorProduct(HdriftCav, IdentQ) + MOs::TensorProduct(IdentCav, HdriftQ)) + Hcouple;
	cout << "Hdrift\n" << Hdrift << endl;
	
	U_rot = ExpM::EigenMethod(Hdrift-Hcouple, -ii*tgate);
	U_rot.SetOutputStyle(Matrix);

//	dimQ=0;
//	U_desired(dimQ+1,0)=1;//std::complex<double>(0,1/sqrt(2));
//	U_desired(0,dimQ+1)=1;//std::complex<double>(0,1/sqrt(2));
//	U_desired(0,0)=U_desired(dimQ+1,dimQ+1)=0;//1/sqrt(2);
//	U_desired= ExpM::EigenMethod(Hdrift, -ii*tgate)*U_desired;
//	U_desired(dimQ+1,dimQ+1)=-1;
	U_desired =  U_rot*MOs::TensorProduct(IdentCav2, MOs::TensorProduct(IdentCav, IdentQ));


	
	HcontrolX = MOs::TensorProduct(IdentCav2, MOs::TensorProduct(IdentCav, 0.5*(o + od)));
	HcontrolY = MOs::TensorProduct(IdentCav2, MOs::TensorProduct(IdentCav, nq));//+MOs::TensorProduct(0.5*(-ii*a+ii*ad), IdentQ);
	HcontrolZ = MOs::TensorProduct(IdentCav2, (MOs::TensorProduct(n, IdentQ)));
	//HcontrolZ = MOs::TensorProduct(IdentCav2, (MOs::TensorProduct(n-0.5*IdentCav, IdentQ) + MOs::TensorProduct(IdentCav,nq- 0.5*IdentQ)));
	HcontrolZ2 = MOs::TensorProduct(0.5*n2, MOs::TensorProduct(IdentCav, IdentQ));
	
		
	AnalyticControl u0(HcontrolX, dt,num_time,rwasub,&sys, "cavcontrol0"); 
	AnalyticControl	u1(HcontrolY, dt, num_time, rwasub,&sys, "cavcontrol1");
	AnalyticControl	u2(HcontrolZ, dt, num_time, rwasub, &sys, "cavcontrol2");
	AnalyticControl	u3(HcontrolZ2, dt, num_time, rwasub,&sys, "cavcontrol3");

	sys.SetHdrift(Hdrift);
	sys.SetRhoDesired(U_desired);


	sys.UnitaryTransfer(); 
	
	return 0;
	
}