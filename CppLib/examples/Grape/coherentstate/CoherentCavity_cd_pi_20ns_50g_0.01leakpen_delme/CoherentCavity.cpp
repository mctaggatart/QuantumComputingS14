/*
Name: GrapeCabity.cpp
Author: fmotzoi

Dependences: GrapeUnitary.hpp  
Brief Discription: This demonstrates grape for lab frame vs rotating frame with large pixels

Version History
	v0: Feb  11th, 2009.
	v1: Oct 8, 2010
*/
#include <SweepOptimization.hpp>

using namespace std;

int main (int argc, char const *argv[]){

	verbose=no;
	cout << "Running program " << argv[0] << endl;

	char sourcename[80], cmd[80];
	strcpy(sourcename, argv[0]);
	strcpy(cmd, "cp ");
	sourcename[strlen(argv[0])-1]='\0';
	strcat(cmd, sourcename);
	strcat(cmd, "cpp ");
	cout << cmd << endl;
	
	

	//Grape inputs
	double wa=0*2*M_PI;  
	double wr = -1*2*M_PI; 		
	double g = 0.05*2*M_PI;
	size_t max_iter=1000000;
	complex<double> ii(0.0,1.0);
	double tolerance=0.0000000000000000000001, fidelity, base_a=2.0, epsilon=1, sigma=1, tgate=20, dt;
	size_t rwasub=40, mult=1;
	size_t num_time=tgate*rwasub*mult, dimQ=2, dimCav=2, dim= dimQ*dimCav, num_controls =3;
	dt=tgate/double(num_time);
	cout << dim << endl;
	SweepOptimization sys(dim, num_controls, num_time, dt, "CoherentCavity_cd_pi_20ns_50g_0.01leakpen_delme");
	strcat(cmd, sys.filename_);
	system(cmd);
	sys.SetNumericalParameters(fidelity=0.000000000000000000000999999999, base_a, epsilon, tolerance, max_iter);
	sys.Phi = &OptimizeEvolution::Phi4Sub2;
	sys.gradPhi = &OptimizeEvolution::GradPhi4Sub2;
	sys.Penalty = &OptimizeEvolution::LeakagePenalty;
	sys.gradPenalty = &OptimizeEvolution::GradLeakagePenalty;
	
		
				
	matrix<complex<double> > HcontrolZ(dim,dim), HcontrolX(dim,dim), HcontrolY(dim,dim);
	matrix<complex<double> > U_desired(dim,dim),  U_rot(dim,dim);
	U_desired.SetOutputStyle(Matrix);
	MOs::Identity(U_desired);
	
	
	matrix<complex<double> > Hdrift(dim, dim), Hcouple(dim, dim), HdriftQ(dimQ, dimQ),  HdriftCav(dimCav, dimCav);  
	matrix<complex<double> > IdentQ(dimQ, dimQ), IdentCav(dimCav, dimCav);
	MOs::Identity(IdentQ); MOs::Identity(IdentCav);
	
	matrix<complex<double> > a(dimCav,dimCav), ad(dimCav,dimCav), n(dimCav,dimCav);
	matrix<complex<double> > o(dimQ,dimQ), od(dimQ,dimQ), nq(dimQ,dimQ);
	
	Hdrift.SetOutputStyle(Matrix); HdriftCav.SetOutputStyle(Matrix); HdriftQ.SetOutputStyle(Matrix);
	HcontrolX.SetOutputStyle(Matrix);HcontrolY.SetOutputStyle(Matrix);HcontrolZ.SetOutputStyle(Matrix);
	
	MOs::Destroy(a); ad = MOs::Dagger(a); n=ad*a;
	MOs::Destroy(o); od = MOs::Dagger(o); nq=od*o;
		
	HdriftCav = wr*n;
	double Delta = -0.330 * 2 * M_PI;
	HdriftQ = wa*nq + Delta/2*( nq*nq - 1.0*nq);
	Hcouple = ( MOs::TensorProduct(a, od) + MOs::TensorProduct(ad, o) );
	Hdrift =  MOs::TensorProduct(HdriftCav, IdentQ) + MOs::TensorProduct(IdentCav, HdriftQ) + 50.0*Hcouple;
	
	U_rot = ExpM::EigenMethod(Hdrift-Hcouple, -ii*tgate);
	U_rot.SetOutputStyle(Matrix);

//	dimQ=0;
//	U_desired(dimQ+1,0)=1;//std::complex<double>(0,1/sqrt(2));
//	U_desired(0,dimQ+1)=1;//std::complex<double>(0,1/sqrt(2));
//	U_desired(0,0)=U_desired(dimQ+1,dimQ+1)=0;//1/sqrt(2);
//	U_desired= ExpM::EigenMethod(Hdrift, -ii*tgate)*U_desired;
//	U_desired(dimQ+1,dimQ+1)=-1;
	U_desired =  U_rot*MOs::TensorProduct(IdentCav, o + od);
//	U_desired =  U_rot*MOs::TensorProduct(IdentCav, IdentQ);

	
	HcontrolX = MOs::TensorProduct(0.5*(a + ad), IdentQ);
	HcontrolY = MOs::TensorProduct(0.5*(-ii*a+ii*ad), IdentQ);
	HcontrolZ = MOs::TensorProduct(IdentCav, nq);	
		
	AnalyticControl u0(HcontrolX, dt, num_time, rwasub, &sys, "cavcontrol0"); 
	//u0.ShiftedGaussian(M_PI, 2.5, NULL);
	u0.setparams(&AnalyticControl::RandomControl, 0, 0.01,  0, NULL, NULL);
	//u0.RandomControl( 0, 1);
	u0.pInterpolate = &Control::GaussianFilter;//CubicInterpolate;//
	u0.pSubGradients = &Control::GaussianFilterGradient;//setLinearSplineGradients;//
	AnalyticControl	u1(HcontrolY, dt, num_time, rwasub, &sys, "cavcontrol1");
	u1.setparams(&AnalyticControl::RandomControl, -0.01, 0.01,  0, NULL, NULL);
	u1.pInterpolate = &Control::GaussianFilter;//CubicInterpolate;//
	u1.pSubGradients = &Control::GaussianFilterGradient;//setLinearSplineGradients;//
	AnalyticControl	u2(HcontrolZ, dt, num_time, num_time, &sys, "cavcontrol2");
	u2.pInterpolate = &Control::GaussianFilter;//CubicInterpolate;//
	u2.pSubGradients = &Control::GaussianFilterGradient;//setLinearSplineGradients;//
	
	sys.SetHdrift(Hdrift);
	sys.SetRhoDesired(U_desired);
	
	sys.sweepinitialconditions(); 		
	sys.SetRhoDesired(o + od);	

	sys.U_[num_time-1].SetOutputStyle(Matrix);
	cout << "cav traced" << endl;
	cout << sys.U_[num_time-1] << endl;

	matrix<complex<double> > UA = (1.0/3)*MOs::TraceOutA(sys.U_[num_time-1], 3);
	
	UA.SetOutputStyle(Matrix);
	cout << "cav traced" << endl;
	cout << UA << endl;

	
	//cout << endl << U_desired << endl;
	//cout << endl << o + od << endl;
	
	cout << "fidelity traced " << sys.Phi4TrCav() << endl;
	
	
	
	return 0;
	
}