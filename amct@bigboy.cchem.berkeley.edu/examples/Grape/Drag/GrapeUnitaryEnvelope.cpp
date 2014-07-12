/*
Name: GrapeUnitaryEnvelope.cpp
Author: felix motzoi

Dependences: GrapeUnitaryEnvelope.hpp  
Brief Discription: This demonstrates grape for lab frame vs rotating frame with large pixels

Version History
	v0: Using AnalyticControl class, Oct. 6 2009.

*/
#include "OptimizeEvolution.hpp"

using namespace std;

int main (int argc, char const *argv[]){

	verbose=no;
	cout << "Running program " << argv[0] << endl;
	
	//System/Optimization parameters
	double qbfreq=51.5, Delta=-2.5, lambda=sqrt(3/2);  //physical params
	double lamb2 =  1/sqrt(2);
	double rab = M_PI;
	double Delta2 = Delta;
	double tolerance=0.00000000000000001, fidelity=0.9, base_a=1.3, epsilon=0.000005, tgate=1, dt;  //search params
	size_t max_iter=500000;
	size_t dim=5;
	double to_ms =1000.0/CLOCKS_PER_SEC;
	int clo;
	int buffer = 1;

	
	//Typical operators	
	complex<double> ii(0.0,1.0);
	matrix<complex<double> > a(dim,dim), ad(dim,dim), n(dim,dim), UNOT(dim,dim), Ident(dim, dim), twophot(dim, dim);
	matrix<complex<double> > Urot(dim,dim);
	MOs::Destroy(a);
	ad = MOs::Dagger(a);
	n=ad*a;
	MOs::Identity(Ident);
	MOs::Null(twophot);
	twophot(0,2)=1;//-ii;
	twophot(2,0)=1;//ii;


//	Control ctrl(a, 1, 400, 1, NULL, "T0.dat");
//	ctrl.readfile("T0.txt");
//	ctrl.DFT();
//	ctrl.writefile();
	
/*	
	size_t subpix = 2000/5*1;
	Control ctrl(a, 1.0/subpix, subpix*(5+2*buffer), subpix, NULL, "T6_pix.dat");
	ctrl.buffer = buffer;
	ctrl.u_[3*subpix/2]=0.15;
	ctrl.u_[5*subpix/2]=0.10;
	ctrl.u_[7*subpix/2]=0.05;
	ctrl.u_[9*subpix/2]=0.176;
	ctrl.u_[11*subpix/2]=-0.05;
	
	ctrl.pInterpolate =  &AnalyticControl::Pixelate;//
	ctrl.Interpolate();	
	
	ctrl.writefile();
	
	ctrl.pInterpolate =  &AnalyticControl::GaussianFilter; //CubicInterpolate;//Pixelate;//

	ctrl.Interpolate();	
	strcpy(ctrl.filename, "T6_filter.dat");
	ctrl.writefile();

	delete [] ctrl.u_filt;
	ctrl.u_filt = ctrl.u_;

//	ctrl.u_[subpix/2]=0.12;
//	ctrl.u_[3*subpix/2]=0.08;
//	ctrl.u_[5*subpix/2]=0.04;
//	ctrl.u_[7*subpix/2]=0.14;
//	ctrl.u_[9*subpix/2]=-0.04;
	ctrl.u_[3*subpix/2]=0.15;
	ctrl.u_[5*subpix/2]=0.10;
	ctrl.u_[7*subpix/2]=0.05;
	ctrl.u_[9*subpix/2]=0.176;
	ctrl.u_[11*subpix/2]=-0.05;
//	ctrl.u_[subpix/2]=0.15*0.5+0.05;
//	ctrl.u_[3*subpix/2]=0.10*0.5+0.05;
//	ctrl.u_[5*subpix/2]=0.05*0.5+0.05;
//	ctrl.u_[7*subpix/2]=0.176*0.5+0.05;
//	ctrl.u_[9*subpix/2]=-0.05*0.5+0.05;

	
	strcpy(ctrl.filename, "T6_spline.dat");
	ctrl.pInterpolate =  &AnalyticControl::Pixelate;//CubicInterpolate;//Pixelate;//		
	ctrl.Interpolate();	
	ctrl.pInterpolate =  &AnalyticControl::CubicInterpolate;//CubicInterpolate;//Pixelate;//		
	ctrl.Interpolate();	
	
	ctrl.writefile();
	/*

	Control ctrl1(a, 1.0/subpix, subpix*(5+2*buffer), subpix, NULL, "T6.dat");
	ctrl1.buffer = buffer;
	ctrl1.readfile("T6.txt");
	ctrl1.writefile();


	return 0;
*/

	//RWA fixed parameters
	const size_t rwasub=20;
	size_t rwa_num_time=tgate*rwasub, num_controls =4;

	//RW fixed parameters
	//const int nsub=20;
	//size_t num_time=rwa_num_time*nsub/rwasub;
	//dt=tgate/double(num_time);
		system("mkdir hello");
	
		//RWA frame	
		double rwa_dt = tgate/double(rwa_num_time);
		OptimizeEvolution sys(dim, num_controls, rwa_num_time, rwa_dt, "anharmonic_ladder");
		
		sys.Phi = &OptimizeEvolution::Phi4Sub2;
		sys.gradPhi = &OptimizeEvolution::GradPhi4Sub2;
				
		//sys.SetHdrift(Delta/2*(n*n-n));			//first transition
		matrix<complex<double> > HDrift(dim, dim);
		//level 0
		HDrift = Delta/2*( n*n - n);
		MOs::GenPauliX(UNOT,0,1);

//change lambda here
//		a(1,2)=a(2,1)=0;//1.75;
//		ad = MOs::Dagger(a);

	ptrSetTransition pSetTransition = &AnalyticControl::setdriving01;
	//level1
	//	HDrift =  Delta/2*( n*n - 3.0*n);
//		HDrift(0,0) = Delta2-Delta;
//		HDrift = HDrift + Delta * Ident;
//		MOs::GenPauliX(UNOT,2,1);
//		sys.setlevel0index(1);
//		a=(1/sqrt(2))*a;
		
		sys.SetHdrift(HDrift);	
		HDrift.SetOutputStyle(Matrix);	
		
	//2nd transition
		
	//	a(0,1)=a(1,0)=0;
	//	a(3,4)=a(4,3)=0;
	//	a(2,3)=a(3,2)=sqrt(2)/2;
	//	lambda = sqrt(3)/sqrt(2);//abs(a(2,3));
	//	lamb2 = 1/sqrt(2);//abs(a(0,1));
	
		sys.SetRhoDesired(UNOT);	
		
		matrix<complex<double> > Shift(dim, dim);
		MOs::Null(Shift);
		matrix<complex<double> > Shift2(dim, dim);
		MOs::Null(Shift2);
			
		Shift(0,1) = Shift(1,0) =1;
		Shift2(0,1) =-ii; 
		Shift2(1,0) =ii;

		ad = MOs::Dagger(a);

		AnalyticControl refcontrol(0.5*(a+ad),rwa_dt, rwa_num_time,rwasub,NULL,"refctrl");
		refcontrol.setparams(&AnalyticControl::ShiftedGaussian, rab, 1.5, 1, pSetTransition, NULL);	
		refcontrol.pInterpolate = &AnalyticControl::Void;

//Gaussian
		sys.setFilename("gauss");
		sys.SetNumericalParameters(fidelity=0.00099999999999999999999,  base_a, epsilon, tolerance, max_iter);
		AnalyticControl u0(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&sys, "gaussian"), 
				u1(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub,&sys),
				u2(n,rwa_dt,rwa_num_time,rwasub,&sys), 
				u3(twophot,rwa_dt,rwa_num_time,rwasub,&sys);
	
		u0.setparams(&AnalyticControl::ShiftedGaussian, rab, 3, 1, pSetTransition, NULL);					
		u0.pInterpolate = u1.pInterpolate = u2.pInterpolate = u3.pInterpolate =  &AnalyticControl::Void;
		sys.sweeptimes();
//		u0.Hcontrol_(1,0) = u0.Hcontrol_(0,1) = u1.Hcontrol_(1,0) = u1.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0.Hcontrol_, u1.Hcontrol_, "lamb5nsHarm2.5_G3.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );

//GRAPE
		sys.setFilename("grape2chan_pix_1_d_ex.dat");
		sys.SetNumericalParameters(fidelity=0.9999999999999999999,  base_a, epsilon, tolerance, max_iter);
		AnalyticControl u0G(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&sys, "GRAPE"), 
				u1G(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub,&sys), 
				u2G(n,rwa_dt,rwa_num_time,rwasub,&sys),
				u3G(0.0*twophot,rwa_dt,rwa_num_time,rwasub,&sys);

		u0G.setparams(&AnalyticControl::ShiftedGaussian, rab, 1.5, 1, pSetTransition, &refcontrol); 
		u1G.setparams(&AnalyticControl::Derivative,  -1/Delta/2, 0, 0, pSetTransition, &refcontrol);

		u0G.pInterpolate = u1G.pInterpolate = u2G.pInterpolate = u3G.pInterpolate =  &AnalyticControl::Pixelate;//CubicInterpolate;//
		u0G.pSubGradients = u1G.pSubGradients = u2G.pSubGradients = u3G.pSubGradients = &AnalyticControl::setPixelsGradients;//setLinearSplineGradients;//

//		sys.sweeptimes(static_cast<ptrPropagate> (&OptimizeEvolution::UnitaryTransfer));
//		u0.Hcontrol_(1,0) = u0.Hcontrol_(0,1) = u1.Hcontrol_(1,0) = u1.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0.Hcontrol_, u1.Hcontrol_, "lamb5nsHarm2.5_G3.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );

		cout << "drag\n";
//1Drag
	sys.setFilename("timesSS");
	sys.SetNumericalParameters(fidelity=0.000000000000009999999999,base_a, epsilon, tolerance, max_iter);		
		AnalyticControl u0C(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&sys, "gauss1"), 
				u1C(0.0*0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub,&sys, "deriv1"), 
				u2C(n,rwa_dt,rwa_num_time,rwa_num_time,&sys, "delta1"),
				u3C(0.0*twophot,rwa_dt,rwa_num_time,rwasub,&sys);

		u0C.setparams(&AnalyticControl::GaussianDRAGA , rab, 1.5, 0, pSetTransition,  &refcontrol);
		u1C.setparams(&AnalyticControl::Derivative, -0/Delta/2, 0, 0, pSetTransition, &refcontrol); 
		u2C.setparams(&AnalyticControl::StarkShiftDRAGA, 0, 0, 0, pSetTransition, &refcontrol);

		u0C.pInterpolate = u1C.pInterpolate = u3C.pInterpolate =  &AnalyticControl::Void;//CubicInterpolate;//LinearInterpolate;//
		u2C.pInterpolate = &AnalyticControl::Pixelate;
		u0C.pSubGradients = u1C.pSubGradients = u3C.pSubGradients = &AnalyticControl::setAMGradient;//
		u2C.pSubGradients = &AnalyticControl::setPixelsGradients;//setLinearSplineGradients;//
		sys.sweeptimes((static_cast<ptrPropagate> (&OptimizeEvolution::UnitaryTransfer)));
//		u0C.Hcontrol_(1,0) = u0C.Hcontrol_(0,1) = u1C.Hcontrol_(1,0) = u1C.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0C.Hcontrol_, u1C.Hcontrol_, "lamb5nsHarm2.5_drgC_G2.dat");

//1DragB
	sys.setFilename("dragC2");
	sys.SetNumericalParameters(fidelity=0.9999999999000000000000009999999999,base_a, epsilon, tolerance, max_iter);		
		AnalyticControl u0B(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&sys, "gaussB"), 
				u1B(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub,&sys, "derivB"), 
				u2B(n,rwa_dt,rwa_num_time,rwa_num_time,&sys, "deltaB"),
				u3B(0.0*twophot,rwa_dt,rwa_num_time,rwasub,&sys);

		u0B.setparams(&AnalyticControl::ShiftedGaussian , rab, 1.5, 1, pSetTransition,  &refcontrol);
		u1B.setparams(&AnalyticControl::Derivative, -1/Delta, 0, 0, pSetTransition, &refcontrol); 
		u2B.setparams(&AnalyticControl::StarkShiftDRAGC, 0, 0, 0, pSetTransition, &refcontrol);

		u0B.pInterpolate = u1B.pInterpolate = u3B.pInterpolate =  &AnalyticControl::Void;//CubicInterpolate;//LinearInterpolate;//
		u2B.pInterpolate = &AnalyticControl::Pixelate;
		u0B.pSubGradients = u1B.pSubGradients = u3B.pSubGradients = &AnalyticControl::setAMGradient;//
		u2B.pSubGradients = &AnalyticControl::setPixelsGradients;//setLinearSplineGradients;//
	//	sys.sweeptimes((static_cast<ptrPropagate> (&OptimizeEvolution::UnitaryTransfer)));
//		u0C.Hcontrol_(1,0) = u0C.Hcontrol_(0,1) = u1C.Hcontrol_(1,0) = u1C.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0C.Hcontrol_, u1C.Hcontrol_, "lamb5nsHarm2.5_drgC_G2.dat");



//1DragA
	sys.setFilename("1DragA");
	sys.SetNumericalParameters(fidelity=0.00099999999999999999999, base_a, epsilon, tolerance, max_iter);
	const int nsub=30;	
		AnalyticControl u0A(0.5*(a+ad), rwa_dt/nsub, rwa_num_time*nsub, nsub,&sys), 
				u1A(0.0*(a+ad), rwa_dt/nsub, rwa_num_time*nsub, nsub,&sys), 
				u2A(0.0*n, rwa_dt/nsub, rwa_num_time*nsub, nsub,&sys),
				u3A(0.25*n, rwa_dt/nsub, rwa_num_time*nsub, nsub,&sys);	

		u0A.setparams(&AnalyticControl::ShiftedGaussian, rab, 3, 0, pSetTransition, &refcontrol);
		u1A.setparams(&AnalyticControl::Derivative, -1/Delta/2, 0, 0, pSetTransition, &refcontrol); 
		u2A.setparams(&AnalyticControl::StarkShiftDRAGAFreq, 0, 0, 0, pSetTransition, &refcontrol);
		u3A.setparams(&AnalyticControl::StarkShiftDRAGA, 0, 0, 0, pSetTransition, &refcontrol);		

		u1A.relphase_ = M_PI/2;
	    u0A.drivefreq_ = u1A.drivefreq_ = qbfreq; 
		u0A.framefreq_ = u1A.framefreq_ = qbfreq; 
		u0A.detctrl = u1A.detctrl = &u2A;
		u0A.pSetHcontrol = u1A.pSetHcontrol = &AnalyticControl::fmdrivenControl;	

		u0A.pInterpolate = u1A.pInterpolate = u2A.pInterpolate = u3A.pInterpolate =  &AnalyticControl::Void;
//		sys.sweeptimes();
//		u0A2.Hcontrol_(1,0) = u0A.Hcontrol_(0,1) = u1A2.Hcontrol_(1,0) = u1A.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0A.Hcontrol_, u1A.Hcontrol_, "lamb5nsHarm2.5_drgA_G3.dat");
u0A.writefile("gaussiandet.dat");
u1A.writefile("gaussiandet1.dat");
u2A.writefile("gaussiandet2.dat");
u3A.writefile("gaussiandet3.dat");

//1DragO3
		sys.setFilename("1DragO3");
		sys.SetNumericalParameters(fidelity=0.00099999999999999999999, base_a, epsilon, tolerance, max_iter);
		AnalyticControl u05(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&sys), 
				u15(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub,&sys), 
				u25(n,rwa_dt,rwa_num_time,rwasub,&sys),
				u35(twophot,rwa_dt,rwa_num_time,rwasub,&sys);

		u05.setparams(&AnalyticControl::Gaussian5DRAGC, rab, 1.5, 0, pSetTransition, &refcontrol); 
		u15.setparams(&AnalyticControl::Derivative7DRAGC, 0, 0, 0, pSetTransition, &refcontrol); 
		u25.setparams(&AnalyticControl::StarkShift4DRAGC, 0, 0, 0, pSetTransition, &refcontrol);
				
		u05.pInterpolate = u15.pInterpolate = u25.pInterpolate = u35.pInterpolate =  &AnalyticControl::Void;
		sys.sweeptimes();
//		u0C.Hcontrol_(1,0) = u0C.Hcontrol_(0,1) = u1C.Hcontrol_(1,0) = u1C.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0C.Hcontrol_, u1C.Hcontrol_, "lamb5nsHarm2.5_drgC_G2.dat");

//1DragTrig
		sys.setFilename("1DragTrig.dat");
		sys.SetNumericalParameters(fidelity=0.00099999999999999999999, base_a, epsilon, tolerance, max_iter);
		AnalyticControl u0T(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&sys,"gausstrig"), 
				u1T(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub,&sys,"dertrig"), 
				u2T(n,rwa_dt,rwa_num_time,rwasub,&sys,"detrig"),
				u3T(twophot,rwa_dt,rwa_num_time,rwasub,&sys,"twoptrig");
				
		u0T.setparams(&AnalyticControl::GaussianTrig, rab, 1.5, 0, pSetTransition, &refcontrol); 
		u1T.setparams(&AnalyticControl::DerivativeTrig, 0, 0, 0, pSetTransition, &refcontrol); 
		u2T.setparams(&AnalyticControl::StarkShiftTrig, 0, 0, 0, pSetTransition, &refcontrol);
		u3T.setparams(&AnalyticControl::TwophotonTrig, 0, 0, 0, pSetTransition, &refcontrol);
				
		u0T.pInterpolate = u1T.pInterpolate = u2T.pInterpolate = u3T.pInterpolate =  &AnalyticControl::Void;
//		sys.sweeptimes(&Evolution::forwardpropagate);
//		u0C.Hcontrol_(1,0) = u0C.Hcontrol_(0,1) = u1C.Hcontrol_(1,0) = u1C.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0C.Hcontrol_, u1C.Hcontrol_, "lamb5nsHarm2.5_drgC_G2.dat");
u0T.writefile("gaussiantrig.dat");
u1T.writefile("gaussiantrig1.dat");
u2T.writefile("gaussiantrig2.dat");
u3T.writefile("gaussiantrig3.dat");		
		
		//1Drag2Phot
		sys.setFilename("1Drag2Phot");
		sys.SetNumericalParameters(fidelity=0.99999999999999999999, base_a, epsilon, tolerance, max_iter);
		AnalyticControl u0P(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&sys, "gaussP"), 
				u1P(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub,&sys, "derivP"), 
				u2P(n,rwa_dt,rwa_num_time,rwasub,&sys, "deltaP"),
				u3P(twophot,rwa_dt,rwa_num_time,rwasub,&sys,"twophP");
								
		u0P.setparams(&AnalyticControl::ShiftedGaussian, rab, 1.5, 1, pSetTransition, &refcontrol); 
		u1P.setparams(&AnalyticControl::Derivative,  -1/Delta, 0, 0, pSetTransition, &refcontrol); 
		u2P.setparams(&AnalyticControl::StarkShiftDRAGC, 0, 0, 0, pSetTransition, &refcontrol);
		u3P.setparams(&AnalyticControl::ShiftedGaussian, rab/20, 1.5, 0, pSetTransition, &refcontrol);
				
				
		u0P.pInterpolate = u1P.pInterpolate = u3P.pInterpolate =  &AnalyticControl::Void;//CubicInterpolate;//LinearInterpolate;//
		u2P.pInterpolate = &AnalyticControl::Pixelate;
		u0P.pSubGradients = u1P.pSubGradients = u3P.pSubGradients = &AnalyticControl::setAMGradient;//
		u2P.pSubGradients = &AnalyticControl::setPixelsGradients;//setLinearSplineGradients;//
		cout << u3P.Hcontrol_ <<endl;	
		sys.sweeptimes((static_cast<ptrPropagate> (&OptimizeEvolution::UnitaryTransfer)));
//		u0C.Hcontrol_(1,0) = u0C.Hcontrol_(0,1) = u1C.Hcontrol_(1,0) = u1C.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0C.Hcontrol_, u1C.Hcontrol_, "lamb5nsHarm2.5_drgC_G2.dat");
		

		//2Drag
		sys.setFilename("2Drag.dat");
		sys.SetNumericalParameters(fidelity=0.00099999999999999999999,base_a, epsilon, tolerance, max_iter);
		AnalyticControl u0C2(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&sys), 
				u1C2(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub,&sys), 
				u2C2(n,rwa_dt,rwa_num_time,rwasub,&sys),
				u3C2(0.0*twophot,rwa_dt,rwa_num_time,rwasub,&sys);
				
		u0C2.setparams(&AnalyticControl::GaussianDRAGC2, rab, 2, 0, pSetTransition, &refcontrol); 
		u1C2.setparams(&AnalyticControl::Derivative, -1/Delta, 2, 0, pSetTransition, &refcontrol); 
		u2C2.setparams(&AnalyticControl::StarkShiftDRAGC2, 0, 0, 0, pSetTransition, &refcontrol);		
				
		u0C2.pInterpolate = u1C2.pInterpolate = u2C2.pInterpolate = u3C2.pInterpolate =  &AnalyticControl::Void;
//		sys.sweeptimes( "dragC2chan.dat");
//		u0C2.Hcontrol_(1,0) = u0C2.Hcontrol_(0,1) = u1C2.Hcontrol_(1,0) = u1C2.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0C2.Hcontrol_, u1C2.Hcontrol_, "lamb5nsHarm2.5_drgC2_G2.dat");

		
		//2Adi
		sys.setFilename("2Adi.dat");
		sys.SetNumericalParameters(fidelity=0.00099999999999999999999, base_a, epsilon, tolerance, max_iter);
		AnalyticControl u0A2(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&sys), 
				u1A2(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub,&sys), 
				u2A2(n,rwa_dt,rwa_num_time,rwasub,&sys),
				u3A2(0.0*twophot,rwa_dt,rwa_num_time,rwasub,&sys);
				
		u0A2.setparams( &AnalyticControl::GaussianDRAGA2, rab, 3, 0, pSetTransition, &refcontrol); 
		u1A2.setparams(&AnalyticControl::Derivative, 0, 0, 0, pSetTransition, &refcontrol);
		u2A2.setparams(&AnalyticControl::StarkShiftDRAGA2, 0, 0, 0, pSetTransition, &refcontrol);		
				
//		sys.sweeptimes( "what3.dat");
//		u0A2.Hcontrol_(1,0) = u0A2.Hcontrol_(0,1) = u1A2.Hcontrol_(1,0) = u1A2.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0A2.Hcontrol_, u1A2.Hcontrol_, "lamb5nsHarm2.5_drgA2_G3.dat");
		
		//2Drg/2
		sys.setFilename("2Drgdiv2.dat");
		sys.SetNumericalParameters(fidelity=0.00099999999999999999999, base_a, epsilon, tolerance, max_iter);
		AnalyticControl u0B2(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&sys), 
				u1B2(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub,&sys), 
				u2B2(n,rwa_dt,rwa_num_time,rwasub,&sys),
				u3B2(0.0*twophot,rwa_dt,rwa_num_time,rwasub,&sys);
				
		u0B2.setparams(&AnalyticControl::GaussianDRAGB2, rab, 3, 0, pSetTransition, &refcontrol); 
		u1B2.setparams(&AnalyticControl::DerivativeDRAGB2, 0, 0, 0, pSetTransition, &refcontrol);
				
//		sys.sweeptimes( "whathalf.dat");
//		u0B2.Hcontrol_(1,0) = u0B2.Hcontrol_(0,1) = u1B2.Hcontrol_(1,0) = u1B2.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0B2.Hcontrol_, u1B2.Hcontrol_, "lamb5nsHarm2.5_drgB2_G3.dat");
		
//		u1.Hcontrol_.SetOutputStyle(Matrix);
//		cout << u1.Hcontrol_ << endl;
		
//		u0.Hcontrol_.SetOutputStyle(Matrix);		
//		cout << u0.Hcontrol_ << endl;
//		u0.ShiftedGaussian(rab,tgate/1.5,&u0, NULL);
//		u0.Normalize(rab);
				
//		matrix<complex<double> > twoquantum(dim,dim);
//		MOs::Null(twoquantum);
//		twoquantum(1,3)=-ii;
//		twoquantum(3,1)=ii;
//		AnalyticControl u3(0.0*n,rwa_dt,rwa_num_time,rwasub, &AnalyticControl::StarkShift, 0*( -lamb2*lamb2/Delta2 +(lambda*lambda - 4)/Delta)/4, 2, 0);
		
//		sys.Setucontrol(&u3,3); //3nd contro	
//		u0.Hcontrol_(1,0) = u0.Hcontrol_(0,1) = 0;
//		Shift2(0,1) =-ii; 
//		Shift2(1,0) =ii;
//		u1.Hcontrol_(1,0) = u1.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0.Hcontrol_, u1.Hcontrol_, argv[1], static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
						
		
		//run grape
		clo= clock();	
	//	sys.sweeptimes( argv[1], static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
		
		//sweep delta
	//	HDrift(0,0)=0;
	//	Shift(0,0) = 1;
	//	sys.sweepenergies(Delta*1.1, Delta, 0.0*Shift, 0.0*Shift, HDrift,  HDrift, argv[1], static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );

		//sweep lambda
	
		//grape
	//	sys.UnitaryTransfer(static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) ); 
		cout << "time taken: " << to_ms*(clock() - clo) << endl;
	
	return 0;
	////////////////////////////
	
	
}