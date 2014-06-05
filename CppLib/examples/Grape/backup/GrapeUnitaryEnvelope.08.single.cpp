/*
Name: GrapeUnitaryEnvelope.cpp
Author: felix motzoi

Dependences: GrapeUnitaryEnvelope.hpp  
Brief Discription: This demonstrates grape for lab frame vs rotating frame with large pixels

Version History
	v0: Using Control class, Oct. 6 2009.

*/
#include "GrapeUnitaryQubit.hpp"

using namespace std;

int main (int argc, char const *argv[]){

	verbose=no;
	cout << "Running program " << argv[0] << endl;
	
	//System/Optimization parameters
	double qbfreq=51.5, Delta=-2.5, lambda=sqrt(3/2);  //physical params
	double lamb2 =  1/sqrt(2);
	double rab = M_PI;
	double Delta2 = Delta;
	double tolerance=0.00000000000000000000000000001, fidelity=0.9, base_a=1.5, epsilon=5, tgate=1, dt;  //search params
	size_t max_iter=100000;
	size_t dim=5;
	double to_ms =1000.0/CLOCKS_PER_SEC;
	int clo;
	
	//Typical operators	
	complex<double> ii(0.0,1.0);
	matrix<complex<double> > a(dim,dim), ad(dim,dim), n(dim,dim), UNOT(dim,dim), Ident(dim, dim),twophot(dim, dim);
	matrix<complex<double> > Urot(dim,dim);
	MOs::Destroy(a);
	ad = MOs::Dagger(a);
	n=ad*a;
	MOs::Identity(Ident);
	MOs::Null(twophot);
	twophot(0,2)=1;//-ii;
	twophot(2,0)=1;//ii;

	//RWA fixed parameters
	const size_t rwasub=10;
	size_t rwa_num_time=tgate*rwasub, num_controls =4;
	
	//RW fixed parameters
    //const int nsub=20;
	//size_t num_time=rwa_num_time*nsub/rwasub;
	//dt=tgate/double(num_time);
	
	system("mkdir hello");
	
		//RWA frame	
		double rwa_dt = tgate/double(rwa_num_time);
		GrapeUnitaryQubit sys(dim, num_controls, rwa_num_time);
		
		//sys.SetHdrift(Delta/2*(n*n-n));			//first transition
		matrix<complex<double> > HDrift(dim, dim);
		//level 0
		HDrift = Delta/2*( n*n - n);
		MOs::GenPauliX(UNOT,0,1);

//change lambda here
//		a(1,2)=a(2,1)=0;//1.75;
//		ad = MOs::Dagger(a);


	//level1
//		HDrift =  Delta/2*( n*n - 3.0*n);
//		HDrift(0,0) = Delta2-Delta;
//		HDrift = HDrift + Delta * Ident;
//		MOs::GenPauliX(UNOT,2,1);
//		sys.setlevel0index(1);
//		a=(1/sqrt(2))*a;

		
		sys.SetHdrift(HDrift);		//2nd transition
		
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
					
		sys.SetPhysicalParameters(fidelity=0.00099999999999999999999, rwa_dt, base_a, epsilon, tolerance, max_iter, qbfreq, Delta);
//Gaussian
		Control u0(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub, &Control::ShiftedGaussian, rab, 3, 1), 
				u1(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub, &Control::Null, 0, 0, 0), 
				u2(n,rwa_dt,rwa_num_time,rwasub, &Control::Null, 0, 0, 0),
				u3(twophot,rwa_dt,rwa_num_time,rwasub, &Control::Null, 0, 0, 0);
		sys.Setucontrol(&u0, 0); //0th control
		sys.Setucontrol(&u1,1); //1st control
		sys.Setucontrol(&u2,2); //2nd control
		sys.Setucontrol(&u3,3); //3rd control
//		sys.sweeptimes( "what2.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
//		u0.Hcontrol_(1,0) = u0.Hcontrol_(0,1) = u1.Hcontrol_(1,0) = u1.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0.Hcontrol_, u1.Hcontrol_, "lamb5nsHarm2.5_G3.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );


	sys.SetPhysicalParameters(fidelity=0.99999999, rwa_dt, base_a, epsilon, tolerance, max_iter, qbfreq, Delta);
//GRAPE
		Control u0G(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub, &Control::ShiftedGaussian, rab, 1.5, 1), 
				u1G(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub, &Control::Derivative,  -rab/Delta, 0, 0), 
				u2G(0.0*n,rwa_dt,rwa_num_time,rwasub, &Control::StarkShiftDRAGC, rab/Delta/Delta, 0, 0),
				u3G(0.0*twophot,rwa_dt,rwa_num_time,rwasub, &Control::Null, 0, 0, 0);
		sys.Setucontrol(&u0G, 0); //0th control
		sys.Setucontrol(&u1G,1); //1st control
		sys.Setucontrol(&u2G,2); //2nd control
		sys.Setucontrol(&u3G,3); //3rd control
		sys.sweeptimes( "grapeXYlinladder.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
//		u0.Hcontrol_(1,0) = u0.Hcontrol_(0,1) = u1.Hcontrol_(1,0) = u1.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0.Hcontrol_, u1.Hcontrol_, "lamb5nsHarm2.5_G3.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
u0G.writefile("ctrlgrape.dat");
u1G.writefile("ctrlgrape1.dat");
u2G.writefile("ctrlgrape2.dat");
u3G.writefile("ctrlgrape3.dat");

sys.SetPhysicalParameters(fidelity=0.00099999999999999999999, rwa_dt, base_a, epsilon, tolerance, max_iter, qbfreq, Delta);
		
//1Drag
		Control u0C(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub, &Control::GaussianDRAGC, rab, 1.5, 0), 
				u1C(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub, &Control::Derivative, -1/Delta, 0, 0), 
				u2C(n,rwa_dt,rwa_num_time,rwasub, &Control::StarkShiftDRAGC, 0, 0, 0),
				u3C(0.0*twophot,rwa_dt,rwa_num_time,rwasub, &Control::Null, 0, 0, 0);
		sys.Setucontrol(&u0C, 0); //0th control
		sys.Setucontrol(&u1C,1); //1st control
		sys.Setucontrol(&u2C,2); //2nd control
		sys.Setucontrol(&u3C,3); //3rd control
//		sys.sweeptimes( "whatnow2.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
//		u0C.Hcontrol_(1,0) = u0C.Hcontrol_(0,1) = u1C.Hcontrol_(1,0) = u1C.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0C.Hcontrol_, u1C.Hcontrol_, "lamb5nsHarm2.5_drgC_G2.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
u0C.writefile("gaussian.dat");
u1C.writefile("gaussian1.dat");
u2C.writefile("gaussian2.dat");
u3C.writefile("gaussian3.dat");


//1DragA
	const int nsub=30;	
	Control u0A(1.0*(a+ad), rwa_dt/nsub, rwa_num_time*nsub, nsub, &Control::ShiftedGaussian, rab, 3, 1), 
				u1A(0.0*(a+ad), rwa_dt/nsub, rwa_num_time*nsub, nsub, &Control::Derivative, -1/Delta/2, 0, 0), 
				u2A(0.0*n, rwa_dt/nsub, rwa_num_time*nsub, nsub, &Control::StarkShiftDRAGAFreq, 0, 0, 0),
				u3A(0.25*n, rwa_dt/nsub, rwa_num_time*nsub, nsub, &Control::StarkShiftDRAGA, 0, 0, 0);	
		u1A.relphase_ = M_PI/2;
	    u0A.drivefreq_ = u1A.drivefreq_ = qbfreq; 
		u0A.framefreq_ = u1A.framefreq_ = qbfreq; 
		u0A.detctrl = u1A.detctrl = &u2A;
		u0A.pSetHcontrol = u1A.pSetHcontrol = &Control::fmdrivenControl;	
		sys.Setucontrol(&u0A, 0); //0th control
		sys.Setucontrol(&u1A,1); //1st control
		sys.Setucontrol(&u2A,2); //2nd control
		sys.Setucontrol(&u3A,3); //3d control
		sys.SetPhysicalParameters(fidelity=0.00099999999999999999999, rwa_dt/nsub, base_a, epsilon, tolerance, max_iter, qbfreq, Delta);
//		sys.sweeptimes( "whata5.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
//		u0A2.Hcontrol_(1,0) = u0A.Hcontrol_(0,1) = u1A2.Hcontrol_(1,0) = u1A.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0A.Hcontrol_, u1A.Hcontrol_, "lamb5nsHarm2.5_drgA_G3.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
u0A.writefile("gaussiandet.dat");
u1A.writefile("gaussiandet1.dat");
u2A.writefile("gaussiandet2.dat");
u3A.writefile("gaussiandet3.dat");		

//1DragO3
		Control u05(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub, &Control::Gaussian5DRAGC, rab, 1.5, 0), 
				u15(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub, &Control::Derivative7DRAGC, 0, 0, 0), 
				u25(n,rwa_dt,rwa_num_time,rwasub, &Control::StarkShift4DRAGC, 0, 0, 0),
				u35(twophot,rwa_dt,rwa_num_time,rwasub, &Control::Null, 0, 0, 0);
		sys.Setucontrol(&u05, 0); //0th control
		sys.Setucontrol(&u15,1); //1st control
		sys.Setucontrol(&u25,2); //2nd control
		sys.Setucontrol(&u35,3); //3rd control
		sys.SetPhysicalParameters(fidelity=0.000099999999999999999999, rwa_dt, base_a, epsilon, tolerance, max_iter, qbfreq, Delta);
//		sys.sweeptimes( "lwhatnow2.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
//		u0C.Hcontrol_(1,0) = u0C.Hcontrol_(0,1) = u1C.Hcontrol_(1,0) = u1C.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0C.Hcontrol_, u1C.Hcontrol_, "lamb5nsHarm2.5_drgC_G2.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );

//1DragTrig
		Control u0T(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub, &Control::GaussianTrig, rab, 1.5, 0), 
				u1T(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub, &Control::DerivativeTrig, 0, 0, 0), 
				u2T(n,rwa_dt,rwa_num_time,rwasub, &Control::StarkShiftTrig, 0, 0, 0),
				u3T(twophot,rwa_dt,rwa_num_time,rwasub, &Control::TwophotonTrig, 0, 0, 0);
		sys.Setucontrol(&u0T, 0); //0th control
		sys.Setucontrol(&u1T,1); //1st control
		sys.Setucontrol(&u2T,2); //2nd control
		sys.Setucontrol(&u3T,3); //3rd control		
//		sys.sweeptimes( "whatrig2.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
//		u0C.Hcontrol_(1,0) = u0C.Hcontrol_(0,1) = u1C.Hcontrol_(1,0) = u1C.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0C.Hcontrol_, u1C.Hcontrol_, "lamb5nsHarm2.5_drgC_G2.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
u0T.writefile("gaussiantrig.dat");
u1T.writefile("gaussiantrig1.dat");
u2T.writefile("gaussiantrig2.dat");
u3T.writefile("gaussiantrig3.dat");		
		
		//1Drag2Phot
		Control u0P(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub, &Control::GaussianDRAGC, rab, 1.5, 0), 
				u1P(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub, &Control::Derivative,  -1/Delta, 0, 0), 
				u2P(n,rwa_dt,rwa_num_time,rwasub, &Control::StarkShiftDRAGC, 0, 0, 0),
				u3P(twophot,rwa_dt,rwa_num_time,rwasub, &Control::Twophoton, 0, 0, 0);
		sys.Setucontrol(&u0P, 0); //0th control
		sys.Setucontrol(&u1P,1); //1st control
		sys.Setucontrol(&u2P,2); //2nd control
		sys.Setucontrol(&u3P,3); //3rd control
//		sys.sweeptimes( "lwhatnext.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
//		u0C.Hcontrol_(1,0) = u0C.Hcontrol_(0,1) = u1C.Hcontrol_(1,0) = u1C.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0C.Hcontrol_, u1C.Hcontrol_, "lamb5nsHarm2.5_drgC_G2.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
		

	

		//2Drag
		Control u0C2(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub, &Control::GaussianDRAGC, rab, 2, 0), 
				u1C2(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub, &Control::Derivative, -1/Delta, 2, 0), 
				u2C2(n,rwa_dt,rwa_num_time,rwasub, &Control::StarkShiftDRAGC2, 0, 0, 0),
				u3C2(0.0*twophot,rwa_dt,rwa_num_time,rwasub, &Control::Null, 0, 0, 0);
		sys.Setucontrol(&u0C2, 0); //0th control
		sys.Setucontrol(&u1C2,1); //1st control
		sys.Setucontrol(&u2C2,2); //2nd control
		sys.Setucontrol(&u3C2,3); //3d control
//		sys.sweeptimes( "what2.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
//		u0C2.Hcontrol_(1,0) = u0C2.Hcontrol_(0,1) = u1C2.Hcontrol_(1,0) = u1C2.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0C2.Hcontrol_, u1C2.Hcontrol_, "lamb5nsHarm2.5_drgC2_G2.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );

		
		//2Adi
		Control u0A2(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub, &Control::GaussianDRAGA2, rab, 3, 0), 
				u1A2(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub, &Control::Derivative, 0, 0, 0), 
				u2A2(n,rwa_dt,rwa_num_time,rwasub, &Control::StarkShiftDRAGA2, 0, 0, 0),
				u3A2(0.0*twophot,rwa_dt,rwa_num_time,rwasub, &Control::Null, 0, 0, 0);
		sys.Setucontrol(&u0A2, 0); //0th control
		sys.Setucontrol(&u1A2,1); //1st control
		sys.Setucontrol(&u2A2,2); //2nd control
		sys.Setucontrol(&u3A2,3); //3d control
//		sys.sweeptimes( "what3.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
//		u0A2.Hcontrol_(1,0) = u0A2.Hcontrol_(0,1) = u1A2.Hcontrol_(1,0) = u1A2.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0A2.Hcontrol_, u1A2.Hcontrol_, "lamb5nsHarm2.5_drgA2_G3.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
		
		//2Drg/2
				Control u0B2(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub, &Control::GaussianDRAGB2, rab, 3, 0), 
				u1B2(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub, &Control::DerivativeDRAGB2, 0, 0, 0), 
				u2B2(n,rwa_dt,rwa_num_time,rwasub, &Control::StarkShift, 0, 0, 0),
				u3B2(0.0*twophot,rwa_dt,rwa_num_time,rwasub, &Control::Null, 0, 0, 0);
		sys.Setucontrol(&u0B2, 0); //0th control
		sys.Setucontrol(&u1B2,1); //1st control
		sys.Setucontrol(&u2B2,2); //2nd control
//		sys.sweeptimes( "what4.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
//		u0B2.Hcontrol_(1,0) = u0B2.Hcontrol_(0,1) = u1B2.Hcontrol_(1,0) = u1B2.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0B2.Hcontrol_, u1B2.Hcontrol_, "lamb5nsHarm2.5_drgB2_G3.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
		
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
//		Control u3(0.0*n,rwa_dt,rwa_num_time,rwasub, &Control::StarkShift, 0*( -lamb2*lamb2/Delta2 +(lambda*lambda - 4)/Delta)/4, 2, 0);
		
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