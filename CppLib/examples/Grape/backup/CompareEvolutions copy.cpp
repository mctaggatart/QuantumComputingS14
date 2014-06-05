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
	srand ( time(NULL) );
	
	//System/Optimization parameters
	double qbfreq=51.5, Delta=-2.75, lambda=sqrt(3/2);  //physical params
	double lamb2 =  1/sqrt(2);
	double Delta2 = Delta;
	double tolerance=0.00000000000000001, fidelity=0.9, base_a=1.7, epsilon=0.000005, tgate=5, dt;  //search params
	double rab = M_PI/tgate;
	size_t max_iter=500000;
	size_t dim=4;
	double to_ms =1000.0/CLOCKS_PER_SEC;
	
	
	//Typical operators	
	complex<double> ii(0.0,1.0);
	matrix<complex<double> > a(dim,dim), ad(dim,dim), n(dim,dim), Ident(dim, dim), twophot(dim, dim), UNOT(dim,dim) ;
	MOs::Destroy(a);
	ad = MOs::Dagger(a);
	n=ad*a;
	MOs::Identity(Ident);
	MOs::GenPauliX(UNOT,0,1);
	MOs::Null(twophot);
	twophot(0,2)=1;//-ii;
	twophot(2,0)=1;//ii;

	//RWA fixed parameters
	
	const size_t rwasub=500;
	size_t mult = 250;
	
	size_t rwa_num_time= tgate*rwasub, num_controls =2;
	
	//RW fixed parameters
    //const int nsub=20;
	//size_t num_time=rwa_num_time*nsub/rwasub;
	//dt=tgate/double(num_time);
	
//	system("mkdir hello");
	
		//RWA frame	
		double rwa_dt = tgate/double(rwa_num_time);		
		
		OptimizeEvolution sys(dim, num_controls, rwa_num_time/mult, rwa_dt*mult), baselinesys(dim, num_controls, rwa_num_time, rwa_dt);
		
		baselinesys.Phi = sys.Phi = &OptimizeEvolution::Phi4Sub2;
		baselinesys.gradPhi = sys.gradPhi = &OptimizeEvolution::GradPhi4Sub2;
		
		//sys.SetHdrift(Delta/2*(n*n-n));			//first transition
		matrix<complex<double> > HDrift(dim, dim);
		//level 0
		HDrift = Delta/2*( n*n - n);
		//HDrift(0,0) = 5;
	
		ptrSetTransition pSetTransition = &AnalyticControl::setdriving01;
	
		HDrift.SetOutputStyle(Matrix);	
		cout << HDrift << endl;
		baselinesys.SetHdrift(HDrift);
		sys.SetHdrift(HDrift);	
		
		ad = MOs::Dagger(a);

		cout << "declaring controls\n";
		
		AnalyticControl refcontrol(0.5*(a+ad), rwa_dt, rwa_num_time, rwasub, NULL, "refcontrol");
		refcontrol.setparams(&AnalyticControl::ShiftedGaussian, rab, 1.5, 1);
		refcontrol.pInterpolate =  &AnalyticControl::Void;

		sys.SetNumericalParameters(fidelity=0.00000099999, base_a, epsilon, tolerance, max_iter);		
		baselinesys.SetNumericalParameters(fidelity=0.99999000000099900099, base_a, epsilon, tolerance, max_iter);	
		baselinesys.SetRhoDesired(UNOT);

		sys.pPropagate = &Evolution::forwardpropagate;

		sys.evolutions = &baselinesys;	

		AnalyticControl u0C(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&baselinesys, "truecontrol0"); 
		AnalyticControl	u1C(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub,&baselinesys, "truecontrol1");
	//	AnalyticControl		u2C(n,rwa_dt,rwa_num_time,rwasub,&baselinesys, "truecontrol2");
	//			u3C(0.0*twophot,rwa_dt,rwa_num_time,rwasub,&baselinesys, "truecontrol3");

		AnalyticControl *sysControls[num_controls];

		u0C.setparams(&AnalyticControl::ShiftedGaussian, rab, 2, 1, pSetTransition, &refcontrol);
//		u0C.setparams(&AnalyticControl::SquarePulse, rab, -rab, 0, pSetTransition, &refcontrol);
		u1C.setparams(&AnalyticControl::Null, rab, -rab, 0, pSetTransition, &refcontrol); 
	//	u2C.setparams(&AnalyticControl::Null, rab, -rab, 0, pSetTransition, &refcontrol);
	//	u3C.setparams(&AnalyticControl::Null, rab, -rab, 0, pSetTransition, &refcontrol); 

	//	u3C.pInterpolate = u2C.pInterpolate = 
		u1C.pInterpolate = u0C.pInterpolate =  &AnalyticControl::CubicInterpolate;//LinearInterpolate;//Pixelate;//GaussianFilter;//Void;////LinearInterpolate;//
		u1C.pSubGradients = u0C.pSubGradients = &AnalyticControl::setLinearSplineGradients;//setPixelsGradients;//
		
		cout << "duplicating controls\n";
		char filenames[num_controls][80];
		for(int k=0; k<num_controls; k++)
		{	sysControls[k] = new AnalyticControl(*(AnalyticControl*)baselinesys.GetuControl(k), rwasub/mult);
			sysControls[k]->setparams(&AnalyticControl::Replicate,0,0,0,pSetTransition,(AnalyticControl*)baselinesys.GetuControl(k));
			
			sysControls[k]->pInterpolate = &AnalyticControl::Void; //:LinearInterpolate;
			sysControls[k]->pSubGradients = &AnalyticControl::setSubVoidGradients; //setLinearSplineGradients;
			
			sysControls[k]->evol = &sys;
			strcpy(filenames[k],"sampled ");
			filenames[k][7] = '0'+k;
			sysControls[k]->filename =filenames[k];
			sys.Setucontrol(sysControls[k]);
		}
	
		cout << "beginning sweep\n";
//		sys.sweeptimesandcompare( "timecslice40.dat");
		sys.sweepfinegrainingandcompare( "com_evolrand_frq10_Hd5_Hc3.dat");

		for(int k=0; k<num_controls; k++)
			delete sysControls[k];

			return 0;
	////////////////////////////
	
	
}