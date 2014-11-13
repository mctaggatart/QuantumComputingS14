/*

EDITS UNTESTED SUMMER 2014, simply done to compile
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
	double qbfreq=2.5*2*M_PI, Delta=-0.5*2*M_PI, lambda=sqrt(3/2);  //physical params
	double lamb2 =  1/sqrt(2);
	double Delta2 = Delta;
	double tolerance=0.00000000000000001, fidelity=0.9, base_a=1.7, epsilon=0.000005, tgate=4 , dt;  //search params
	double rab = M_PI;
	size_t max_iter=500000;
	size_t dim=3, refdim=60;
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
	
	const size_t rwasub=320;
	
	
	size_t rwa_num_time= tgate*rwasub, num_controls =2;
	
	//RW fixed parameters
    //const int nsub=20;
	//size_t num_time=rwa_num_time*nsub/rwasub;
	//dt=tgate/double(num_time);
	
//	system("mkdir hello");
	
		//RWA frame	
		double rwa_dt = tgate/double(rwa_num_time);		
		size_t numDis=1, typeDis=2;
  matrix<complex<double> >* dis;
	dis= new matrix<std::complex<double> >[numDis*typeDis];
	for(int k=0; k<typeDis*numDis; ++k){
	  dis[k].initialize(dim,dim); //NEED TO MAKE A A GLOBAL VARIABLE SET!!
}	
		MOs::Destroy(dis[0]);
	dis[1]=(MOs::Dagger(dis[0]))*(dis[0]);
 

	OptimizeEvolution baselinesys(dim, num_controls, rwa_num_time, rwa_dt,dis, numDis, typeDis, "EnvelopeCompare_test");
		baselinesys.Phi = &OptimizeEvolution::Phi4Sub2;
		baselinesys.gradPhi = &OptimizeEvolution::GradPhi4Sub2;
		baselinesys.pPropagate = &Evolution::forwardpropagate;
		baselinesys.pGetGradient = 	&OptimizeEvolution::computegradient;
		baselinesys.SetOppsDesired(dis);
		baselinesys.gatesteptime = 1;
		//baselinesys.ngatetimes = 6;  //times
		baselinesys.ngatetimes = 30;  //subpix
		baselinesys.startfid = 0.99;
		
		matrix<complex<double> > HDriftRef(dim, dim);
		//level 0
		HDriftRef = Delta/2*( n*n - n);
		HDriftRef.SetOutputStyle(Matrix);	
		HDriftRef.resize(refdim, refdim);
		for(size_t d = 3; d<refdim; d++)
			HDriftRef(d,d) = HDriftRef(2,2)*(1.0+sqrt(d)/10.0);
		
		matrix<complex<double> > HDrift(dim, dim);
		//level 0
		HDrift = Delta/2*( n*n - n);		
		
		ptrSetTransition pSetTransition = &AnalyticControl::setdriving01;
	
		HDriftRef.SetOutputStyle(Matrix);	
		//cout << HDriftRef << endl;
		baselinesys.SetHdrift(HDrift);

		baselinesys.SetNumericalParameters(fidelity=0.999999999999999999, base_a, epsilon, tolerance, max_iter);	
		//	baselinesys.SetRhoDesired(UNOT);
		matrix<complex<double> > init;
		init.initialize(dim,dim);
		init(1,1)=1;
baselinesys.SetRhoInitial(init);
	baselinesys.SetUDesired(UNOT);
	baselinesys.SetTrueRhoDesired(UNOT);
	baselinesys.SetOmega(5.0);
		cout << "declaring controls...\n";
		
		ad = MOs::Dagger(a);
		AnalyticControl refcontrol(0.5*(a+ad), rwa_dt, rwa_num_time, rwasub, NULL, "refcontrol");
		refcontrol.setparams(&AnalyticControl::ShiftedGaussian, rab, 2, 1);
		refcontrol.pInterpolate =  &AnalyticControl::Void;

		
		AnalyticControl u0C(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&baselinesys, "pixcontrol0"); 
		AnalyticControl	u1C(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub,&baselinesys, "pixcontrol1");
	//	AnalyticControl	u2C(n,rwa_dt,rwa_num_time,rwasub,&baselinesys, "pixcontrol2");
	//	AnalyticControl	u3C(twophot,rwa_dt,rwa_num_time,rwasub,&baselinesys, "pixcontrol3");

		u0C.setparams(&AnalyticControl::ShiftedGaussian, rab, 2, 1, pSetTransition, NULL);
//		u0C.setparams(&AnalyticControl::SquarePulse, rab, -rab, 0, pSetTransition, &refcontrol);
		u1C.setparams(&AnalyticControl::Derivative, -1/Delta/2, -rab, 0, pSetTransition, &refcontrol); 
	//	u2C.setparams(&AnalyticControl::Null, rab, -rab, 0, pSetTransition, &refcontrol);
	//	u3C.setparams(&AnalyticControl::Null, rab, -rab, 0, pSetTransition, &refcontrol); 

	//	u3C.pInterpolate = u2C.pInterpolate = 
		u1C.pInterpolate = u0C.pInterpolate =  &AnalyticControl::GaussianFilter; ////;//;//Void;////LinearInterpolate;//
		u1C.pSubGradients = u0C.pSubGradients = &AnalyticControl::GaussianFilterGradient; //;

		
		cout << "creating alternate evolutions...\n";
		size_t num_evol=5;
		OptimizeEvolution **sys = new OptimizeEvolution*[num_evol];
		baselinesys.evolutions = sys;	
		baselinesys.num_evol = num_evol;
		size_t mult[] = {320,80, 320, 1, 320}; //sweeppixels
		//size_t mult[] = {160, 80, 160, 16, 320}; //sweeptime
		
		for (size_t s=0; s<num_evol; s++)
		{	
		
			sys[s] = new OptimizeEvolution(dim, num_controls, rwa_num_time/mult[s], rwa_dt*mult[s],  dis, numDis, typeDis, "altevol");
			sys[s]->SetNumericalParameters(fidelity=0.999999999999999999, base_a, epsilon, tolerance, max_iter);
			sys[s]->Phi = &OptimizeEvolution::Phi4Sub2;
			sys[s]->gradPhi = &OptimizeEvolution::GradPhi4Sub2;
			sys[s]->pPropagate = &Evolution::forwardpropagate;
			sys[s]->pGetGradient = 	&OptimizeEvolution::computegradient;
			sys[s]->SetHdrift(HDrift);						
		}
		
	//	sys[1]->pPropagate = &Evolution::forwardpropagate_Commute;	
	//	sys[1]->pGetGradient = 	&OptimizeEvolution::computegradient_Commute;			
		
		cout << "duplicating controls...\n";
		AnalyticControl *sysControls[num_evol][num_controls];
		Control *refdimcon[num_controls];
		
		char filenames[num_controls][80];
		for(int k=0; k<num_controls; k++)
		{	for (size_t s=0; s<num_evol; s++)
			{
				sysControls[s][k] = new AnalyticControl(*(AnalyticControl*)baselinesys.GetuControl(k), rwasub/mult[s]);				
				sysControls[s][k]->setparams(&AnalyticControl::ShiftedGaussian,rab,2,1,pSetTransition,(AnalyticControl*)baselinesys.GetuControl(k));
				
				sysControls[s][k]->evol = sys[s];
				if(s==4){ strcpy(filenames[k],"pixels ");
				sysControls[s][k]->pInterpolate = &AnalyticControl::Pixelate; //:LinearInterpolate;
				sysControls[s][k]->pSubGradients = &AnalyticControl::setPixelsGradients; //setLinearSplineGradients;
				}
				else if(s==0 || s==1){  if(s==0) strcpy(filenames[k],"spline "); else strcpy(filenames[k],"splint ");
				sysControls[s][k]->pInterpolate = &AnalyticControl::CubicInterpolate;//GaussianFilter;//:LinearInterpolate;
				sysControls[s][k]->pSubGradients = &AnalyticControl::setCubicSplineGradients;//setLinearSplineGradients;//GaussianFilterGradient;//;				
				}
				else{  strcpy(filenames[k],"filter ");
				sysControls[s][k]->pInterpolate = &AnalyticControl::GaussianFilter;//;//:LinearInterpolate;
				sysControls[s][k]->pSubGradients = &AnalyticControl::GaussianFilterGradient;//;//;				
				}
				filenames[k][6] = '0'+k;
				strcpy(sysControls[s][k]->filename,filenames[k]);
				sys[s]->Setucontrol(sysControls[s][k]);
			}
			
			refdimcon[k] = new Control(*baselinesys.GetuControl(k), rwasub/mult[0]);
			refdimcon[k]->Hcontrol_ = baselinesys.GetuControl(k)->Hcontrol_;
			refdimcon[k]->Hcontrol_.resize(refdim,refdim);
			for(size_t d = 3; d<refdim; d++)
				refdimcon[k]->Hcontrol_(1,d) = refdimcon[k]->Hcontrol_(d,1) = 1.0;		
		}
	
		cout << "beginning sweep...\n";
//		baselinesys.sweeptimesandcompare(static_cast<ptrPropagate> (&OptimizeEvolution::UnitaryTransfer));
		baselinesys.UnitaryTransfer();
//		baselinesys.sweepdimandcompare( "com_evolrand_frq10_Hd5_Hc3.dat", HDriftRef, refdimcon);

		cout <<"deleting controls\n";		
		for (size_t s=0; s<num_evol; s++)
		for(int k=0; k<num_controls; k++)
			delete sysControls[s][k];
		for(int k=0; k<num_controls; k++)
			delete refdimcon[k];

		cout <<"deleting alternate evolutions\n";		
		for (size_t s=0; s<num_evol; s++) delete sys[s];
		delete [] sys;
		
		cout <<"done\n";		
		return 0;
	////////////////////////////
	
	
}
