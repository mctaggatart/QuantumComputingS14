/*
Name: Grape.cpp
Author: Jay Gambetta

Dependences: GrapeUnitary.hpp  
Brief Discription: This demonstrates grape for unitary using Phi3 and doing a pi pulse on a qubit
Limitations: None

Version History
	v0: Feb  11th, 2009.

time ./GrapeUnitary1.e GrapeUnitary1.dat

*/
#include <OptimizeEvolution.hpp>
#include <random>
#include <iostream>
#include<cmath>
#include <math.h>
using namespace std;

int main (int argc, char const *argv[]){
	
	verbose=no;
	cout << "Running program " << argv[0] << endl;
	
	//Grape inputs
	size_t numDis=1, typeDis=2, num_time=100000, dim = 2, refdim=60, num_controls =1, max_iter=1000;//10000
	double tolerance=std::numeric_limits<double>::min(), fidelity=.999999, base_a=1.5, epsilon=5000, tgate=4, dt, alpha, gamma;
	dt=tgate/double(num_time);
	gamma=1.0/((1.0/dt)*11.0);
	cout<<gamma<<"gamma"<<endl;
	if(sqrt((double)(exp(-dt/gamma)))>.005){
	  cout<<"You have chosen an unsupported value for gamma and dt. Please make sure your Kraus Opperator is small enough for the function to show through the dissipation"<<endl;
	}
	//cout<<"GENDIS"<<dt<<endl;
	matrix<complex<double> >** dis;
	dis= new matrix<std::complex<double> >*[num_time];
	  for(size_t l=0; l<num_time; ++l){
	    dis[l]=new matrix<std::complex<double> > [typeDis*numDis];
	for(size_t k=0; k<typeDis*numDis; ++k){
	   dis[l][k].initialize(dim,dim); //NEED TO MAKE A A GLOBAL VARIABLE SET!!
	  }
}	
	  
	  default_random_engine randgen;
	  normal_distribution<double> normaldis(.002, .0005);	//sets to matrix of each size
	
	  /*	for(int k=0; k<num_time; ++k){
	  //	  MOs::Destroy(dis[k][0]);
	for(int i=0; i<dim; ++i){
	   dis[k][0](i,i)=1.0;
	   dis[k][0](i,i)=dis[k][0](i,i)*(normaldis(randgen));
	  for(int j=0; j<dim; ++j){
	    //	    dis[k][0](i,j)=dis[k][0](i,j)*(normaldis(randgen));
	    
	    }}}*/
	
	  for(int k=0; k<num_time; ++k){
	  //MOs::Destroy(dis[k][0]);
	    //dis[k][0]=(MOs::Dagger(dis[k][0]))*(dis[k][0]);
	    dis[k][0](1,1)=sqrt(exp((-dt/gamma)));
			cout<<dis[k][0](1,1)<<endl;
	for(int i=0; i<dim; ++i){
	  for(int j=0; j<dim; ++j){
	    //	    	    dis[k][0](i,j)=(dis[k][0](i,j))*(normaldis(randgen));

	  }}
	dis[k][1](0,0)=1.0;
	dis[k][1](1,1)=sqrt((1.0-exp(-dt/gamma)));
		cout<<dis[k][1](1,1)<<endl;
}

	//cout<<dis[1][1]<<"dis test"<<endl;
	OptimizeEvolution sys(dim, num_controls, num_time, dt, dis, numDis, typeDis, "Unitary1");
	sys.SetNumericalParameters(fidelity, base_a, epsilon, tolerance, max_iter);

	matrix<complex<double> >  a(dim,dim), ad(dim,dim), n(dim,dim), Hdrift(dim,dim), Hcontrol(dim,dim);
	double delta=0;
	MOs::Destroy(a); ad = MOs::Dagger(a);
	n=ad*a;
	Hdrift = delta*n;
	Hcontrol=0.5*(a+ad);
	sys.SetHdrift(Hdrift);
 	
	Control u0(Hcontrol, dt, num_time, 1, &sys, "u1_control");
	u0.ShiftedGaussian(M_PI, alpha=2.5, NULL);
	u0.Normalize(M_PI/10);
	//Initial condition
	matrix<complex<double> > U_desired(dim,dim);
	U_desired.SetOutputStyle(Matrix);
	MOs::Null(U_desired);
		U_desired(1,0)=std::complex<double>(0,-1);
		U_desired(0,1)=std::complex<double>(0,-1);
		// NOT phase preserving version of U Des. 
	//	U_desired(0,1)=std::complex<double>(1,0);
	//U_desired(1,0)=std::complex<double>(1,0);
	       
		//gamma=1/(omega/10) eg, for omega=50, gamma=1/(50/10)=1/5=.2
		//This means T1A=Omega/10, T2A=omega/20. this can be adjusted inside evolution. 
		sys.SetOmega(50.0);
		//rho initial is set inside Evolution, it should be set with a calling method before this.

matrix<complex<double> > init;
		init.initialize(dim,dim);
		init(0,0)=1;
		sys.SetOppsDesired(dis);
		sys.SetRhoInitial(init);
		sys.SetUDesired(U_desired);
	sys.SetTrueRhoDesired(U_desired);
	sys.Phi=&OptimizeEvolution::Phi2;
	sys.gradPhi=&OptimizeEvolution::GradPhi2;
	sys.pPropagate=&Evolution::forwardpropagate_Density;
	sys.pGetGradient=&OptimizeEvolution::computegradient_Density;
	//	sys.Phi=&OptimizeEvolution::robustfidelity;
	// 	sys.pPropagate=static_cast<ptrPropagate>(&OptimizeEvolution::robustforwardpropagate);
	//	sys.pGetGradient=&OptimizeEvolution::robustgradient;

	size_t num_evol = 9;
	OptimizeEvolution **evolsys= new OptimizeEvolution*[num_evol];

	sys.evolutions=evolsys;
	sys.num_evol=num_evol;
		for(size_t s=0; s<num_evol; s++){
	  evolsys[s]=new OptimizeEvolution(dim, num_controls, num_time, dt,dis, numDis, typeDis, "alt evol");
	evolsys[s]->SetNumericalParameters(fidelity=.999999, base_a, epsilon, tolerance, max_iter);
	evolsys[s]->SetOmega(50.0);
	evolsys[s]->SetOppsDesired(dis);
	evolsys[s]->SetUDesired(U_desired);
	evolsys[s]->SetTrueRhoDesired(U_desired);
	evolsys[s]->Phi=&OptimizeEvolution::Phi2;
	evolsys[s]->gradPhi=&OptimizeEvolution::GradPhi2;
	evolsys[s]->pPropagate=&Evolution::forwardpropagate_Density;
	evolsys[s]->pGetGradient=&OptimizeEvolution::computegradient_Density;

	evolsys[s]->SetHdrift(Hdrift);
	//	evolsys[s]->u_=&(sys.u_);
 	}

		
		cout << "duplicating controls...\n";
		AnalyticControl *sysControls[num_evol][num_controls];
	Control *refdimcon[num_controls];
		srand(time(0));
	
		
		char filenames[num_controls][80];
		for(int k=0; k<num_controls; k++)
		{	for (size_t s=0; s<num_evol; s++)
			{
				sysControls[s][k] = new AnalyticControl(*(AnalyticControl*)sys.GetuControl(k), 100);				
				sysControls[s][k]->setparams(&Control::ShiftedGaussian,M_PI,2,1,&AnalyticControl::setdriving01,(AnalyticControl*)sys.GetuControl(k));
 
					sysControls[s][k]->pSetHcontrol =  &Control::drivenControl;//CubicInterpolate;//Pixelate;//GaussianFilter;//Void;////LinearInterpolate;//
					sysControls[s][k]->pSetHgradient = &Control::drivenControlGradient;//setPixelsGradients;//				
				
					//		sysControls[s][k]->pSetHcontrol = &Control::linearControl;
					// 		sysControls[s][k]->pSetHgradient =  &Control::linearGradient;
					//	sysControls[s][k]->drivefreq_= M_PI;
					//	sysControls[s][k]->framefreq_=M_PI;
				
				sysControls[s][k]->evol = evolsys[s];
				if(s) strcpy(filenames[k],"sampled ");
				else  strcpy(filenames[k],"splited ");
				filenames[k][7] = '0'+k;
				strcpy(sysControls[s][k]->filename,filenames[k]);
				evolsys[s]->Setucontrol(sysControls[s][k]);
	
				//	if(!k) 
				sysControls[s][k]->Hcontrol_ = (a+ad);
				//for(int i=0; i<num_time; i++){		
				//	sysControls[s][k]->u_=(sys.theControls_[k]->u_);
	//	}
				//	
				//	if(k==1){  //sysControls[s][1]->Hcontrol_ = 0.5*(-ii*a+ii*ad);
					  //	     sysControls[s][1]->relphase_ = M_PI/2;
				//		}
					//			
					//			if(!k) {	sysControls[s][k]->phase_ = M_PI*(rand()/(double)RAND_MAX);
			//				cout << s << " " << sysControls[s][k]->phase_ << endl;
			//			}
			//	else 	
				//		sysControls[s][k]->phase_ = phases[s];	
			}
			
		  //		refdimcon[k] = new Control(*sys.GetuControl(k), 100);
		  //		refdimcon[k]->Hcontrol_ = sys.GetuControl(k)->Hcontrol_;//
		  //		refdimcon[k]->Hcontrol_.resize(60,60);
		  //		for(size_t d = 2; d<60; d++)
		  //	refdimcon[k]->Hcontrol_(1,d) = refdimcon[k]->Hcontrol_(d,1) = 1.0;		
		  		}
	
	//run grape	
	sys.UnitaryTransfer();
	
	
		
	return 0;
}
