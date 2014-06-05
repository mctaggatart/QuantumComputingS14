/*
Name: GrapeUnitaryEnvelope.cpp
Author: felix motzoi

Dependences: GrapeUnitaryEnvelope.hpp  
Brief Discription: This demonstrates grape for lab frame vs rotating frame with large pixels

Version History
	v0: Using Control class, Oct. 6 2009.

*/
#include "OptimizeEvolution.hpp"
using namespace std;

int main (int argc, char const *argv[]){

	verbose=no;
	cout << "Running program " << argv[0] << endl;

	//System/Optimization parameters
	double qbfreq=2*M_PI*4, Delta=-2.03, lambda=sqrt(2);  //physical params
	double tolerance=0.00000000000000000000001, fidelity=0.9, base_a=1.5, epsilon=5, tgate=0.25, dt;  //search params
	size_t max_iter=20000;
	size_t dim=2;
	double to_ms =1000.0/CLOCKS_PER_SEC;
	int clo;

	//Typical operators	
	complex<double> ii(0.0,1.0);
	matrix<complex<double> > a(dim,dim), ad(dim,dim), n(dim,dim), UNOT(dim,dim);
	matrix<complex<double> > Urot(dim,dim);
	matrix<complex<double> >HDrift(dim,dim);
	MOs::GenPauliX(UNOT,0,1);
	MOs::Destroy(a);
	//a(2,3)=a(3,2)=0;
	ad = MOs::Dagger(a);
	n=ad*a;

	//Grape inputs
	string const outfile = argv[1];
	ofstream dataout;
	UFs::OpenFile(outfile,dataout, 16);

	//RWA fixed parameters
    const int npert=8;
	const size_t rwasub=1;
	size_t rwa_num_time=tgate*rwasub*npert, num_controls =3;

	//RW fixed parameters
    const int nsub=150;
	size_t num_time=rwa_num_time*nsub/rwasub;
	dt=tgate/double(num_time);

	
//	for(size_t itime=0; itime<20; itime++, tgate+=0.25, num_time+=nsub*npert/4, rwa_num_time+=0.25*npert)
	{	
		//RWA frame	
		double rwa_dt = tgate/double(rwa_num_time);
		OptimizeEvolution sys(dim, num_controls, rwa_num_time, rwa_dt);
		sys.SetNumericalParameters(fidelity=0.00009999999, base_a, epsilon, tolerance, max_iter);

		sys.SetHdrift(HDrift=0.0*(n));
		
		AnalyticControl u0(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub,&sys, "labcontrol0"),
						u1(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time, rwasub, &sys, "labcontrol1"),
						u2(n,rwa_dt,rwa_num_time,rwasub, &sys, "labcontrol2");

		u0.setparams(&AnalyticControl::ShiftedGaussian, M_PI, 1.5, 1);
		u1.setparams(&AnalyticControl::Derivative,  1/qbfreq/2, 0, 0, NULL, &u0);
		u2.setparams(&AnalyticControl::StarkShift,  1.0/4.0/qbfreq, 0, 0, NULL, &u0);

		u0.init(); 
		u1.init(); 
		u2.init(); 

		//double norm = u0.Normalize(M_PI);
		//u1.GaussianDerivative(0, tgate, M_PI/4/qbfreq, tgate/4);
		for(size_t j = 0; j < rwa_num_time; ++j)
		{ 	//u2.u_[j]=-3.42657/(tgate*tgate*qbfreq); 
//			u2.u_[j]= (u0.u_[j]*u0.u_[j])/4/qbfreq;
			u0.u_[j]-= u0.u_[j]*u0.u_[j]*u0.u_[j]/qbfreq/qbfreq/16;  
		}		

		u0.writefile("RWA0.dat");	
		u1.writefile("RWA1.dat");	
		u2.writefile("RWA2.dat");	
		sys.Setucontrol(&u0,0); //0th control
		sys.Setucontrol(&u1,1); //1st control
		sys.Setucontrol(&u2,2); //2nd control

		sys.SetRhoDesired(UNOT);
		//run grape
		clo= clock();	

	//	sys.UnitaryTransfer( static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2),  static_cast< GradFid >(&GrapeUnitaryQubit::GradPhi4Sub2)); 

		dt=tgate/double(num_time);	

		AnalyticControl uR0(u0,nsub), uR1(u1,nsub), uR2(u2,nsub);
		matrix<complex<double> > Hc = a + ad;
		uR0.Hcontrol_ = uR1.Hcontrol_ = Hc;
	    uR0.drivefreq_ = uR1.drivefreq_ = qbfreq; 
		uR0.framefreq_ = uR1.framefreq_ = qbfreq; 
		uR0.pSetHcontrol = uR1.pSetHcontrol = &Control::drivenControl;		
		uR0.pSetHgradient = uR1.pSetHgradient = &Control::drivenControlGradient;
		uR0.pInterpolate = uR1.pInterpolate = uR2.pInterpolate =  &Control::LinearInterpolate;//Pixelate;//
		uR0.pSubGradients = uR1.pSubGradients = uR2.pSubGradients =  &Control::setLinearSplineGradients;//setPixelsGradients;//
		uR1.relphase_ = M_PI/2;
		
		for(size_t j = 0; j < rwa_num_time; ++j)
		{	//u2.u_[j]=-3.42657/(tgate*tgate*qbfreq); 
			//u0.u_[j]=u0.u_[j]*cos(qbfreq*j*rwa_dt)+u1.u_[j]*sin(qbfreq*j*rwa_dt);
			//u1.u_[j]*= sin(qbfreq*j*rwa_dt);
		}		
		
		uR0.writefile("RW0.dat");
		OptimizeEvolution sysRWC(dim, num_controls, num_time, dt);	
		sysRWC.SetNumericalParameters(fidelity=0.001,  base_a, epsilon, tolerance, max_iter);
		sysRWC.Phi = &OptimizeEvolution::Phi4Sub2;
		sysRWC.gradPhi = &OptimizeEvolution::GradPhi4Sub2;
		sysRWC.SetRhoDesired(UNOT);	
		sysRWC.nconfigs_=20;		
		sysRWC.SetHdrift(HDrift=0.0*(n));
		uR0.config_ =uR1.config_ = &sysRWC.config_;
		uR0.nconfigs_ = uR1.nconfigs_ = &sysRWC.nconfigs_;
						//
//		uR0.settotaltime(num_time);
//		uR1.settotaltime(num_time);
//		uR2.settotaltime(num_time);
//	//	
//		uR0.ShiftedGaussian( 0, tgate, M_PI, tgate/2);
//		norm = uR0.Normalize(M_PI);
//		for(size_t j = 0; j < num_time; ++j)
//			uR2.u_[j]-=uR0.u_[j]*uR0.u_[j]*cos(2*qbfreq*(j+0.5)*dt)/2/qbfreq;
//				
		
		uR0.init(); 
		uR1.init(); 
		uR2.init(); 		
		for(size_t j = 0; j < num_time; ++j)
		{ 	//u2.u_[j]=-3.42657/(tgate*tgate*qbfreq); 
//			u2.u_[j]= (u0.u_[j]*u0.u_[j])/4/qbfreq;
			uR0.u_[j]-= uR0.u_[j]*uR0.u_[j]*uR0.u_[j]/qbfreq/qbfreq/16;  
		}	
								
		sysRWC.Setucontrol(&uR0,0);	
		sysRWC.Setucontrol(&uR1,1);	
		sysRWC.Setucontrol(&uR2,2);
			//	sysRWC.SetHcontrol(Hc,0);
	//	sysRWC.SetHcontrol(Hc,1);
		sysRWC.UnitaryTransfer();
		 uR1.writefile("RW1.dat"); uR2.writefile("RW2.dat");
		 
		double fidRWC = sysRWC.GetTopFidelity();
		
		uR0.writefile("RWA0int.dat");	
		uR1.writefile("RWA1int.dat");	
		uR2.writefile("RWA2int.dat");	
		
	
		uR0.Interpolate();
		uR1.Interpolate();
		uR2.Interpolate();
		
		sysRWC.SetNumericalParameters(fidelity=0.99999999999, base_a, epsilon, tolerance, max_iter);
		sysRWC.nconfigs_=7;		
		
		sysRWC.UnitaryTransfer( static_cast< Fid >(&GrapeUnitaryQubitLab::Phi4Sub2),  static_cast< GradFid >(&GrapeUnitaryQubitLab::GradPhi4Sub2));
		
	//	sysRWC.sweeptimes("rwagrape.dat");
		
		
		uR0.writefile("RWA0fin.dat");	
		uR1.writefile("RWA1fin.dat");	
		uR2.writefile("RWA2fin.dat");	
		
//		cout << tgate << '\t' << 1- fidRWC<< '\t' << 1-sys.GetTopFidelity() << '\t' << 1 - sysRWC.GetTopFidelity() << '\n';
//		
//		dataout << tgate << '\t' << 1-fidRWC << '\t' << 1-sys.GetTopFidelity() <<  '\t' << abs(1-sys.GetTopFidelity()) << '\t' << 1 - sysRWC.GetTopFidelity()<< '\n';



	}
	
	dataout.close();
	return 0;
	////////////////////////////
	
	
//Lab frame
	size_t nsublab=nsub*10;
	//tgate+=2; //pulse broadening
	
	//use rw controls
	num_time=(size_t)(rwa_num_time)*nsublab;
	dt=tgate/double(num_time);
	vector<double> ucontrollab0(num_time), ucontrollab1(num_time), ucontrollab2(num_time);
//	for(size_t j = 0; j < old_num_time; ++j)
//	{	for(size_t l = 0; l < nlabsubpixels[0]/nRWCsubpixels[0]; ++l)
		{	//ucontrollab0[j*nsublab/nsub+l]=ucontrol0RWC[j];//USs::Gaussian(j*dt, 0, tgate-1, M_PI, tgate/2);	
			//ucontrollab1[j*nsublab/nsub+l]=ucontrol1RWC[j];			
			//ucontrollab2[j*nsublab/nsub+l]=ucontrol2RWC[j];			
		}	
//	}
/*	for(size_t j = 0; j < old_num_time; ++j)
	{	for(size_t l = 0; l < nlabsubpixels[0]; ++l)
		{	ucontrollab0[j*nlabsubpixels[0]+l]=newcontrols0[j];//USs::Gaussian(j*dt, 0, tgate-1, M_PI, tgate/2);	
			ucontrollab1[j*nlabsubpixels[0]+l]=newcontrols1[j];
			ucontrollab2[j*nlabsubpixels[0]+l]=newcontrols2[j];			
		}	
	}*/
	
	//Filter
//	double cntrl, cntrl2, cntrl3;
//	double sig = 0.6;
//	double normalizn;
//	normalizn=0;
	
//	for(size_t l = 0; l < num_time; ++l)					
//	{														
//		normalizn += exp(-abs((int)num_time/2-(int)l)*dt/sig)*dt/sig/2;	
//	}
//	std::cout << normalizn << "\n";
//	vector<double> ucontrollab0filt(num_time), ucontrollab1filt(num_time), ucontrollab2filt(num_time);
//	for(size_t j = 0; j < num_time; ++j)
//	{	
//		cntrl=cntrl2=cntrl3=0;
//		for(size_t k = 0; k < num_time; ++k)					
//		{														
//			cntrl += ucontrollab0[k] * exp(-abs((int)(k)-(int)j)*dt/sig)*dt/sig/2/normalizn;	
//		}														
//		ucontrollab0filt[j]= cntrl;	
										
		//cout << j << " " << ucontrollab0filt[j] << "\n";

//	}
/*	
	GrapeUnitaryQubitLab sysLab(dim, num_controls, num_time);
	//sysLab=sysRWC;
	//cout << fidelity<< dt<< base_a<< epsilon<< tolerance<< max_iter<< nlabsubpixels[2]<< qbfreq<< Delta << endl;
	sysLab.SetPhysicalParameters(fidelity=0.9, dt, base_a, epsilon, tolerance, max_iter, qbfreq, Delta);
	Urot(0,0) = 1;	
	Urot(1,1) = exp(-complex<double>(0.0,tgate*qbfreq));
	Urot(2,2) = exp(-complex<double>(0.0,2.0*tgate*qbfreq));
	Urot(3,3) = exp(-complex<double>(0.0,3.0*tgate*qbfreq));
	
	sysLab.SetRhoDesired(Urot*UNOT);		
	sysLab.SetHdrift(HDrift=qbfreq*n + Delta/2*(n*n-n));
		
	sysLab.Setucontrol(ucontrollab0,0);
	sysLab.Setucontrol(ucontrollab1,1);
	sysLab.Setucontrol(ucontrollab2,2);
//	sysLab.SetHcontrol(a+ad,0);
//	sysLab.SetHcontrol(a+ad,1);
//	sysLab.SetHcontrol(n,2);		
	
	clo = clock();	
	//sysLab.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2Filter, &Grape::SetLab2QuadFiltered, &Grape::SetGradLab2QuadFiltered);
//	sysLab.UnitaryTransfer( static_cast< Fid >(&Grape::Phi4Sub2),  static_cast< GradFid >(&Grape::GradPhi4Sub2), &Grape::SetLab2Quadratures, &Grape::SetGradLab2Quadratures);
	clo = to_ms*(clock() - clo);
	cout << "time " << clo << endl;
		
	*/
	
//	for(size_t j =0; j < rwa_num_time; j++)
//		dataout<< j*rwa_dt << " " << abs(sys.rho_[j](0,0))*abs(sys.rho_[j](0,0)) << " " << abs(sys.rho_[j](1,0))* abs(sys.rho_[j](1,0))  << " " << abs(sys.rho_[j](2,0))*abs(sys.rho_[j](2,0)) << endl;
	
	return 0;
}