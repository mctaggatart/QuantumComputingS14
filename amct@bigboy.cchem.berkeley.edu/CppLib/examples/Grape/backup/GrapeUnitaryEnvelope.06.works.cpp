/*
Name: GrapeUnitaryEnvelope.cpp
Author: felix motzoi

Dependences: GrapeUnitaryEnvelope.hpp  
Brief Discription: This demonstrates grape for lab frame vs rotating frame with large pixels

Version History
	v0: Using Control class, Oct. 6 2009.

*/
#include <Grape.hpp>
using namespace std;

int main (int argc, char const *argv[]){

	verbose=no;
	cout << "Running program " << argv[0] << endl;

	//System/Optimization parameters
	double qbfreq=51.5, Delta=-200.03, lambda=sqrt(2);  //physical params
	double tolerance=0.00000000000000001, fidelity=0.9, base_a=2.0, epsilon=1, tgate=1, dt;  //search params
	size_t max_iter=1000;
	size_t dim=4;
	double to_ms =1000.0/CLOCKS_PER_SEC;
	int clo;
	
	//Typical operators	
	complex<double> ii(0.0,1.0);
	matrix<complex<double> > a(dim,dim), ad(dim,dim), n(dim,dim), UNOT(dim,dim);
	matrix<complex<double> > Urot(dim,dim);
	MOs::GenPauliX(UNOT,0,1);
	MOs::Destroy(a);
	ad = MOs::Dagger(a);
	n=ad*a;

	//Grape inputs
	string const outfile = argv[1];
	ofstream dataout;
	UFs::OpenFile(outfile,dataout, 16);
	
	
	//RWA fixed parameters
	const size_t rwasub=1;
	size_t rwa_num_time=tgate*rwasub*1, num_controls =3;
	
	//RW fixed parameters
    const int nsub=2000;
	size_t num_time=rwa_num_time*nsub/rwasub;
	dt=tgate/double(num_time);
	
	
	for(size_t itime=0; itime<14; itime++, tgate++, num_time+=nsub, rwa_num_time+=rwasub)
	{
	
		//RWA frame	
		double rwa_dt = tgate/double(rwa_num_time);
		Grape sys(dim, num_controls, rwa_num_time);
		sys.SetNumericalParameters(fidelity=0.09999999, rwa_dt, base_a, epsilon, tolerance, max_iter, qbfreq, Delta);
		
		sys.SetHdrift(Delta/2*(n*n-n));
		
		Control u0(0.5*(a+ad), rwa_dt,rwa_num_time,rwasub), u1(0.5*(-ii*a+ii*ad),rwa_dt, rwa_num_time,rwasub), u2(n,rwa_dt,rwa_num_time,rwasub);
		u0.ShiftedGaussian( 0, tgate, M_PI, tgate/2);
		double norm = u0.Normalize(M_PI);
		u1.GaussianDerivative(0, tgate, -M_PI/Delta*norm, tgate/2);
		for(size_t j = 0; j < rwa_num_time; ++j)
			u2.u_[j]=(lambda*lambda-4)*(u0.u_[j]*u0.u_[j])/4/Delta;
		
		u0.settotaltime(rwa_num_time);
		u0.writefile("RWA0.dat");	
		u1.settotaltime(rwa_num_time);
		u2.settotaltime(rwa_num_time);
		sys.Setucontrol(u0, 0); //0th control
		sys.Setucontrol(u1,1); //1st control
		sys.Setucontrol(u2,2); //2nd control
			
		sys.SetRhoDesired(UNOT);	
		
		//run grape
		clo= clock();	
		sys.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2); 
		
		dt=tgate/double(num_time);		
		Control uR0(u0,nsub), uR1(u1,nsub), uR2(u2,nsub);
	    uR0.writefile("RW0.dat");
		Grape sysRWC(dim, num_controls, num_time);	
		sysRWC.SetNumericalParameters(fidelity=0.1, dt, base_a, epsilon, tolerance, max_iter, qbfreq, Delta);
		sysRWC.SetRhoDesired(UNOT);	
		sysRWC.nconfigs_=9;		
		sysRWC.SetHdrift(Delta/2*(n*n-n));
						//
//		uR0.settotaltime(num_time);
//		uR1.settotaltime(num_time);
//		uR2.settotaltime(num_time);
//	//	
//		uR0.ShiftedGaussian( 0, tgate, M_PI, tgate/2);
//		norm = uR0.Normalize(M_PI);
//		uR1.GaussianDerivative(0, tgate, -M_PI/Delta*norm, tgate/2);
//		for(size_t j = 0; j < num_time; ++j)
//			uR2.u_[j]=(lambda*lambda-4)*(uR0.u_[j]*uR0.u_[j])/4/Delta;
//		
		sysRWC.Setucontrol(uR0,0);	
		sysRWC.Setucontrol(uR1,1);	
		sysRWC.Setucontrol(uR2,2);
		sysRWC.SetHcontrol(a+ad,0);
		sysRWC.SetHcontrol(a+ad,1);

		sysRWC.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2,&Grape::SetCounterRotating, &Grape::SetGradCounterRotating);
		 uR1.writefile("RW1.dat"); uR2.writefile("RW2.dat");
		cout << tgate << '\t' << 1-sysRWC.top_fidelity_ << '\t' << 1-sys.top_fidelity_ <<  '\t' << abs(sys.top_fidelity_-sysRWC.top_fidelity_) << '\n';
		dataout << tgate << '\t' << 1-sysRWC.top_fidelity_ << '\t' << 1-sys.top_fidelity_ <<  '\t' << abs(sys.top_fidelity_-sysRWC.top_fidelity_) << '\n';
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
	
	Grape sysLab(dim, num_controls, num_time);
	//sysLab=sysRWC;
	//cout << fidelity<< dt<< base_a<< epsilon<< tolerance<< max_iter<< nlabsubpixels[2]<< qbfreq<< Delta << endl;
	sysLab.SetNumericalParameters(fidelity=0.9, dt, base_a, epsilon, tolerance, max_iter, qbfreq, Delta);
	Urot(0,0) = 1;	
	Urot(1,1) = exp(-complex<double>(0.0,tgate*qbfreq));
	Urot(2,2) = exp(-complex<double>(0.0,2.0*tgate*qbfreq));
	Urot(3,3) = exp(-complex<double>(0.0,3.0*tgate*qbfreq));
	
	sysLab.SetRhoDesired(Urot*UNOT);		
	sysLab.SetHdrift(qbfreq*n + Delta/2*(n*n-n));
		
	sysLab.Setucontrol(ucontrollab0,0);
	sysLab.Setucontrol(ucontrollab1,1);
	sysLab.Setucontrol(ucontrollab2,2);
	sysLab.SetHcontrol(a+ad,0);
	sysLab.SetHcontrol(a+ad,1);
	sysLab.SetHcontrol(n,2);		
	
	clo = clock();	
	//sysLab.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2Filter, &Grape::SetLab2QuadFiltered, &Grape::SetGradLab2QuadFiltered);
	sysLab.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2, &Grape::SetLab2Quadratures, &Grape::SetGradLab2Quadratures);
	clo = to_ms*(clock() - clo);
	cout << "time " << clo << endl;
		
	
	
//	for(size_t j =0; j < rwa_num_time; j++)
//		dataout<< j*rwa_dt << " " << abs(sys.rho_[j](0,0))*abs(sys.rho_[j](0,0)) << " " << abs(sys.rho_[j](1,0))* abs(sys.rho_[j](1,0))  << " " << abs(sys.rho_[j](2,0))*abs(sys.rho_[j](2,0)) << endl;
	
	return 0;
}