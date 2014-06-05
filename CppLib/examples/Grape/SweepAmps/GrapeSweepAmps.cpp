/*
Name: Grape.cpp
Author: Jay Gambetta

Dependences: GrapeUnitary.hpp  
Brief Discription: This demonstrates grape for lab frame vs rotating frame with large pixels

Version History
	v0: Feb  11th, 2009.

time ./GrapeUnitary4.e GrapeUnitary4.dat


n = 10000
real	0m8.287s
user	0m7.917s
sys	0m0.134s

real	0m0.027s
user	0m0.023s
sys	0m0.003s


*/
#include <Grape.hpp>
using namespace std;

int main (int argc, char const *argv[]){

	verbose=no;
	cout << "Running program " << argv[0] << endl;

	//Grape inputs
	double qbfreq=51.5;
	double delta=0.0, Delta=-2.03;
	size_t max_iter=1000;
	double tolerance=std::numeric_limits<double>::min(), fidelity=0.9999, base_a=2.0, epsilon=1, sigma=1, tgate=4, dt;
	size_t num_time=tgate, dim = 4, num_controls =2;
	size_t rwa_num_time=num_time;
	dt=tgate/double(num_time);
	size_t nsubpixels[]={1,1,(size_t)tgate};
	// cout << dt << endl;
	Grape sys(dim, num_controls, num_time);
	sys.SetNumericalParameters(fidelity, dt, base_a, epsilon, tolerance, max_iter, nsubpixels, qbfreq, Delta);
	
	
	matrix<complex<double> > a(dim,dim), ad(dim,dim), n(dim,dim), Hdrift(dim,dim),  HcontrolZ(dim,dim), HcontrolX(dim,dim), HcontrolY(dim,dim),  Htemp(dim,dim);
	matrix<complex<double> > Urot(dim,dim);
	vector<double> ucontrol0(num_time), ucontrol1(num_time);
	matrix<complex<double> > U_desired(dim,dim);
	U_desired.SetOutputStyle(Matrix);
	MOs::Null(U_desired);
	Hdrift.SetOutputStyle(Matrix);
	HcontrolX.SetOutputStyle(Matrix);
	HcontrolY.SetOutputStyle(Matrix);
	MOs::Destroy(a);
	ad = MOs::Dagger(a);
	n=ad*a;

	matrix<complex<double> > Hperf(dim,dim);
	Hperf(0,1)=Hperf(1,0)=Hperf(2,2)=1;

	//RWA frame
	//Hdrift = (delta-0.5*Delta)*n+ 0.5*Delta*n*n;	
	sys.SetHdrift(Hdrift);
	cout << Hdrift << endl;
	
	double lambda1 = 2.8;
	double lambda2 = 6;
		//0th control
	HcontrolX=(a+ad)*0.5;
	HcontrolX(1,0) = HcontrolX(0,1) = -1;
	HcontrolX(1,2) = HcontrolX(2,1) = -lambda1;
	HcontrolX(3,2) = HcontrolX(2,3) = -lambda2;
	
	Hdrift(1,1) =  delta;	
	Hdrift(2,2) =  Delta+delta*2;	
	Hdrift(3,3) = 3.0*Delta+delta*3;	
	
	cout << HcontrolX << endl;
	for(size_t j = 0; j < num_time; ++j)
	{	ucontrol0[j]=USs::Gaussian(j*dt, 0, tgate-1, M_PI, tgate/4)*4;
		ucontrol1[j]=USs::GaussianDerivative(j*dt, 0, tgate-1, M_PI, tgate/4)*4/2;
	}
	sys.Setucontrol(ucontrol0,0);
	//sys.Normalizeucontrol(0, 3.14159265358);
	sys.SetHcontrol(HcontrolX,0);
	//1st control
	HcontrolY=(-complex<double>(0.0,1.0)*a+complex<double>(0.0,1.0)*ad)*0.5;
	HcontrolY(1,2) = complex<double>(0.0,lambda1);
	HcontrolY(2,1) = -complex<double>(0.0,lambda1);
	HcontrolY(3,2) = -complex<double>(0.0,lambda2);
	HcontrolY(2,3) = complex<double>(0.0,lambda2);	

	cout << HcontrolY << endl;	sys.Setucontrol(ucontrol1,1);
	sys.SetHcontrol(HcontrolY,1);
	
	HcontrolZ = n;
	//sys.Setucontrol(ucontrol1,2);
	//sys.SetHcontrol(HcontrolZ,2);
				
	//Initial condition
	//PI pulse
	U_desired(1,0)=std::complex<double>(1,0);
	U_desired(0,1)=std::complex<double>(1,0);
	U_desired(2,2)=std::complex<double>(1,0);
	MOs::Identity(U_desired);
	U_desired = ExpM::EigenMethod(Hperf,std::complex<double>(0,3.1415926*4/2));
	sys.SetRhoDesired(U_desired);	
	
		
	//run grape
	double to_ms =1000.0/CLOCKS_PER_SEC;
	int clo ;
	clo = clock();	
	sys.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2); //, &Grape::SetCounterRotating, &Grape::SetCounterRotating);
	double* newcontrols0 = sys.Getucontrol(0);
	double* newcontrols1 = sys.Getucontrol(1);
	//double* newcontrols2 = sys.Getucontrol(2);

	clo = to_ms*(clock() - clo);
	cout << "time " << clo << endl;
	
	vector<double> ucontrol0final, ucontrol1final;
	// The finial control
	ucontrol0final.assign(sys.Getucontrol(0), sys.Getucontrol(0) + num_time );
	ucontrol1final.assign(sys.Getucontrol(1), sys.Getucontrol(1) + num_time );

	newcontrols0[0]=2.1549284941378857e+00;
	newcontrols1[0]=2.1536126521533370e+00;
	newcontrols0[1]=5.6067694780564850e+00;
	newcontrols1[1]=3.2114511481479249e+00;
	newcontrols0[2]=5.4047370482193218e+00;
	newcontrols1[2]=-3.6284422906264950e+00;
	newcontrols0[3]=2.0259643975393358e+00;
	newcontrols1[3]=-2.3947284073168538e+00;
	size_t rwsub=200;
	size_t nRWCsubpixels[]={rwsub,rwsub,(size_t)tgate*rwsub};
	size_t old_num_time = num_time;
	num_time=(size_t)rwa_num_time*rwsub;
	dt=tgate/double(num_time);
	cout << dt << ' ' << tgate << " " << num_time << "\n";
	vector<double> ucontrolRWC0(num_time), ucontrolRWC1(num_time), ucontrolRWC2(num_time), ucontrolRWC0filt(num_time), ucontrolRWC1filt(num_time) ;
	for(size_t j = 0; j < old_num_time/nsubpixels[0]; ++j)
	{	for(size_t l = 0; l < nRWCsubpixels[0]; ++l)
		{	ucontrolRWC0[j*nRWCsubpixels[0]+l]=newcontrols0[j*nsubpixels[0]];//USs::Gaussian(j*dt, 0, tgate-1, M_PI, tgate/2);	
			ucontrolRWC1[j*nRWCsubpixels[0]+l]=newcontrols1[j*nsubpixels[0]];
		//	ucontrolRWC2[j*nRWCsubpixels[0]+l]=newcontrols2[j*nsubpixels[0]];			
		}	
	}
	double cntrl, cntrl2;
	double sig = 0.5;
	double normalizn;
	normalizn=0;
	for(size_t l = 0; l < num_time; ++l)					
	{														
		normalizn += exp(-abs((int)num_time/2-(int)l)*dt/sig)*dt/sig/2;	
	}

	for(size_t j = 0; j < num_time; ++j)
	{	
		cntrl=cntrl2=0;
		for(size_t k = 0; k < num_time; ++k)					
		{														
			cntrl += ucontrolRWC0[k] * exp(-abs((int)(k)-(int)j)*dt/sig)*dt/sig/2/normalizn;	
			cntrl2 +=  ucontrolRWC1[k] * exp(-abs((int)(k)-(int)j)*dt/sig)*dt/sig/2/normalizn;			
		}														
		ucontrolRWC0filt[j]= cntrl;						
		ucontrolRWC1filt[j]= cntrl2;								
		//cout << j << " " << ucontrollab0filt[j] << "\n";

	}

	string const outfile = argv[1];
	ofstream dataout;
	UFs::OpenFile(outfile,dataout, 16);
	
		
	matrix<complex<double> > unitary(dim,dim);
	/*for(double amplif=0.05; amplif<2; amplif+=0.05)
	{
	
		MOs::Identity(unitary);
		for(size_t j = 0; j < num_time; ++j){
			Htemp = Hdrift;			         
			for(size_t q=0; q< dim; ++q) {
				for(size_t p=0; p<dim; ++p)
						Htemp(q,p) += amplif*ucontrol0[j]*HcontrolX(q,p);
			}
			unitary = ExpM::EigenMethod(Htemp,-std::complex<double>(0,dt))*unitary;
		}	
		U_desired = ExpM::EigenMethod(Hperf,std::complex<double>(0,amplif*3.1415926*3/2));
		unitary = unitary * U_desired;
		std::cout <<  amplif << " " << real(( unitary(0,0)+unitary(1,1) ) * conj( unitary(0,0)+unitary(1,1) )/4.0) << "\n";		
		dataout << real(( unitary(0,0)+unitary(1,1) ) * conj( unitary(0,0)+unitary(1,1) )/4.0)<< endl;
		//dataout << amplif << " " << abs(unitary(0,0)) << " " << abs(unitary(1,0)) << " " <<  abs(unitary(2,0)) << " " << endl;
		//std::cout << amplif << " " << abs(unitary(0,0)) << " " << abs(unitary(1,0)) << " " << abs (unitary(2,0)) << " " << endl << endl;
	}
	
	cout <<  "\n\n";
	*/
	//matrix<complex<double> > unitary(dim,dim);
	for(double amplif=0.001; amplif<1.3; amplif+=0.02)
	{
		MOs::Identity(unitary);
		for(size_t j = 0; j < num_time; ++j){
			Htemp = Hdrift;			         
			for(size_t q=0; q< dim; ++q) {
				for(size_t p=0; p<dim; ++p)
						Htemp(q,p) += amplif*ucontrolRWC0filt[j]*HcontrolX(q,p)+ amplif*ucontrolRWC1filt[j]*HcontrolY(q,p);
			}
			
			unitary = ExpM::EigenMethod(Htemp,-std::complex<double>(0,dt))*unitary;
			//std::cout <<unitary << "\n";
		}	
		U_desired = ExpM::EigenMethod(Hperf, std::complex<double>(0,amplif*3.1415926*4/2));
		//std::cout << unitary << "\n";
		//unitary = unitary * U_desired;
		//std::cout << amplif << " " << real(( unitary(0,0)+unitary(1,1) ) * conj( unitary(0,0)+unitary(1,1) )/4.0)<< "\n";
		//dataout << amplif << " " << real(( unitary(0,0)+unitary(1,1) ) * conj( unitary(0,0)+unitary(1,1) )/4.0)<< endl;
		dataout << amplif << " " << abs(unitary(0,0)) << " " << abs(unitary(1,0)) << " " <<  abs(unitary(2,0)) << " " << endl;
		std::cout << amplif << " " << abs(unitary(0,0)) << " " << abs(unitary(1,0)) << " " <<  abs(unitary(2,0)) << " " << endl << endl;
	}
	
	
	dataout.close();
	return 0;
}