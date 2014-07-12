/*
Name: GrapeSweepAmps2.cpp
Author: Felix Motzoi

Dependences: GrapeSweepAmps2.hpp  
Brief Discription: This is used for testing sweeping various amplitude parameters / benchmarking
Version History
	v0: Mar  31th, 2009.


*/
#include <Grape.hpp>
using namespace std;

int main (int argc, char const *argv[]){

	verbose=no;
	cout << "Running program " << argv[0] << endl;
	
	//constants
	double to_ms =1000.0/CLOCKS_PER_SEC;
	int clo ;
	
	
	double fidel;
	complex<double> Udiag[2][4];
	
	double prefacX1[] = {0,0,1,1,0,0,-1,-1,0,0,	0.5,0.5,0,0,-0.5,-0.5,0,0,	1,1,1,1,1,1,1,1,	0,0,0,0,0,0,0,0,	1,1,1,1,1,1,1,1,			0,0,0,0,0,0,0,0,			0.5,0.5,0.5,0.5,0,0,0,0,	0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,	0,0,0,0,0,0,0,0};
	double prefacY1[] = {0,0,0,0,1,1,0,0,-1,-1,	0,0,0.5,0.5,0,0,-0.5,-0.5,	0,0,0,0,0,0,0,0,	1,1,1,1,1,1,1,1,	0,0,0,0,0,0,0,0,			1,1,1,1,1,1,1,1,			0,0,0,0,0.5,0.5,0.5,0.5,	0,0,0,0,0,0,0,0,				0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5};
	double prefacX2[] = {0,0,0,0,0,0,0,0,0,0,	0,0,0,0,0,0,0,0,			1,1,-1,-1,0,0,0,0,	1,1,-1,-1,0,0,0,0,	0.5,0.5,-0.5,-0.5,0,0,0,0,	0.5,0.5,-0.5,-0.5,0,0,0,0,	1, 1, 0, 0, 1, 1, 0, 0,		0.5,0.5,-0.5,-0.5,0,0,0,0,			0.5,0.5,-0.5,-0.5,0,0,0,0};
	double prefacY2[] = {0,0,0,0,0,0,0,0,0,0,	0,0,0,0,0,0,0,0,			0,0,0,0,1,1,-1,-1,	0,0,0,0,1,1,-1,-1,	0,0,0,0,0.5,0.5,-0.5,-0.5,	0,0,0,0,0.5,0.5,-0.5,-0.5,	0, 0, 1, 1, 0, 0, 1, 1,		0,0,0,0,0.5,0.5,-0.5,-0.5,			0,0,0,0,0.5,0.5,-0.5,-0.5};	
	
	double Xamps[]={0,0,1,-1,0,0,0,0,0.5,-0.5,0,0};	
	double Yamps[]={0,0,0,0,1,-1,0.5,-0.5,0,0,0,0};
	double fidelitybins[17];
	double timebins[17];
	double countbins[17];
	double clifffidelitybins[4][17];
	double clifftimebins[4][17];
	double cliffcountbins[4][17];	
	double xavg=0, x2avg=0, lnyavg=0, xlnyavg=0, xavglny=0;
	double xyavg=0, x2yavg=0, ylnyavg=0, xylnyavg=0, xyavgylny=0, yavgxylny, yavgx2y;
	
	ifstream benching;
	benching.open("benching.dat");
	double Xseqs[544*100];
	double Yseqs[544*100];
	size_t finalstates[544];
	unsigned seqnum=0;
	char separator, space, quad, sign;
	unsigned pulsenum;
	double ampX, ampY;
	unsigned nsteps=0;
	while(benching.good())
	{
		finalstates[seqnum] = benching.get() - '0';
		//cout <<	finalstate	<< "\n";			
		if(!benching.good()) break;
		separator=' ';
		separator=benching.get();
		while (separator!='}')
		{
			ampX=ampY=0;
			if(!benching.good()) break;
			quad = benching.get();
			//cout << 'n';
			if(quad=='\n' || quad=='\r')
			{	//cout << 'l';
				if(!benching.good()) break;
				space =  benching.get();
				if(!benching.good()) break;
				quad = benching.get();
			}
			if(quad=='X') ampX=1;
			else if(quad=='Y') ampY=1;
			else space =  benching.get();
			if(!benching.good()) break;
			sign =  benching.get();
			//cout << quad;
			if(sign=='9')
			{	ampX/=2;
				ampY/=2;
			//cout << sign;
			if(!benching.good()) break;
				space =  benching.get();
			if(!benching.good()) break;
				sign =  benching.get();		  
			}
			//cout << sign;
			if(sign=='m')
			{	ampX*=-1;
				ampY*=-1;		  
			}
			//cout << ' ';
			if(!benching.good()) break;
			separator =  benching.get();
			if(!benching.good()) break;
			space =  benching.get();
			//cout << 'm';
			Xseqs[nsteps]=ampX;
			Yseqs[nsteps++]=ampY;
		}
		Xseqs[nsteps]=55;
		Yseqs[nsteps++]=55;
		seqnum++;
	}	
	benching.close();
	Xseqs[nsteps]=110;
	Yseqs[nsteps++]=110;
	
	
	string const outfile = argv[1];
	ofstream dataout;
	UFs::OpenFile(outfile,dataout, 16);
	
	//Grape inputs
	double detune=0;//0.07/6;  //driving frequency detuning
	double qbfreq=51.5-detune;  //frame frequency
	//double lambda1=1.88, lambda2=3.47;
	double lambda1=2.45, lambda2=5.5;
	double delta=detune, Delta=-2.03;  //relative energies
	size_t max_iter=2000, dim = 4, num_controls =3;
	double tolerance=0.000000000000000000000000001, fidelity=0.0009999, base_a=2.0, epsilon=1, sigma=1, tgate, dt;
	double nPI=1;  //area of pulse
	double *newcontrols0, *newcontrols1, *newcontrols2, *endcontrols0, *endcontrols1, *endcontrols2;
	
	matrix<complex<double> > a(dim,dim), ad(dim,dim), n(dim,dim), Hdrift(dim,dim), Hdrift_lab(dim,dim), HcontrolZ(dim,dim), HcontrolX(dim,dim), HcontrolY(dim,dim), Htemp(dim,dim);
	matrix<complex<double> > Urot(dim,dim);
	matrix<complex<double> > U_desired(dim,dim), U_desired_lab(dim,dim);
	matrix<complex<double> > unitary(dim,dim);
	unitary.SetOutputStyle(Matrix);
	U_desired.SetOutputStyle(Matrix);
	MOs::Null(U_desired);
	Hdrift.SetOutputStyle(Matrix);
	Hdrift_lab.SetOutputStyle(Matrix);
	HcontrolX.SetOutputStyle(Matrix);
	HcontrolY.SetOutputStyle(Matrix);
	HcontrolZ.SetOutputStyle(Matrix);
	MOs::Destroy(a);
	ad = MOs::Dagger(a);
	n=ad*a;
	
	//Pick desired Gate
	matrix<complex<double> > Hperf(dim,dim);
	Hperf(0,1)=Hperf(1,0)=Hperf(2,2)=1;
	U_desired = ExpM::EigenMethod(Hperf, std::complex<double>(0,M_PI*round(nPI)/2));	
	
	//0th control
	HcontrolX=(a+ad);
	HcontrolY=(-complex<double>(0.0,1.0)*a+complex<double>(0.0,1.0)*ad)*0.5;		
	HcontrolZ = n;
	//cavity below
	HcontrolX(1,0) = HcontrolX(0,1) = -1;
	HcontrolX(1,2) = HcontrolX(2,1) = -lambda1;
	HcontrolX(3,2) = HcontrolX(2,3) = -lambda2;
	HcontrolY(1,0) = -complex<double>(0.0,1.0);
	HcontrolY(0,1) = complex<double>(0.0,1.0);
	HcontrolY(1,2) = complex<double>(0.0,lambda1);
	HcontrolY(2,1) = -complex<double>(0.0,lambda1);
	HcontrolY(3,2) = -complex<double>(0.0,lambda2);
	HcontrolY(2,3) = complex<double>(0.0,lambda2);	
	HcontrolX = HcontrolX*0.5;
	HcontrolY = HcontrolY*0.5;
	HcontrolZ = n;
	cout << HcontrolX << endl;	
	cout << HcontrolY << endl;
	cout << HcontrolZ << endl;
	
	vector<double> ucontrol0(21*12), ucontrol1(21*12), ucontrol2(21*12);
	size_t num_time, rwa_num_time;
	size_t nsubpixels[num_controls];
	double rwa_dt;
	Grape* sys;
	for(double old_tgate=4; old_tgate<10; old_tgate++)
	{
		tgate=old_tgate;
		num_time=tgate*12;
		rwa_num_time = num_time;
		dt=tgate/double(num_time);
		rwa_dt = dt;
		nsubpixels[0]=nsubpixels[1]=1;
		nsubpixels[2]=num_time;
		// cout << dt << endl;
		sys = new Grape(dim, num_controls, num_time);
		cout << dt << ' ' << tgate << " " << num_time << "\n";
		sys->SetNumericalParameters(fidelity, dt, base_a, epsilon, tolerance, max_iter, nsubpixels, qbfreq, Delta);

		//RWA frame
		//Hdrift = (delta-0.5*Delta)*n+ 0.5*Delta*n*n;	
		Hdrift(1,1) =  delta;	
		Hdrift(2,2) =  Delta+delta*2;	
		Hdrift(3,3) = 3.0*Delta+delta*3;	

		sys->SetHdrift(Hdrift);
		cout << Hdrift << endl;	

		for(size_t j = 0; j < num_time; ++j){
			ucontrol0[j]=USs::ShiftedGaussian((j+0.5)*dt, 0, tgate, M_PI, tgate/4)*nPI;
	//		ucontrol0[j]=USs::Gaussian((j+0.5)*dt, 0, tgate, M_PI, tgate/4)*nPI;
			ucontrol1[j]=-USs::GaussianDerivative((j+0.5)*dt, 0, tgate, M_PI, tgate/4)*nPI/Delta*1.15;//*2.95;//*USs::Gaussian(j*dt, 0, tgate-1, M_PI, tgate/4)*nPI ; //*0.58;
		//	ucontrol2[j]=-0.2;
		//	cout << ucontrol0[j] << " " << ucontrol1[j] << endl;
		}
	
		//0th control
		sys->Setucontrol(ucontrol0,0);
		sys->Normalizeucontrol(0, 3.14159265358*nPI);
		sys->SetHcontrol(HcontrolX,0);
		//1st control
		sys->Setucontrol(ucontrol1,1);
		sys->SetHcontrol(HcontrolY,1);	
		//2nd control
		sys->Setucontrol(ucontrol2,2);
		sys->SetHcontrol(0.0*HcontrolZ,2);						
		sys->SetRhoDesired(U_desired);	
	
		//run grape
		clo = clock();	
		sys->UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2); 
		clo = to_ms*(clock() - clo);
		cout << "RWA time " << clo << endl;
		
		newcontrols0 = sys->Getucontrol(0);
		newcontrols1 = sys->Getucontrol(1);
//		newcontrols2 = sys->Getucontrol(2);
	
		//for(size_t j = 0; j < num_time; ++j)
		//	cout << newcontrols0[j] << " " << newcontrols1[j] << endl;
		MOs::Identity(unitary);	
		for(size_t j = 0; j < rwa_num_time; ++j){
			
//			Udiag[0][0]=Udiag[1][0]=1;
//			Udiag[0][1]=exp(-std::complex<double>(0.0,(double)j*qbfreq*dt));
//			Udiag[0][2]=exp(-std::complex<double>(0.0,2.0*(double)j*qbfreq*dt));
//			Udiag[0][3]=exp(-std::complex<double>(0.0,3.0*(double)j*qbfreq*dt));
//			Udiag[1][1]=exp(std::complex<double>(0.0,(double)j*qbfreq*dt));
//			Udiag[1][2]=exp(std::complex<double>(0.0,2.0*(double)j*qbfreq*dt));	
//			Udiag[1][3]=exp(std::complex<double>(0.0,3.0*(double)j*qbfreq*dt));	
			
			Htemp = Hdrift;			  
			 for(size_t q=0; q< dim; ++q) {
				for(size_t p=0; p<dim; ++p){
						//Htemp(p,q) += Xamps[pulsenum]*Udiag[0][q]*Udiag[1][p]*(ctrl0[j]*HcontrolX(p,q)*cos(((double)j)*qbfreq*dt) + ctrl1[j]*HcontrolX(p,q)*sin(((double)j)*qbfreq*dt) );
						//Htemp(p,q) += Yamps[pulsenum]*Udiag[0][q]*Udiag[1][p]*(ctrl1[j]*HcontrolX(p,q)*cos(((double)j)*qbfreq*dt) - ctrl0[j]*HcontrolX(p,q)*sin(((double)j)*qbfreq*dt) );
						Htemp(p,q) += (newcontrols0[j]*HcontrolX(p,q) + newcontrols1[j]*HcontrolY(p,q));
						//Htemp(q,p) += amplif*(ucontrollab0filt[j]*HcontrolX(q,p)+ ucontrollab1filt[j]*HcontrolY(q,p));
				}
			 }
			 unitary = ExpM::EigenMethod(Htemp,-std::complex<double>(0,dt))*unitary; ///ddt		
			 for(size_t k = 0; k < 10; ++k){
		///		pop = unitary*pop*MOs::Dagger(unitary); + (ddt/1000)*(a*pop*ad -0.5*ad*a*pop -0.5*pop*ad*a);
				//pop +=  -std::complex<double>(0,ddt)*Htemp*pop + std::complex<double>(0,ddt)*pop*Htemp +  (ddt/1000)*(a*pop*ad -0.5*ad*a*pop -0.5*pop*ad*a);
			 }
		}										
		fidel = abs( unitary(0,1) * unitary(1,0) );
		cout << "avg error " << 1.0-  fidel << endl;	
		//dataout << 	old_tgate << '\t' <<1.0-  fidel << "\t" << 1-sys->fidelity_ <<endl;
		//delete(sys);
		//continue;
/*	//YALE amplitude sweep test
	newcontrols0[0]=2.1549284941378857e+00;
	newcontrols1[0]=2.1536126521533370e+00;
	newcontrols0[1]=5.6067694780564850e+00;
	newcontrols1[1]=3.2114511481479249e+00;
	newcontrols0[2]=5.4047370482193218e+00;
	newcontrols1[2]=-3.6284422906264950e+00;
	newcontrols0[3]=2.0259643975393358e+00;
	newcontrols1[3]=-2.3947284073168538e+00;
*/
	
//RW with Counter-Rotating
/*	size_t rwsub=1;
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
			ucontrolRWC2[j*nRWCsubpixels[0]+l]=newcontrols2[j*nsubpixels[0]];			
		}	
	}
	*/
		//Filter
/*	double cntrl, cntrl2;
	double sig = 0.03;
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

	}	*/
	
//	nRWCsubpixels[0]/=2;
//	nRWCsubpixels[1]/=2;
/*	
	for(size_t j = 0; j < num_time; ++j){
		ucontrolRWC0filt[j]=ucontrolRWC0filt[(size_t)((j/nRWCsubpixels[0]+0.5)*nRWCsubpixels[0])];
		ucontrolRWC1filt[j]=ucontrolRWC1filt[(size_t)((j/nRWCsubpixels[0]+0.5)*nRWCsubpixels[0])];
	}
	fidelity=0.9;
	//fidelity=0.050;
	Grape sysRWC(dim, num_controls, num_time);
	sysRWC.SetNumericalParameters(fidelity, dt, base_a, epsilon, tolerance, max_iter, nRWCsubpixels, qbfreq, Delta);
	sysRWC.SetRhoDesired(U_desired);	
	
	sysRWC.SetHdrift(Hdrift);
	sysRWC.SetHcontrol(2.0*HcontrolX,0);
	sysRWC.SetHcontrol(2.0*HcontrolX,1);
	sysRWC.SetHcontrol(0.0*HcontrolZ,2);		
	sysRWC.Setucontrol(ucontrolRWC0,0);
	sysRWC.Setucontrol(ucontrolRWC1,1);
	sysRWC.Setucontrol(ucontrolRWC2,2);
	
	//sysRWC.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2RWFilter, &Grape::SetCounterRotatingFilt, &Grape::SetGradCounterRotatingFilt);
	//sysRWC.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2, &Grape::SetCounterRotating, &Grape::SetGradCounterRotating);
	
	vector<double> ucontrol0RWC, ucontrol1RWC, ucontrol2RWC;
	ucontrol0RWC.assign(sysRWC.Getucontrol(0), sysRWC.Getucontrol(0) + num_time );
	ucontrol1RWC.assign(sysRWC.Getucontrol(1), sysRWC.Getucontrol(1) + num_time );
	ucontrol2RWC.assign(sysRWC.Getucontrol(2), sysRWC.Getucontrol(2) + num_time );
	
	cout << ucontrol0RWC[0] << " " << ucontrol0RWC[201] << " " << ucontrol0RWC[401]  << "\n"; 
	
	cout << "detuning " << ucontrol2RWC[1] << endl;
*/	
//Lab frame
/*	old_num_time = num_time;
	tgate=old_tgate+2.0; //pulse broadening
	size_t nlabsubpixels[]={rwsub,rwsub,(size_t)tgate*rwsub};
	num_time=(size_t)(rwsub*rwa_num_time/old_tgate*(double)tgate);
	dt=tgate/double(num_time);
	vector<double> ucontrollab0(num_time), ucontrollab1(num_time), ucontrollab2(num_time);
	
	for(size_t j = 0; j < old_num_time; ++j)
	{	for(size_t l = 0; l < nlabsubpixels[0]/nRWCsubpixels[0]; ++l)
		{	ucontrollab0[j*nlabsubpixels[0]/nRWCsubpixels[0]+l+2*nlabsubpixels[0]*(rwa_num_time/(old_tgate))]=ucontrolRWC0[j];//USs::Gaussian(j*dt, 0, tgate-1, M_PI, tgate/2);	
			ucontrollab1[j*nlabsubpixels[0]/nRWCsubpixels[0]+l+2*nlabsubpixels[0]*(rwa_num_time/(old_tgate))]=ucontrolRWC1[j];			
			ucontrollab2[j*nlabsubpixels[0]/nRWCsubpixels[0]+l+2*nlabsubpixels[0]*(rwa_num_time/(old_tgate))]=ucontrolRWC2[j];			
		}	
	}
	cout << dt << ' '<< tgate << " " << num_time << " " << (qbfreq/2/M_PI/8/dt) << "\n";
	
	//Filter
	vector<double> ucontrollab0filt(num_time), ucontrollab1filt(num_time), ucontrollab2filt(num_time);
/*	std::cout << normalizn << "\n";
	for(size_t j = 0; j < num_time; ++j)
	{	
		cntrl=cntrl2=0;
		double endtime = j+300;
		if(endtime>num_time) endtime = num_time;
		for(size_t k = (j>300)*(j-300); k < endtime; ++k)					
		{														
			cntrl += ucontrollab0[k] * exp(-abs((int)(k)-(int)j)*dt/sig)*dt/sig/2/normalizn;	
			cntrl2 +=  ucontrollab1[k] * exp(-abs((int)(k)-(int)j)*dt/sig)*dt/sig/2/normalizn;			
		}														
		ucontrollab0filt[j]= cntrl;						
		ucontrollab1filt[j]= cntrl2;								
		//cout << j << " " << ucontrollab0filt[j] << "\n";
	} */
/*	for(size_t j = 0; j < num_time; ++j){
		ucontrollab0filt[j]=ucontrollab0filt[(size_t)((j/nlabsubpixels[0]+0.5)*nlabsubpixels[0])];
		ucontrollab1filt[j]=ucontrollab1filt[(size_t)((j/nlabsubpixels[0]+0.5)*nlabsubpixels[0])];
	}
	Grape sysLab(dim, num_controls, num_time);
	sysLab.SetNumericalParameters(0.9999, dt, base_a, epsilon, tolerance, max_iter, nlabsubpixels, qbfreq, Delta);
	Urot(0,0) = 1;
	Urot(1,1) = exp(-complex<double>(0.0,tgate*(qbfreq- ucontrollab2[nlabsubpixels[0]])));
	Urot(2,2) = exp(-complex<double>(0.0,2.0*tgate*(qbfreq- ucontrollab2[nlabsubpixels[0]])));
	Urot(3,3) = exp(-complex<double>(0.0,3.0*tgate*(qbfreq- ucontrollab2[nlabsubpixels[0]])));
	
	U_desired_lab=Urot*U_desired;
	
	sysLab.SetRhoDesired(U_desired);
	
	Hdrift_lab(1,1) = qbfreq + delta;
	Hdrift_lab(2,2) = 2.0*qbfreq + Delta + delta*2;	
	Hdrift_lab(3,3) = 3.0*qbfreq + 3*Delta + delta*3;

	sysLab.SetHdrift(Hdrift);
	cout << Hdrift << endl;

	sysLab.Setucontrol(ucontrollab0,0);
	sysLab.Setucontrol(ucontrollab1,1);
	sysLab.Setucontrol(ucontrollab2,2);

	sysLab.SetHcontrol(2.0*HcontrolX,0);
	sysLab.SetHcontrol(2.0*HcontrolX,1);	
	sysLab.SetHcontrol(0.0*HcontrolZ,2);	
	
	//sysLab.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2, &Grape::SetCounterRotating, &Grape::SetGradCounterRotating);
	//sysLab.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2RWFilter, &Grape::SetCounterRotatingFilt, &Grape::SetGradCounterRotatingFilt);
	//sysLab.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2RWFilter, &Grape::SetCounterRotatingFilt, &Grape::SetGradCounterRotatingFilt);
	//sysLab.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2Filter, &Grape::SetLab2QuadFiltered, &Grape::SetGradLab2QuadFiltered);
	//sysLab.UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2, &Grape::SetLab2Quadratures, &Grape::SetGradLab2Quadratures);
	
	endcontrols0 = sysLab.Getucontrol(0);
	endcontrols1 = sysLab.Getucontrol(1);
	endcontrols2 = sysLab.Getucontrol(2);*/
		
	//cavity below
	HcontrolX(1,0) = HcontrolX(0,1) = -1;
	HcontrolX(1,2) = HcontrolX(2,1) = -lambda1;
	HcontrolX(3,2) = HcontrolX(2,3) = -lambda2;
	HcontrolY(1,0) = -complex<double>(0.0,1.0);
	HcontrolY(0,1) = complex<double>(0.0,1.0);
	HcontrolY(1,2) = complex<double>(0.0,lambda1);
	HcontrolY(2,1) = -complex<double>(0.0,lambda1);
	HcontrolY(3,2) = -complex<double>(0.0,lambda2);
	HcontrolY(2,3) = complex<double>(0.0,lambda2);

//RWA
	HcontrolX = HcontrolX*0.5;
	HcontrolY = HcontrolY*0.5;
	
	cout << HcontrolY << endl;	
	//cavity above
	//HcontrolX(1,0) = HcontrolX(0,1) = 1;
	//HcontrolX(1,2) = HcontrolX(2,1) = 1.067;
	//HcontrolX(3,2) = HcontrolX(2,3) = 1.05;
	//HcontrolY(1,0) = HcontrolYX(0,1) = 1;
	//HcontrolY(1,2) = HcontrolY(2,1) = 1.067;
	//HcontrolY(3,2) = HcontrolY(2,3) = 1.05;	
	
	matrix<complex<double> > blank(dim,dim), pop(dim,dim);
	a.SetOutputStyle(Matrix);
	pop.SetOutputStyle(Matrix);
	
	for(int j=0; j<17; j++)
	{	fidelitybins[j]=0;
		timebins[j]=0;
		countbins[j]=0;
	}
	for(int k=0; k<4; k++)
	for(int j=0; j<17; j++)
	{	clifffidelitybins[k][j]=0;
		clifftimebins[k][j]=0;
		cliffcountbins[k][j]=0;
	}
	
	xavg=x2avg=lnyavg=xlnyavg=xavglny=0;
	xyavg=x2yavg=ylnyavg=xylnyavg=xyavgylny=yavgxylny=yavgx2y=0;
	
	double* ctrl0 = newcontrols0;
	//vector<double>
	double* ctrl1 = newcontrols1;
	unsigned index=0;
	pop.SetOutputStyle(Matrix);
	unsigned gatenum=0;
	seqnum=0;
	unsigned finalstate;
	while(Xseqs[gatenum]!=110)
	{
		finalstate=finalstates[seqnum++];			
		nsteps=0;
		pop = blank;
		pop(0,0)=1;
		while (Xseqs[gatenum]<55)
		{
			Htemp = Hdrift;
			MOs::Identity(unitary);		
			
			double ddt;
			ddt = rwa_dt/10.0;
			for(size_t j = 0; j < rwa_num_time; ++j){
			
			Udiag[0][0]=Udiag[1][0]=1;
			Udiag[0][1]=exp(-std::complex<double>(0.0,(double)j*qbfreq*dt));
			Udiag[0][2]=exp(-std::complex<double>(0.0,2.0*(double)j*qbfreq*dt));
			Udiag[0][3]=exp(-std::complex<double>(0.0,3.0*(double)j*qbfreq*dt));
			Udiag[1][1]=exp(std::complex<double>(0.0,(double)j*qbfreq*dt));
			Udiag[1][2]=exp(std::complex<double>(0.0,2.0*(double)j*qbfreq*dt));	
			Udiag[1][3]=exp(std::complex<double>(0.0,3.0*(double)j*qbfreq*dt));	
			
			Htemp = Hdrift;			  
			 for(size_t q=0; q< dim; ++q) {
				for(size_t p=0; p<dim; ++p){
						//Htemp(p,q) += Xamps[pulsenum]*Udiag[0][q]*Udiag[1][p]*(ctrl0[j]*HcontrolX(p,q)*cos(((double)j)*qbfreq*dt) + ctrl1[j]*HcontrolX(p,q)*sin(((double)j)*qbfreq*dt) );
						//Htemp(p,q) += Yamps[pulsenum]*Udiag[0][q]*Udiag[1][p]*(ctrl1[j]*HcontrolX(p,q)*cos(((double)j)*qbfreq*dt) - ctrl0[j]*HcontrolX(p,q)*sin(((double)j)*qbfreq*dt) );
						Htemp(p,q) += Xseqs[gatenum]*(ctrl0[j]*HcontrolX(p,q) + ctrl1[j]*HcontrolY(p,q));
						Htemp(p,q) += Yseqs[gatenum]*(-ctrl1[j]*HcontrolX(p,q) + ctrl0[j]*HcontrolY(p,q));
						//Htemp(q,p) += amplif*(ucontrollab0filt[j]*HcontrolX(q,p)+ ucontrollab1filt[j]*HcontrolY(q,p));
				}
			 }
			 unitary = ExpM::EigenMethod(Htemp,-std::complex<double>(0,ddt));		
			 nsteps++;		
			 for(size_t k = 0; k < 10; ++k){
				pop = unitary*pop*MOs::Dagger(unitary) + (ddt/1000)*(a*pop*ad -0.5*ad*a*pop -0.5*pop*ad*a);
				//pop +=  -std::complex<double>(0,ddt)*Htemp*pop + std::complex<double>(0,ddt)*pop*Htemp +  (ddt/1000)*(a*pop*ad -0.5*ad*a*pop -0.5*pop*ad*a);
			 }
			}
			//cout <<  " " << abs(pop(0,0) + pop(1,1) + pop(2,2) +pop(3,3))  << " " << ampX << " " << ampY << "\n";
			//return 0;
			nsteps+=2/rwa_dt;
			for(size_t j = 0; j < 2/ddt; ++j){
				pop += (ddt/1000)*(a*pop*ad -0.5*ad*a*pop -0.5*pop*ad*a);
			}			
			gatenum++;
		}
		gatenum++;
		if(abs(pop(finalstate,finalstate))<=0.5) pop(finalstate,finalstate)=0.5000000000000001;
		//cout << endl;
		xavg+=nsteps*rwa_dt;
		x2avg+=(nsteps*rwa_dt)*(nsteps*rwa_dt);
		lnyavg+=log(abs(pop(finalstate,finalstate))-0.5);
		xlnyavg+=(nsteps*rwa_dt)*log(abs(pop(finalstate,finalstate))-0.5);
		xyavg+=(nsteps*rwa_dt)*(abs(pop(finalstate,finalstate))-0.5);
		x2yavg+=(nsteps*rwa_dt)*(nsteps*rwa_dt)*(abs(pop(finalstate,finalstate))-0.5);
		ylnyavg+=(abs(pop(finalstate,finalstate))-0.5)*log(abs(pop(finalstate,finalstate))-0.5);
		xylnyavg+=(nsteps*rwa_dt)*(abs(pop(finalstate,finalstate))-0.5)*log(abs(pop(finalstate,finalstate))-0.5);
		
		timebins[index%17] = nsteps*rwa_dt;
		countbins[index%17]++;
		fidelitybins[index%17] += abs(pop(finalstate,finalstate))-0.5;
		clifffidelitybins[index/17/8][index%17]+= abs(pop(finalstate,finalstate))-0.5;
		cliffcountbins[index/17/8][index%17]++;
		clifftimebins[index/17/8][index%17] = nsteps*rwa_dt;
		if(index%17==16)
			cout <<index<< " " << timebins[index%17] << '\t' << abs(pop(0,0)) << '\t' << abs(pop(1,1)) <<'\t' << abs(pop(2,2)) << '\t'<< abs(pop(3,3)) <<'\t'<<  finalstate << '\t' << 0.5+fidelitybins[index%17]/countbins[index%17] << endl; 
		index++;
		//dataout << nsteps*rwa_dt << '\t' << abs(pop(finalstate,finalstate)) << endl;
	}	
	for(size_t j=0;j<17;j++)
	{	xavglny += timebins[j]*lnyavg;
		yavgxylny+= fidelitybins[j]*xylnyavg;
		yavgx2y+= fidelitybins[j]*x2yavg;
	}
	xyavgylny=ylnyavg*xyavg;
	double b = (32*17*xlnyavg - xavglny)/(32*17*x2avg-xavg*xavg);
	double b2 = (yavgxylny-xyavgylny)/(yavgx2y-xyavg*xyavg);
	cout << b2*(tgate) << endl;
	cout << x2avg << " " << xavg << " " << xlnyavg << " " << xavglny << " " << lnyavg<< " " << b << endl;					
	cout << 	old_tgate << '\t' <<  abs(b2*(tgate)) << '\t' << abs(b*(tgate)) << '\t' <<  1.0-fidel << '\t' <<  1.0-sys->fidelity_ <<endl;
	dataout << 	old_tgate << '\t' <<  abs(b2*(tgate)) << '\t' << abs(b*(tgate)) << '\t' <<  1.0-fidel << '\t' <<  1.0-sys->fidelity_ <<endl;
	cout << b2*(tgate) << endl;
		
	delete(sys);
	//cout << "done\n";	
			
//	for(int j=0; j<17; j++)
//	{	//cout << timebins[j] << '\t' << fidelitybins[j]/countbins[j] << endl; 
//		//dataout << timebins[j] << '\t' << fidelitybins[j]/countbins[j] << endl;
//	
//	dataout << clifftimebins[0][j] << '\t' << 0.5+0.5*exp(b2*clifftimebins[0][j]) << '\t' << 0.5+fidelitybins[j]/countbins[j] << '\t' << 0.5+fidelitybins[j]/countbins[j] << '\t' << 0.5+clifffidelitybins[0][j]/cliffcountbins[0][j]<< '\t' << 0.5+clifffidelitybins[1][j]/cliffcountbins[1][j]<< '\t' << 0.5+clifffidelitybins[2][j]/cliffcountbins[2][j]<< '\t' << 0.5+clifffidelitybins[3][j]/cliffcountbins[3][j] <<  endl;
//
//	}
	
	}
	dataout.close();
	
/*	for(unsigned amplif=0; amplif<74; amplif+=1)
	//for(double amplif=0.0001; amplif<1.5; amplif+=0.05)
	{
		MOs::Identity(unitary);
		for(size_t j = 0; j < num_time; ++j){
			Udiag[0][0]=Udiag[1][0]=1;
			Udiag[0][1]=exp(-std::complex<double>(0.0,(double)j*qbfreq*dt));
			Udiag[0][2]=exp(-std::complex<double>(0.0,2.0*(double)j*qbfreq*dt));
			Udiag[0][3]=exp(-std::complex<double>(0.0,3.0*(double)j*qbfreq*dt));
			Udiag[1][1]=exp(std::complex<double>(0.0,(double)j*qbfreq*dt));
			Udiag[1][2]=exp(std::complex<double>(0.0,2.0*(double)j*qbfreq*dt));	
			Udiag[1][3]=exp(std::complex<double>(0.0,3.0*(double)j*qbfreq*dt));	
			
			Htemp = Hdrift;			  
			for(size_t q=0; q< dim; ++q) {
				for(size_t p=0; p<dim; ++p){
						Htemp(p,q) += prefacX1[amplif]*Udiag[0][q]*Udiag[1][p]*(ctrl0[j]*HcontrolX(p,q)*cos(((double)j)*qbfreq*dt) + ctrl1[j]*HcontrolX(p,q)*sin(((double)j)*qbfreq*dt) );
						Htemp(p,q) += prefacY1[amplif]*Udiag[0][q]*Udiag[1][p]*(ctrl1[j]*HcontrolX(p,q)*cos(((double)j)*qbfreq*dt) - ctrl0[j]*HcontrolX(p,q)*sin(((double)j)*qbfreq*dt) );
						//Htemp(q,p) += amplif*(ucontrollab0filt[j]*HcontrolX(q,p)+ ucontrollab1filt[j]*HcontrolY(q,p));
				}
			}
			unitary = ExpM::EigenMethod(Htemp,-std::complex<double>(0,dt))*unitary;
			
		}	
		
		for(size_t j = 0; j < num_time; ++j){
			Udiag[0][0]=Udiag[1][0]=1;
			Udiag[0][1]=exp(-std::complex<double>(0.0,(double)j*qbfreq*dt));
			Udiag[0][2]=exp(-std::complex<double>(0.0,2.0*(double)j*qbfreq*dt));
			Udiag[0][3]=exp(-std::complex<double>(0.0,3.0*(double)j*qbfreq*dt));
			Udiag[1][1]=exp(std::complex<double>(0.0,(double)j*qbfreq*dt));
			Udiag[1][2]=exp(std::complex<double>(0.0,2.0*(double)j*qbfreq*dt));	
			Udiag[1][3]=exp(std::complex<double>(0.0,3.0*(double)j*qbfreq*dt));	
			
			Htemp = Hdrift;			  
			for(size_t q=0; q< dim; ++q) {
				for(size_t p=0; p<dim; ++p){
						Htemp(p,q) += prefacX2[amplif]*Udiag[0][q]*Udiag[1][p]*(ctrl0[j]*HcontrolX(p,q)*cos(((double)j)*qbfreq*dt) + ctrl1[j]*HcontrolX(p,q)*sin(((double)j)*qbfreq*dt) );
						Htemp(p,q) += prefacY2[amplif]*Udiag[0][q]*Udiag[1][p]*(ctrl1[j]*HcontrolX(p,q)*cos(((double)j)*qbfreq*dt) - ctrl0[j]*HcontrolX(p,q)*sin(((double)j)*qbfreq*dt) );
						//Htemp(q,p) += 1.0*(ucontrollab0filt[j]*HcontrolX(q,p)+ ucontrollab1filt[j]*HcontrolY(q,p))/2.0;
				}
			}			
			unitary = ExpM::EigenMethod(Htemp,-std::complex<double>(0,dt))*unitary;
		}	
		
		//unitary = unitary * U_desired;
		//std::cout << amplif << " " << real(( unitary(0,0)+unitary(1,1) ) * conj( unitary(0,0)+unitary(1,1) )/4.0)<< "\n";
		//dataout << amplif << " " << real(( unitary(0,0)+unitary(1,1) ) * conj( unitary(0,0)+unitary(1,1) )/4.0)<< endl;
		//dataout << amplif << " " << abs(unitary(0,0)) << " " << abs(unitary(1,0)) << " " <<  abs(unitary(2,0)) << " " << endl;
		//std::cout << amplif << " " << real(( unitary(0,0)+unitary(1,1) ) * conj( unitary(0,0)+unitary(1,1) )/4.0) << " " << abs(unitary(0,0)) << " " << abs(unitary(1,0)) << " " <<  abs(unitary(2,0)) << " " << endl ;
		dataout << abs(unitary(0,0)*unitary(0,0)) << ' ' << abs(unitary(2,0)*unitary(2,0)) << endl;
		std::cout << abs(unitary(0,0)*unitary(0,0)) << ' ' ;
		
		
	}*/
//	for(size_t j =0; j < num_time; j++){
//			//dataout << j*dt << '\t' << ucontrollab0filt[j] << '\t' << ucontrollab1filt[j] << '\t' << ucontrollab2filt[j]<< '\t' << endcontrols0[j] << '\t' << endcontrols1[j]<< '\t' << endcontrols2[j] << endl;
//			dataout << j*dt << '\t' << ucontrol0RWC[j] << '\t' << ucontrol1RWC[j] << '\t' << ucontrol2RWC[j]<< '\t' << ucontrollab0[j+200] << '\t' << ucontrollab1[j+200]<< '\t' << ucontrollab2[j+200] << endl;
//	}
	for(size_t j =0; j < rwa_num_time; j++){
//		dataout << j*rwa_dt << '\t' << newcontrols0[j] << '\t' << newcontrols1[j] << std::endl; 
	}
	std::cout <<'\n';
	

	return 0;
}