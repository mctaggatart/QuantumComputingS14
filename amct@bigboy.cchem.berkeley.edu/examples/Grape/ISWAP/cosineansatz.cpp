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
	double qbfreq=5*M_PI, Delta=-4.00, lambda=sqrt(3/2);  //physical params
	double lamb2 =  1/sqrt(2);
	double Delta2 = Delta;
	double tolerance=0.00000000000000001, fidelity=0.9, base_a=1.7, epsilon=0.00005, tgate=8, dt;  //search params
	double rab = M_PI;
	size_t max_iter=1000;
	size_t dimQ=3, dimL=2, dimL2=2, dimL3=2;
	size_t dim=dimQ*dimQ;//dimQ*dimQ*dimL*dimL2*dimL3;
	size_t refdim=dimQ*dimQ*dimQ*dimQ*dimQ;
	double to_ms =1000.0/CLOCKS_PER_SEC;
	size_t buffer = 0;
	
	//Typical operators	
	complex<double> ii(0.0,1.0);
	matrix<complex<double> > a(refdim,refdim), ad(refdim,refdim), n(refdim,refdim), Ident(refdim, refdim), twophot(refdim, refdim), UNOT(refdim,refdim) ;
	MOs::Destroy(a);
	ad = MOs::Dagger(a);
	n=ad*a;
	MOs::Identity(Ident);
	MOs::GenPauliX(UNOT,0,1);
	MOs::Null(twophot);
	twophot(0,2)=1;//-ii;
	twophot(2,0)=1;//ii;

	//RWA fixed parameters
	
	const size_t rwasub=20;
	size_t mult = 1;
	
	size_t rwa_num_time= tgate*rwasub*4, num_controls =3;
	
	//RW fixed parameters
    //const int nsub=20;
	//size_t num_time=rwa_num_time*nsub/rwasub;
	//dt=tgate/double(num_time);
	
//	system("mkdir hello");
	
		//RWA frame	
		double rwa_dt = tgate/double(rwa_num_time);		
		
		OptimizeEvolution baselinesys(dim, num_controls, rwa_num_time, rwa_dt, "LandauZener");
		baselinesys.Phi = &OptimizeEvolution::Phi4;
		baselinesys.gradPhi = &OptimizeEvolution::GradPhi4;
		baselinesys.pPropagate = &Evolution::forwardpropagate;
		baselinesys.pGetGradient = 	&OptimizeEvolution::computegradient;
		
		cout << "declaring Hilbert space\n";
		
		matrix<complex<double> > HDriftRef(refdim, refdim);
		//level 0
		
		matrix<complex<double> > HDrift(dimQ, dimQ), Hcouple(dimQ*dimQ, dimQ*dimQ), HdriftQ(dimQ, dimQ),  HdriftQ2(dimQ, dimQ);
		HDrift.SetOutputStyle(Matrix);
		//level 0
		//HDrift = Delta/2*( n*n - n);
		
		matrix<complex<double> > IdentQ(dimQ, dimQ), IdentQ2(dimL, dimL), IdentQ3(dimL2, dimL2), IdentQ4(dimL3, dimL3);
		MOs::Identity(IdentQ); IdentQ2=MOs::TensorProduct(IdentQ,IdentQ); IdentQ3=MOs::TensorProduct( IdentQ,IdentQ2); IdentQ4=MOs::TensorProduct( IdentQ,IdentQ3);
		matrix<complex<double> > o(dimQ,dimQ), od(dimQ,dimQ), nq(dimQ,dimQ);
		
//		MOs::Destroy(a); ad = MOs::Dagger(a); n=ad*a;
		MOs::Destroy(o); 
		//o(1,2)=o(2,1)=0;
		od = MOs::Dagger(o); nq=od*o;
		
		double det = 0*2*2*M_PI;
		double det2 = 0.7*2*M_PI;
		double det3 = 0.8*2*M_PI;
		double det4 = 0.9*2*M_PI;
		Delta = -1 * 0.15 * 2*M_PI;
		//Delta = 25.49/tgate + 0.25;
		
		HdriftQ2 = det*nq + Delta/2*( nq*nq - 1.0*nq);
		HdriftQ = MOs::TensorProduct(IdentQ3, MOs::TensorProduct(HdriftQ2, IdentQ) + MOs::TensorProduct(IdentQ, HdriftQ2));
		Hcouple =  -1.0*( MOs::TensorProduct(o, od) + MOs::TensorProduct(od, o) );
		HDriftRef = HdriftQ  + MOs::TensorProduct(IdentQ3, (M_PI/2/(tgate-2*buffer*rwa_dt*rwasub)) * Hcouple)
				+  1.0*(M_PI/2/(tgate-2*buffer)) * MOs::TensorProduct(IdentQ2, MOs::TensorProduct( Hcouple , IdentQ) )
				+  1.1*(M_PI/2/(tgate-2*buffer)) *MOs::TensorProduct(IdentQ2, MOs::TensorProduct( o,MOs::TensorProduct(IdentQ,od)) + MOs::TensorProduct( od,MOs::TensorProduct(IdentQ,o)))
				+  0.8*(M_PI/2/(tgate-2*buffer)) *MOs::TensorProduct(IdentQ, ( MOs::TensorProduct(o, MOs::TensorProduct( IdentQ2,od)) + MOs::TensorProduct(od, MOs::TensorProduct( IdentQ2,o))))
				+  0.9*(M_PI/2/(tgate-2*buffer)) *MOs::TensorProduct(IdentQ, ( MOs::TensorProduct(Hcouple,IdentQ2)))
				+  0.7*(M_PI/2/(tgate-2*buffer)) *MOs::TensorProduct(IdentQ, ( MOs::TensorProduct(o, MOs::TensorProduct( IdentQ,MOs::TensorProduct(od, IdentQ))) + MOs::TensorProduct(od, MOs::TensorProduct( IdentQ,MOs::TensorProduct(o, IdentQ)))))
				
				+  0.5*(M_PI/2/(tgate-2*buffer)) * MOs::TensorProduct(Hcouple, IdentQ3) 
		//		+   1.1*(M_PI/2/(tgate-2*buffer)) *(MOs::TensorProduct(od,( MOs::TensorProduct(IdentQ, MOs::TensorProduct( o,MOs::TensorProduct(IdentQ,IdentQ))))) 
		//								+   MOs::TensorProduct(o, MOs::TensorProduct(IdentQ, MOs::TensorProduct( od,MOs::TensorProduct(IdentQ,IdentQ)))) )
				+  0.85*(M_PI/2/(tgate-2*buffer)) *(MOs::TensorProduct(o, ( MOs::TensorProduct(IdentQ, MOs::TensorProduct( IdentQ,MOs::TensorProduct(od, IdentQ))))) 
										+ MOs::TensorProduct(od,   MOs::TensorProduct(IdentQ, MOs::TensorProduct( IdentQ,MOs::TensorProduct(o, IdentQ)))))
				;

		HDriftRef.SetOutputStyle(Matrix);	


								
//		cout << HDriftRef << endl;
		HDrift = HDriftRef;
		HDrift.resize(dim, dim);
		
		HDrift.SetOutputStyle(Matrix);	
		cout << "Hdrift \n" << HDrift << endl;
		baselinesys.SetHdrift(HDrift);
		
		
		matrix<complex<double> > U_desired(dimQ*dimQ,dimQ*dimQ);
		
		MOs::Identity(U_desired);
		//ISWAP GATE
		U_desired(dimQ,1)=complex<double>(0,1);//std::complex<double>(0,1/sqrt(2));
		U_desired(1,dimQ)=complex<double>(0,1);//std::complex<double>(0,1/sqrt(2));
		U_desired(1,1)=U_desired(dimQ,dimQ)=0;//1/sqrt(2);
		
		//CNOT GATE
//		U_desired(dimQ,dimQ+1)=complex<double>(1,0);//std::complex<double>(0,1/sqrt(2));
//		U_desired(dimQ+1,dimQ)=complex<double>(1,0);//std::complex<double>(0,1/sqrt(2));
//		U_desired(dimQ,dimQ)=U_desired(dimQ+1,dimQ+1)=0;//1/sqrt(2);
		
		matrix<complex<double> > U_desired_big(refdim, refdim);
		U_desired_big =  MOs::TensorProduct(IdentQ, MOs::TensorProduct(IdentQ,MOs::TensorProduct(IdentQ, U_desired)));
		U_desired.SetOutputStyle(Matrix);
		

		ptrSetTransition pSetTransition = &AnalyticControl::setdriving01;
		
	
		baselinesys.SetNumericalParameters(fidelity=0.9999, base_a, epsilon, tolerance, max_iter);	
		U_desired =U_desired_big;
		U_desired.resize(dim,dim);
		baselinesys.SetRhoDesired(U_desired);
		cout << "u desired: \n" << U_desired << endl;

		cout << "declaring controls...\n";
		
		ad = MOs::Dagger(a);
		AnalyticControl refcontrol(0.5*(o+od), rwa_dt, rwa_num_time, rwa_num_time, NULL, "refcontrol");
		refcontrol.setparams(&AnalyticControl::SquarePulse, det, 0, 0);
		refcontrol.pInterpolate =  &AnalyticControl::Void;
		refcontrol.setparams(&AnalyticControl::SquarePulse, 1, 0, 0, pSetTransition, NULL);
	
		cout << "delcaring base controls\n";
		
		Control *refdimcon[num_controls];
		AnalyticControl *baseControl[num_controls];
		char filenames[num_controls][80];
		matrix<complex<double> > *Hcontrolref;
		for(int k=0; k<num_controls; k++)
		{	refdimcon[k] = new Control(refcontrol, rwasub);
			Hcontrolref = new matrix<complex<double> >(1,1);
			(*Hcontrolref)(0,0) = 1;
			for(int j=0; j<num_controls; j++)
				if(j!=k) *Hcontrolref = MOs::TensorProduct(IdentQ, *Hcontrolref);
				else *Hcontrolref = MOs::TensorProduct(nq, *Hcontrolref);
			refdimcon[k]->Hcontrol_ = *Hcontrolref;
			if(k==2) refdimcon[k]->Hcontrol_=0.0*refdimcon[k]->Hcontrol_;
			Hcontrolref->resize(dim,dim);
			strcpy(filenames[k],"tempcontrol_ ");
			filenames[k][12] = '0'+k;
			//baseControl[k] = new AnalyticControl(Hcouple /* *Hcontrolref*/ , rwa_dt, rwa_num_time,rwasub,&baselinesys, filenames[k]);
			baseControl[k] = new AnalyticControl(Hcouple /* *Hcontrolref*/ , rwa_dt, rwa_num_time,rwa_num_time,&baselinesys, filenames[k]);  //flat controls
			baseControl[k]->pInterpolate = &Control::Pixelate;//GaussianFilter;//CubicInterpolate;//
			baseControl[k]->pSubGradients = &Control::setPixelsGradients; //GaussianFilterGradient;//setLinearSplineGradients;//
			baseControl[k]->buffer = buffer;	
			baseControl[k]->setparams(&AnalyticControl::scanamp, 0, 0,  0, pSetTransition, &refcontrol);
		}

		baseControl[2]->pSetHcontrol = &Control::drivefreqControl;	
		baseControl[2]->pSetHgradient = &Control::drivefreqControlGradient;	
		baseControl[1]->pSetHcontrol = &Control::fmdrivenControl;	
		baseControl[1]->pSetHgradient = &Control::fmdrivenControlGradient;	
		baseControl[1]->framefreq_=0;
		baseControl[1]->drivefreq_=Delta;
	//	baseControl[1]->Hcontrol_ = (M_PI/2/(tgate-2*buffer*rwa_dt*rwasub))*Hcouple; // baseControl[1]->Hcontrol_;
		baseControl[0]->Hcontrol_ = baseControl[1]->Hcontrol_; // baseControl[1]->Hcontrol_;
//		baseControl[2]->pSetHgradient = &Control::fmdrivenControlGradient;	
	//	double SS = -(M_PI/2/(tgate-2*buffer*rwa_dt*rwasub))*(M_PI/2/(tgate-2*buffer*rwa_dt*rwasub))/Delta;
	//	cout << " Stark Shift " << SS << endl;
		//baseControl[1]->setparams(&AnalyticControl::SquarePulse, det, 0, 0, pSetTransition, &refcontrol);
//		baseControl[0]->setparams(&AnalyticControl::SquarePulse, -Delta*0.2398, 0, 0, pSetTransition, &refcontrol);
		baseControl[2]->setparams(&AnalyticControl::scanamp, -Delta/4/2, -3*Delta/4*2,  0, pSetTransition, &refcontrol);
		baseControl[1]->setparams(&AnalyticControl::scanamp, 0, (M_PI/2/(tgate-2*buffer*rwa_dt*rwasub))*3,  0, pSetTransition, &refcontrol);
		baseControl[0]->setparams(&AnalyticControl::scanamp, 0, (M_PI/2/(tgate-2*buffer*rwa_dt*rwasub))*3,  0, pSetTransition, &refcontrol);
		baseControl[2]->ampctrl = baseControl[1];
		baseControl[1]->detctrl = baseControl[2];
		//baseControl[1]->setparams(&AnalyticControl::Sine, 0.3, abs(Delta*1+0*SS), 0, pSetTransition, &refcontrol);
		//baseControl[1]->setparams(&AnalyticControl::RandomControl, 0.5, -0.5, 0, pSetTransition, &refcontrol);
		cout << Delta << " " << Delta/2/M_PI << " is delta\n";	
		
		double freqs[] = {0, 0, det, det, det2, det2, det+det2, det+det2, det3, det3, det+det3, det+det3, det2+det3, det2+det3, det+det2+det3, det+det2+det3,
						   det4, det4, det+det4, det+det4, det2+det4, det2+det4, det+det2+det4, det+det2+det4, det3+det4, det3+det4, det+det3+det4, det+det3+det4, det2+det3+det4, det2+det3+det4, det+det2+det3+det4, det+det2+det3+det4};
//		double freqs[] = {0, 0, 0, det, det, det, 2*det, 2*det, det2, det2, det2, det+det2, det+det2, det+det2, det3, det3, det+det3, det+det3, det2+det3, det2+det3, det+det2+det3, det+det2+det3,
//						   det4, det4, det+det4, det+det4, det2+det4, det2+det4, det+det2+det4, det+det2+det4, det3+det4, det3+det4, det+det3+det4, det+det3+det4, det2+det3+det4, det2+det3+det4, det+det2+det3+det4, det+det2+det3+det4};

		baselinesys.freqs_ = NULL; //freqs;
//		u3C.freqs_ = u2C.freqs_ = u0C.freqs_ = u1C.freqs_ = freqs;
//		u3C.drivefreq_ = u2C.drivefreq_ = u0C.drivefreq_ =u1C.drivefreq_ = det; 
//		u3C.framefreq_ = u2C.framefreq_ = u0C.framefreq_ =u1C.framefreq_ = 0*det; 
//		u3C.pSetHcontrol = u2C.pSetHcontrol = u1C.pSetHcontrol = u0C.pSetHcontrol =  &AnalyticControl::drivenControl;
//		u3C.pSetHgradient = u2C.pSetHgradient = u1C.pSetHgradient = u0C.pSetHgradient = &AnalyticControl::drivenControlGradient;
//		u3C.relphase_ = u1C.relphase_ = M_PI/2;
		
		cout << "creating alternate evolutions...\n";
		OptimizeEvolution **sys = new OptimizeEvolution*[2];
		baselinesys.evolutions = sys;	
		
		matrix<complex<double> > HDriftS(dimQ, dimQ);
		HDriftS = HDrift;
		
		for (size_t s=0; s<2; s++)
		{	sys[s] = new OptimizeEvolution(dim, num_controls, rwa_num_time/mult, rwa_dt*mult, "LandauZener");
			sys[s]->SetNumericalParameters(fidelity=0.0000999, base_a, epsilon, tolerance, max_iter);
			sys[s]->Phi = baselinesys.Phi;
			sys[s]->gradPhi = baselinesys.gradPhi;
			sys[s]->pPropagate = &Evolution::forwardpropagate;
			sys[s]->pGetGradient = 	&OptimizeEvolution::computegradient;
//			HDriftS(2,2) = HDriftS(5,5) = HDriftS(6,6) = HDriftS(7,7) = -0.330 * 2 * M_PI;
//			HDriftS(8,8) = -0.330 * 4 * M_PI;
			sys[s]->SetHdrift(HDriftS);	
			cout << HDriftS << endl;
			sys[s]->freqs_ = baselinesys.freqs_;							
		}

		sys[1]->pPropagate = &Evolution::forwardpropagate_Commute;	
		sys[1]->pGetGradient = 	&OptimizeEvolution::computegradient_Commute;			
		
		cout << "duplicating controls...\n";
		AnalyticControl *sysControls[2][num_controls];
		
		for(int k=0; k<num_controls; k++)
		{	for (size_t s=0; s<2; s++)
			{
				sysControls[s][k] = new AnalyticControl(*(AnalyticControl*)baselinesys.GetuControl(k), rwasub/mult);				
				sysControls[s][k]->setparams(&AnalyticControl::ShiftedGaussian,rab,2,1,pSetTransition,(AnalyticControl*)baselinesys.GetuControl(k));
				
				sysControls[s][k]->pInterpolate = &AnalyticControl::CubicInterpolate;//Pixelate; //Void;//LinearInterpolate;//:
				sysControls[s][k]->pSubGradients = &AnalyticControl::setLinearSplineGradients;//setSubVoidGradients;//setPixelsGradients; //
				
				sysControls[s][k]->pSetHcontrol =  &Control::linearControl;
				sysControls[s][k]->pSetHgradient = &Control::linearGradient;
				
				sysControls[s][k]->evol = sys[s];
				if(s) strcpy(filenames[k],"sampled ");
				else  strcpy(filenames[k],"splited ");
				filenames[k][7] = '0'+k;
				strcpy(sysControls[s][k]->filename,filenames[k]);
				sys[s]->Setucontrol(sysControls[s][k]);
			}						
		}
	
		cout << "beginning sweep...\n";
//		baselinesys.sweepinitialconditions();
		double min = 0.0, max = 2;
		HdriftQ.resize(dim, dim);
		matrix<complex<double> > Shift = HdriftQ;		
		matrix<complex<double> > Shift2 = HdriftQ;
		
//		baselinesys.sweeptimes(static_cast<ptrPropagate> (&OptimizeEvolution::sweepinitialconditions));
		baselinesys.sweepenergies(min, max, Shift, Shift2, HDrift, HDrift, static_cast<ptrPropagate> (&OptimizeEvolution::sweepinitialconditions));
//		baselinesys.sweeptimesandcompare(static_cast<ptrPropagate> (&OptimizeEvolution::UnitaryTransfer));
//		baselinesys.sweepdimandcompare(	static_cast<ptrPropagate> (&OptimizeEvolution::UnitaryTransfer), HDriftRef, U_desired_big, refdimcon);

		cout <<"deleting controls\n";		
		for (size_t s=0; s<2; s++)
		for(int k=0; k<num_controls; k++)
			delete sysControls[s][k];
		for(int k=0; k<num_controls; k++)
			delete refdimcon[k];

		cout <<"deleting alternate evolutions\n";		
		for (size_t s=0; s<2; s++) delete sys[s];
		delete [] sys;
		
		cout <<"done\n";		
		return 0;
	////////////////////////////
	
	
}