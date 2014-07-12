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
	double qbfreq=51.5, Delta=-2.53, lambda=sqrt(3/2);  //physical params
	double lamb2 =  1/sqrt(2);
	double rab = M_PI;
	double Delta2 = Delta;
	double tolerance=0.00000000000000000000000001, fidelity=0.9, base_a=1.5, epsilon=5, tgate=1.00, dt;  //search params
	size_t max_iter=100000;
	size_t dim=2;
	double to_ms =1000.0/CLOCKS_PER_SEC;
	int clo;
	
	//Typical operators	
	complex<double> ii(0.0,1.0);
	matrix<complex<double> > a(dim,dim), ad(dim,dim), n(dim,dim), UNOT(dim,dim), Ident(dim, dim), Hx1(dim*dim,dim*dim), Hy1(dim*dim,dim*dim), Hx2(dim*dim,dim*dim), Hy2(dim*dim,dim*dim);
	matrix<complex<double> > Urot(dim,dim), Udes(dim*dim,dim*dim);
	MOs::Destroy(a);
	ad = MOs::Dagger(a);
	n=ad*a;
	MOs::Identity(Ident);
	
	//RWA fixed parameters
	const size_t rwasub=40;
	size_t rwa_num_time=tgate*rwasub, num_controls =4;
	
	//RW fixed parameters
    //const int nsub=20;
	//size_t num_time=rwa_num_time*nsub/rwasub;
	//dt=tgate/double(num_time);
	
	system("mkdir hello");
	
		//RWA frame	
		double rwa_dt = tgate/double(rwa_num_time);
		GrapeUnitaryQubit sys(dim*dim, num_controls, rwa_num_time);
		
		sys.ptrSetHnonlinearity=&GrapeUnitary::SetOscillation; 
		sys.ptrSetHgradnonlin= &GrapeUnitary::SetGradOscillation;
		sys.SetPhysicalParameters(fidelity=0.0000000000000099999999999999999999, rwa_dt, base_a, epsilon, tolerance, max_iter, qbfreq, Delta);
		
		//sys.SetHdrift(Delta/2*(n*n-n));			//first transition
		matrix<complex<double> > HDrift1(dim, dim), HDrift2(dim,dim), HDrift(dim*dim,dim*dim);
		//level 0
		HDrift1 = 0.0*n;
		HDrift2 = Delta*n;
		MOs::GenPauliX(UNOT,0,1);		
		Udes=MOs::TensorProduct(UNOT, Ident);
		HDrift = MOs::TensorProduct(Ident, HDrift1) + MOs::TensorProduct(HDrift1, Ident);
		Hx2 =  MOs::TensorProduct(0.5*(a+ad), Ident) + MOs::TensorProduct(Ident, 0.5*(a+ad));
		Hx1 = MOs::TensorProduct(0.5*(a+ad), Ident) + MOs::TensorProduct(Ident, 0.5*(a+ad));
		Hy2 =  MOs::TensorProduct(Ident, 0.5*(-ii*a+ii*ad)) + MOs::TensorProduct( 0.5*(-ii*a+ii*ad), Ident); 
		Hy1 = MOs::TensorProduct(Ident, 0.5*(-ii*a+ii*ad)) + MOs::TensorProduct( 0.5*(-ii*a+ii*ad), Ident);
		
		Udes.SetOutputStyle(Matrix);
		Hx1.SetOutputStyle(Matrix);
		Hy1.SetOutputStyle(Matrix);
		HDrift.SetOutputStyle(Matrix);
		cout << Udes << "\n";
		cout << Hx1 << "\n" << Hy1 << "\n" << HDrift << endl;
//change lambda here
//		a(1,2)=a(2,1)=1.75;
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
	
		sys.SetRhoDesired(Udes);	
		
		matrix<complex<double> > Shift(dim, dim);
		MOs::Null(Shift);
		matrix<complex<double> > Shift2(dim, dim);
		MOs::Null(Shift2);
			
		Shift(0,1) = Shift(1,0) =1;
		Shift2(0,1) =-ii; 
		Shift2(1,0) =ii;
		double freq0[] = {0,Delta,0,Delta};
		double freq1[] = {0,-Delta,-2*Delta,-3*Delta};
					
		//Gaussian
		Control u0(Hx1, rwa_dt,rwa_num_time,rwasub, &Control::ShiftedGaussian, rab, 2.5, 1),
				u1(0.0*Hy1,rwa_dt, rwa_num_time, rwasub, &Control::Derivative, 1/Delta/4, 0, 0), 
				u2(0.0*Hx2,rwa_dt,rwa_num_time,rwasub, &Control::StarkShiftSelect, -1/8.0/Delta/Delta, 1.5, 0),
				u3(0.0*Hy2,rwa_dt,rwa_num_time,rwasub, &Control::Derivative, -1.333/Delta, 0, 0);
		u0.freqs_ = freq0; u1.freqs_ = freq0; u2.freqs_ = freq1; u3.freqs_ = freq1;	
		u0.pSetHcontrol = &Control::oscillatoryControl;	u1.pSetHcontrol = &Control::oscillatoryControl;	u2.pSetHcontrol = &Control::oscillatoryControl;	u3.pSetHcontrol = &Control::oscillatoryControl;	
		sys.Setucontrol(&u0, 0); //0th control
		sys.Setucontrol(&u1,1); //1st control
		sys.Setucontrol(&u2,2); //2nd control
		sys.Setucontrol(&u3,3); //3rd control
//		sys.sweeptimes( "sel_gaussian.dat", (&GrapeUnitary::Phi4), ( &GrapeUnitary::GradPhi4) );
//		u0.Hcontrol_(1,0) = u0.Hcontrol_(0,1) = u1.Hcontrol_(1,0) = u1.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0.Hcontrol_, u1.Hcontrol_, "lamb5nsHarm2.5_G3.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );

//1Drag
		Control u0C(Hx1, rwa_dt,rwa_num_time,rwasub, &Control::StarkShiftSelect, -1/8.0/Delta/Delta, 1.5, 0), 
				u1C(Hy1,rwa_dt, rwa_num_time, rwasub, &Control::Derivative, 1/Delta/4, 0, 0),
				u2C(Hx2,rwa_dt,rwa_num_time,rwasub, &Control::StarkShiftSelect, -1/8.0/Delta/Delta, 1.5, 0),
				u3C(Hy2,rwa_dt,rwa_num_time,rwasub, &Control::Derivative, -1.5/Delta, 0, 0);
		u0C.freqs_ = u1C.freqs_ = freq0; u2C.freqs_ = freq1; u3C.freqs_ = freq1;		
		u0C.pSetHcontrol = &Control::oscillatoryControl;	u1C.pSetHcontrol = &Control::oscillatoryControl;	u2C.pSetHcontrol = &Control::oscillatoryControl;	u3C.pSetHcontrol = &Control::oscillatoryControl;	
		sys.Setucontrol(&u0C, 0); //0th control
		sys.Setucontrol(&u1C,1); //1st control
		sys.Setucontrol(&u2C,2); //2nd control
		sys.Setucontrol(&u3C,3); //3rd control
		u0C.pInterpolate = u1C.pInterpolate = u2C.pInterpolate = u3C.pInterpolate =  &Control::Void;//Pixelate;//
//		sys.sweeptimes( "sel_drag2.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4) );

//		u0C.Hcontrol_(1,0) = u0C.Hcontrol_(0,1) = u1C.Hcontrol_(1,0) = u1C.Hcontrol_(0,1) = 0;
//		sys.sweepenergies(0, 2, Shift, Shift2, u0C.Hcontrol_, u1C.Hcontrol_, "lamb5nsHarm2.5_drgC_G2.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4Sub2), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4Sub2) );
//u0C.writefile("gaussian.dat");
//u1C.writefile("gaussian1.dat");
//u2C.writefile("gaussian2.dat");
//u3C.writefile("gaussian3.dat");

		//grp
		Control u0G(Hx1, rwa_dt,rwa_num_time,rwasub, &Control::StarkShiftSelect, -1/8.0/Delta/Delta, 1.5, 0), 
				u1G(Hy1,rwa_dt, rwa_num_time, rwasub, &Control::Derivative, 1/Delta/4, 0, 0),
				u2G(Hx2,rwa_dt,rwa_num_time,rwasub, &Control::StarkShiftSelect, -1/8.0/Delta/Delta, 1.5, 0),
				u3G(Hy2,rwa_dt,rwa_num_time,rwasub, &Control::Derivative, -1.5/Delta, 0, 0);
		u0G.freqs_ = freq0; u1G.freqs_ = freq0; u2G.freqs_ = freq1; u3G.freqs_ = freq1;		
		u0G.pSetHcontrol = &Control::oscillatoryControl;	u1G.pSetHcontrol = &Control::oscillatoryControl;	u2G.pSetHcontrol = &Control::oscillatoryControl;	u3G.pSetHcontrol = &Control::oscillatoryControl;	
		u0G.pInterpolate = u1G.pInterpolate = u2G.pInterpolate = u3G.pInterpolate =  &Control::LinearInterpolate;//Pixelate;//
		u0G.pSubGradients = u1G.pSubGradients = u2G.pSubGradients = u3G.pSubGradients = &Control::setLinearSplineGradients;//setPixelsGradients;//
		sys.Setucontrol(&u0G, 0); //0th control
		sys.Setucontrol(&u1G,1); //1st control
		sys.Setucontrol(&u2G,2); //2nd control
		sys.Setucontrol(&u3G,3); //3rd control
		sys.SetPhysicalParameters(fidelity=0.99999999999999999999, rwa_dt, base_a, epsilon, tolerance, max_iter, qbfreq, Delta);
		sys.sweeptimes( "sel_grape.dat", static_cast< Fid >(&GrapeUnitaryQubit::Phi4), static_cast< GradFid >( &GrapeUnitaryQubit::GradPhi4) );
		
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
		cout << "time take: " << to_ms*(clock() - clo) << endl;
	
	return 0;
	////////////////////////////
	
	
}