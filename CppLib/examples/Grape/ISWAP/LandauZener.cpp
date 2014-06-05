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
	double qbfreq=5*M_PI, Delta=-4.00, lambda=sqrt(3/2), det;  //physical params
	double lamb2 =  1/sqrt(2);
	double Delta2 = Delta;
	double tolerance=0.00000000000000001, fidelity, base_a=1.7, epsilon=0.00005, tgate=8, dt;  //search params
	double rab = M_PI;
	size_t max_iter=10000;
	size_t dimQ=2;
	size_t dim=dimQ*dimQ;//dimQ*dimQ*dimL*dimL2*dimL3;
	double to_ms =1000.0/CLOCKS_PER_SEC;
	size_t buffer = 0;
	
	const size_t rwasub=20;
	size_t mult = 1;
	
	size_t rwa_num_time= tgate*rwasub, num_controls =1;
	
		double rwa_dt = tgate/double(rwa_num_time);		
		
		OptimizeEvolution baselinesys(dim, num_controls, rwa_num_time, rwa_dt, "LandauZener");
		baselinesys.SetNumericalParameters(fidelity=0.9999, base_a, epsilon, tolerance, max_iter);	

		cout << "declaring Hilbert space\n";
		
		matrix<complex<double> > HDrift(dimQ, dimQ), Hcouple(dimQ*dimQ, dimQ*dimQ), HdriftQtot(dimQ, dimQ), HdriftQ1(dimQ, dimQ),  HdriftQ2(dimQ, dimQ);
		matrix<complex<double> > IdentQ(dimQ, dimQ);
		MOs::Identity(IdentQ);
		matrix<complex<double> > o(dimQ,dimQ), od(dimQ,dimQ), nq(dimQ,dimQ);
		
//		MOs::Destroy(a); ad = MOs::Dagger(a); n=ad*a;
		MOs::Destroy(o); 
		//o(1,2)=o(2,1)=0;
		od = MOs::Dagger(o); nq=od*o;
		
		
				
		//Jay - frequency selecition		
		Hcouple =  -(M_PI/2/(tgate-2*buffer*rwa_dt*rwasub)) *( MOs::TensorProduct(o, od) + MOs::TensorProduct(od, o) );
		det = -2*2*M_PI;
		
		//Amira - amplitude modulation
	//	Hcouple = 0.0 * MOs::TensorProduct(IdentQ, IdentQ);
	//	det = 0.0;
		
	
		Delta = -5 * 2*M_PI;
		//Delta = 25.49/tgate + 0.25;
		HdriftQ1 = det*nq + Delta/2*( nq*nq - 1.0*nq);
		HdriftQ2 = Delta/2*( nq*nq - 1.0*nq);
				
		HdriftQtot =MOs::TensorProduct(HdriftQ1, IdentQ) + MOs::TensorProduct(IdentQ, HdriftQ2);
		
		
		HDrift = HdriftQtot  +  Hcouple;

		HDrift.SetOutputStyle(Matrix);	
		cout << "Hdrift \n" << HDrift << endl;
		baselinesys.SetHdrift(HDrift);
		
		
		matrix<complex<double> > U_desired(dimQ*dimQ,dimQ*dimQ);
		
		MOs::Identity(U_desired);
		//ISWAP GATE
		U_desired(dimQ,1)=complex<double>(0,1);//std::complex<double>(0,1/sqrt(2));
		U_desired(1,dimQ)=complex<double>(0,1);//std::complex<double>(0,1/sqrt(2));
		U_desired(1,1)=U_desired(dimQ,dimQ)=0;//1/sqrt(2);
		baselinesys.SetRhoDesired(U_desired);
		U_desired.SetOutputStyle(Matrix);   //NOTE: for dim=3 need to change fidelity from Phi4 to something like Phi4Sub2 for 2 qubits
		cout << "u desired: \n" << U_desired << endl;

	
		cout << "delcaring base controls\n";
		
		//Jay controls
		AnalyticControl control(MOs::TensorProduct(nq, IdentQ), rwa_dt, rwa_num_time, rwasub, &baselinesys, "zcontrol");
		control.setparams(&AnalyticControl::SquarePulse, -det, 0, 0, NULL, NULL);
		
		//Amira controls
//		AnalyticControl control( MOs::TensorProduct(o, od) + MOs::TensorProduct(od, o), rwa_dt, rwa_num_time, rwasub, &baselinesys, "zcontrol");
//		control.setparams(&AnalyticControl::SquarePulse, -(M_PI/2/(tgate-2*buffer*rwa_dt*rwasub)), 0, 0, NULL, NULL);


		//populate controls
		control.pInterpolate = &Control::GaussianFilter;//CubicInterpolate;//
		control.pSubGradients = &Control::GaussianFilterGradient;//setLinearSplineGradients;//
		control.init();
		
		cout << Delta << " " << Delta/2/M_PI << " is delta\n";	
		
/* rotating frame stuff -- broken / not very useful
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
*/			
	
		cout << "beginning sweep...\n";
		baselinesys.UnitaryTransfer();
		
//		baselinesys.sweepinitialconditions();
		double min = 0.0, max = 2;
		HdriftQtot.resize(dim, dim);
		matrix<complex<double> > Shift = HdriftQtot;		
		matrix<complex<double> > Shift2 = HdriftQtot;
		
//		baselinesys.sweeptimes(static_cast<ptrPropagate> (&OptimizeEvolution::sweepinitialconditions));

		//baselinesys.sweepenergies(min, max, Shift, Shift2, HDrift, HDrift, static_cast<ptrPropagate> (&OptimizeEvolution::sweepinitialconditions));


		
		cout <<"done\n";		
		return 0;
	////////////////////////////
	
	
}