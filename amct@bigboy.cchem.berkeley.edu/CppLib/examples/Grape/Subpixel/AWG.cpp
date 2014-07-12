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
	double tolerance=0.00000000000000001, fidelity, base_a=1.7, epsilon=0.00005, tgate=20, dt;  //search params
	double rab = M_PI;
	size_t max_iter=500000;
	size_t dimQ=2, dimL=2, dimL2=2, dimL3=2;
	size_t dim=dimQ*dimQ;//dimQ*dimQ*dimL*dimL2*dimL3;
	size_t refdim=dimQ*dimQ*dimQ*dimQ*dimQ;
	double to_ms =1000.0/CLOCKS_PER_SEC;
	size_t buffer = 1;

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
	
	const size_t rwasub=600;
	size_t mult = 1;
		
	size_t subpix = 2000/5*1;
	Control ctrl(a, 1.0/subpix, subpix*(5+2*buffer), subpix, NULL, "T6_pix.dat");
	ctrl.buffer = buffer;
	ctrl.u_[3*subpix/2]=0.15;
	ctrl.u_[5*subpix/2]=0.10;
	ctrl.u_[7*subpix/2]=0.05;
	ctrl.u_[9*subpix/2]=0.176;
	ctrl.u_[11*subpix/2]=-0.05;
	
	ctrl.pInterpolate =  &AnalyticControl::Pixelate;//
	ctrl.Interpolate();	
	
	ctrl.writefile();
	
	ctrl.pInterpolate =  &AnalyticControl::GaussianFilter; //CubicInterpolate;//Pixelate;//

	ctrl.Interpolate();	
	strcpy(ctrl.filename, "T6_filter.dat");
	ctrl.writefile();

	delete [] ctrl.u_filt;
	ctrl.u_filt = ctrl.u_;

//	ctrl.u_[subpix/2]=0.12;
//	ctrl.u_[3*subpix/2]=0.08;
//	ctrl.u_[5*subpix/2]=0.04;
//	ctrl.u_[7*subpix/2]=0.14;
//	ctrl.u_[9*subpix/2]=-0.04;
	ctrl.u_[3*subpix/2]=0.15;
	ctrl.u_[5*subpix/2]=0.10;
	ctrl.u_[7*subpix/2]=0.05;
	ctrl.u_[9*subpix/2]=0.176;
	ctrl.u_[11*subpix/2]=-0.05;
//	ctrl.u_[subpix/2]=0.15*0.5+0.05;
//	ctrl.u_[3*subpix/2]=0.10*0.5+0.05;
//	ctrl.u_[5*subpix/2]=0.05*0.5+0.05;
//	ctrl.u_[7*subpix/2]=0.176*0.5+0.05;
//	ctrl.u_[9*subpix/2]=-0.05*0.5+0.05;

	
	strcpy(ctrl.filename, "T6_spline.dat");
	ctrl.pInterpolate =  &AnalyticControl::Pixelate;//CubicInterpolate;//Pixelate;//		
	ctrl.Interpolate();	
	ctrl.pInterpolate =  &AnalyticControl::CubicInterpolate;//CubicInterpolate;//Pixelate;//		
	ctrl.Interpolate();	
	
	ctrl.writefile();
	/**/

	Control ctrl1(a, 1.0/subpix, subpix*(5+2*buffer), subpix, NULL, "T6.dat");
	ctrl1.buffer = buffer;
	ctrl1.readfile("T6.txt");
	ctrl1.writefile();


	return 0;

	
}