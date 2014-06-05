/*
Name: GrapeCabity.cpp
Author: fmotzoi

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
#include <GrapeUnitaryCavity.hpp>
using namespace std;

int main (int argc, char const *argv[]){

	verbose=no;
	cout << "Running program " << argv[0] << endl;

	//Grape inputs
	double wa=8.1*2*M_PI;  //9.2
	double wr = 10.7*2*M_PI; 		//g2=0.12 @7GHz
	double DeltaW=2*2*M_PI, g = 0.190*2*M_PI;
	size_t max_iter=100000;
	complex<double> ii(0.0,1.0);
	double tolerance=0.0000000000000000000001, fidelity, base_a=2.0, epsilon=1, sigma=1, tgate=5, dt;
	size_t num_time=tgate*10, dimQ=2, dimCav=2, dim= dimQ*dimCav, num_controls =4;
	dt=tgate/double(num_time);
	cout << dim << endl;
	GrapeUnitaryCavity sys(dim, num_controls, num_time);
	sys.SetPhysicalParameters(fidelity=0.00000000999, dt, base_a, epsilon, tolerance, max_iter, wa, DeltaW);
	
	matrix<complex<double> > HcontrolZ(dim,dim), HcontrolZ2(dim,dim), HcontrolX(dim,dim), HcontrolY(dim,dim);
	matrix<complex<double> > Urot(dim,dim);
	matrix<complex<double> > U_desired(dim,dim);
	U_desired.SetOutputStyle(Matrix);
	MOs::Identity(U_desired);
	
	matrix<complex<double> > Hdrift(dim, dim), Hcouple(dim, dim), HdriftQ(dimQ, dimQ),  HdriftCav(dimCav, dimCav);  
	matrix<complex<double> > IdentQ(dimQ, dimQ), IdentCav(dimCav, dimCav);
	MOs::Identity(IdentQ); MOs::Identity(IdentCav);
	
	matrix<complex<double> > a(dimCav,dimCav), ad(dimCav,dimCav), n(dimCav,dimCav);
	matrix<complex<double> > o(dimQ,dimQ), od(dimQ,dimQ), nq(dimQ,dimQ);
	
	Hdrift.SetOutputStyle(Matrix); HdriftCav.SetOutputStyle(Matrix); HdriftQ.SetOutputStyle(Matrix);
	HcontrolX.SetOutputStyle(Matrix);HcontrolY.SetOutputStyle(Matrix);HcontrolZ.SetOutputStyle(Matrix);HcontrolZ2.SetOutputStyle(Matrix);
	
	MOs::Destroy(a); ad = MOs::Dagger(a); n=ad*a;
	MOs::Destroy(o); od = MOs::Dagger(o); nq=od*o;
	
	//double wa = -Delta/2;
	//double wr = Delta/2;
		
	MOs::Destroy(HdriftQ);
	MOs::Destroy(HdriftCav);
	HdriftCav = wr*n;
	double Delta = -0.330 * 2 * M_PI;
	HdriftQ = wa*nq + Delta/2*( nq*nq - 1.0*nq);
	Hcouple = g * ( MOs::TensorProduct(a, od) + MOs::TensorProduct(ad, o) );
	Hdrift = MOs::TensorProduct(HdriftCav, IdentQ) + MOs::TensorProduct(IdentCav, HdriftQ) + Hcouple;
	cout << Hdrift << endl;
	
	double vals[dim];
	
	HcontrolX = MOs::TensorProduct(IdentCav, 0.5*(o + od));
	HcontrolY = MOs::TensorProduct(IdentCav, 0.5*(-ii*o+ii*od));//+MOs::TensorProduct(0.5*(-ii*a+ii*ad), IdentQ);
	HcontrolZ = (MOs::TensorProduct(0.5*n, IdentQ) + MOs::TensorProduct(IdentCav, 0.5*nq));
	//HcontrolZ = (MOs::TensorProduct(n-0.5*IdentCav, IdentQ) + MOs::TensorProduct(IdentCav,nq- 0.5*IdentQ));
	HcontrolZ2 = (MOs::TensorProduct(0.5*n, IdentQ) + MOs::TensorProduct(IdentCav, nq*0.5));

	//HcontrolY = Hcouple; //MOs::TensorProduct(0.5*(a + ad), IdentQ);
	//cout << HcontrolZ << endl;
	
	// HcontrolX(4,5) = HcontrolX(5,4) = 1;
	
	
	QOs::MatrixElements(Hdrift, dim, vals, HcontrolX);
	if(dimQ==3) HcontrolX(3,2) = HcontrolX(2,3) = 0;
	cout << "...\n";
	QOs::MatrixElements(Hdrift, dim, vals, HcontrolY);
	cout << "...\n";
	QOs::MatrixElements(Hdrift, dim, vals, Hdrift);
//	HcontrolX(1,2) = HcontrolX(2,1) = 0.0;
//	
//	HcontrolX(2,5) = HcontrolX(5,2) = 0;
//	HcontrolX(4,5) = HcontrolX(5,4) = 0.5;
//	HcontrolX(1,0) = HcontrolX(0,1) = 0.25;
//	HcontrolX(4,3) = HcontrolX(3,4) = 0;
	
//	HcontrolY(1,2) = HcontrolY(2,1) = 0.0;
//	HcontrolY(3,2) = HcontrolY(2,3) = 0;
//	HcontrolY(2,5) = HcontrolY(5,2) = 0;
//	HcontrolY(3,0) = HcontrolY(0,3) = 0;
//	HcontrolY(1,4) = HcontrolY(4,1) = 0;
//	HcontrolY(4,5) = HcontrolY(5,4) = 0;
//	HcontrolY(1,0) = HcontrolY(0,1) = 0;
//	HcontrolY(4,3) = HcontrolY(3,4) = 0;

	HcontrolY(2,0) = 0;
	HcontrolY(0,2) = 0;
	HcontrolY(1,3) = 0;
	HcontrolY(3,1) = 0;
	
	cout << HcontrolY << endl;
	complex<double> wd = (Hdrift(dimQ+1,dimQ+1) - Hdrift(0,0))*0.5;
	Hdrift = Hdrift - wd*( MOs::TensorProduct(n, IdentQ) + MOs::TensorProduct(IdentCav, nq) );

	double lambda1 = real(HcontrolX(dimQ,dimQ+1));
	double lambda2 = real(HcontrolX(0,1));
	double lambda3 = real(HcontrolX(0,dimQ));
	double lambda4 = real(HcontrolX(1,dimQ+1));
	double delta1 = real(Hdrift(1,1));
	double delta2 = real(Hdrift(dimQ,dimQ));
	double delta3 = 1;

	double lambda7=0, lambda8=0;
	double delta4 = real(Hdrift(2,2));
	
	double lambda5, lambda6, ddelta;
	if(dimQ==2){
		 lambda5 = 0;	
		 delta3 =10000000;
		 delta4 =10000000;
		 lambda6 = 0;	
		 ddelta=10000000;
	} else {
		 lambda5 = real(HcontrolX(dimQ+1,dimQ+2));	
		 delta3 = real(Hdrift(dimQ+2,dimQ+2));
		 lambda6 = real(HcontrolX(1,2));
		 ddelta = delta4 - delta1;
	}
	
	double rab_amp_th = -lambda1*lambda3/delta2 - lambda2*lambda4/delta1;
	
	double rab_amp2 =  lambda1*lambda1*lambda1*lambda3/6/delta2/delta2/delta2 + lambda1*lambda2*lambda2*lambda3/2/delta1/delta1/delta2 
						+ lambda1*lambda1*lambda2*lambda4/2/delta1/delta2/delta2 + lambda2*lambda2*lambda2*lambda4/6/delta1/delta1/delta1;
	
	double rab_amp2b = +(lambda1*lambda1*lambda1*lambda3)/(2*delta2*delta2*delta2) + (lambda1*lambda2*lambda2*lambda3)/(delta1*delta2*delta2) + (lambda1*lambda2*lambda2*lambda3)/(2*delta1*delta1*delta2)
						+(3*lambda1*lambda3*lambda3*lambda3)/(2*delta2*delta2*delta2) + (lambda1*lambda1*lambda2*lambda4)/(2*delta1*delta2*delta2) +(lambda1*lambda1*lambda2*lambda4)/(delta1*delta1*delta2)
						+(lambda2*lambda2*lambda2*lambda4)/(2*delta1*delta1*delta1) + (lambda2*lambda3*lambda3*lambda4)/(2*delta1*delta2*delta2) +(lambda1*lambda3*lambda4*lambda4)/(2*delta1*delta1*delta2)
						+(3*lambda2*lambda4*lambda4*lambda4)/(2*delta1*delta1*delta1) +(lambda1*lambda3*lambda5*lambda5)/(2*delta2*delta3*delta3)+(lambda2*lambda4*lambda5*lambda5)/(2*delta1*delta3*delta3)
						+(lambda2*lambda4*lambda5*lambda5)/(delta1*delta1*delta3) - (lambda2*lambda4*lambda6*lambda6)/(2*delta1*delta1*ddelta)-(lambda1*lambda2*lambda6*lambda7)/(2*delta1*delta2*ddelta)
						-(lambda2*lambda5*lambda6*lambda8)/(2*delta1*delta3*ddelta);
	
	double rab_amp2c = (lambda1*lambda1*lambda1*lambda3)/delta2/delta2/delta2+(0.5*lambda1*lambda2*lambda2*lambda3)/(delta1*delta2*delta2)+(0.5* lambda1*lambda2*lambda2*lambda3)/(delta1*delta1*delta2)
					+(lambda1*lambda3*lambda3*lambda3)/delta2/delta2/delta2+(0.5*lambda1*lambda1*lambda2*lambda4)/(delta1*delta2*delta2)+(lambda2*lambda2*lambda2*lambda4)/delta1/delta1/delta1
					+(0.5*lambda2*lambda3*lambda3*lambda4)/(delta1*delta2*delta2)+(0.5*lambda2*lambda3*lambda3*lambda4)/(delta1*delta1*delta2)+(0.5*lambda1*lambda3*lambda4*lambda4)/(delta1*delta2*delta2)
					+(0.5*lambda1*lambda3*lambda4*lambda4)/(delta1*delta1*delta2)+(lambda2*lambda4*lambda4*lambda4)/delta1/delta1/delta1+(0.5*lambda1*lambda3*lambda5*lambda5)/(delta2*delta3*delta3)
					+(0.5*lambda1*lambda3*lambda5*lambda5)/(delta2*delta2*delta3)+(0.5*lambda2*lambda4*lambda5*lambda5)/(delta1*delta3*delta3)+(0.5*lambda2*lambda4*lambda5*lambda5)/(delta1*delta1*delta3)
					-(0.5*lambda2*lambda4*lambda6*lambda6)/(delta1*delta1*ddelta)-(0.5*lambda1*lambda2*lambda6*lambda7)/(delta1*delta2*ddelta)-(0.5*lambda2*lambda5*lambda6*lambda8)/(delta1*delta3*ddelta);
	
	double shift1 = lambda1*lambda1/delta2 -lambda2*lambda2/delta1 +lambda4*lambda4/delta1 -lambda3*lambda3/delta2 +lambda5*lambda5/delta3;	
																																					
	double shift2 = - 2*lambda1*lambda2*lambda3*lambda4/delta1/delta2/delta2																	
					- lambda1*lambda1*lambda4*lambda4/delta1/delta2/delta2
					- lambda1*lambda1*lambda4*lambda4/delta1/delta1/delta2 
					- lambda4*lambda4*lambda4*lambda4/delta1/delta1/delta1
					- lambda4*lambda4*lambda5*lambda5/delta1/delta3/delta3
					- lambda4*lambda4*lambda5*lambda5/delta1/delta1/delta3
					 
					- lambda1*lambda1*lambda1*lambda1/delta2/delta2/delta2
					- lambda1*lambda1*lambda3*lambda3/delta2/delta2/delta2
					- lambda1*lambda1*lambda5*lambda5/delta2/delta3/delta3
					- lambda1*lambda1*lambda5*lambda5/delta2/delta2/delta3
					
					- lambda5*lambda5*lambda5*lambda5/delta3/delta3/delta3    
					 
					+ lambda2*lambda2*lambda2*lambda2/delta1/delta1/delta1
					+ lambda2*lambda2*lambda3*lambda3/delta1/delta2/delta2 
					+ lambda2*lambda2*lambda3*lambda3/delta1/delta1/delta2
					+ 2*lambda1*lambda2*lambda3*lambda4/delta1/delta1/delta2
					+ 2*lambda2*lambda2*lambda4*lambda4/delta1/delta1/delta1
					+ 0.25*lambda2*lambda2*lambda6*lambda6/delta1/ddelta/ddelta
					- 0.75*lambda2*lambda2*lambda6*lambda6/delta1/delta1/ddelta
			//		- lambda2*lambda3*lambda6*lambda7/delta1/delta2/ddelta
					
					+ lambda3*lambda3*lambda3*lambda3/delta2/delta2/delta2;
	
	double shift2a	= 0.5*lambda2*lambda2*lambda2*lambda2/delta1/delta1/delta1 - 0.5*lambda1*lambda1*lambda2*lambda2/delta1/delta2/delta2 
					- 0.5*lambda1*lambda1*lambda1*lambda1/delta2/delta2/delta2 + 0.5*lambda1*lambda1*lambda2*lambda2/delta1/delta1/delta2;				
	
	double shift2b	= -1;
	shift2b *= (0.5*lambda1*lambda1*lambda1*lambda1)/delta2/delta2/delta2+(0.5*lambda1*lambda1*lambda2*lambda2)/(delta1*delta2*delta2)-(0.5*lambda1*lambda1*lambda2*lambda2)/(delta1*delta1*delta2)
				-(0.5*lambda2*lambda2*lambda2*lambda2)/delta1/delta1/delta1+(2*lambda1*lambda1*lambda3*lambda3)/delta2/delta2/delta2-(0.5*lambda2*lambda2*lambda3*lambda3)/(delta1*delta2*delta2)
				-(0.5*lambda2*lambda2*lambda3*lambda3)/(delta1*delta1*delta2)-(0.5*lambda3*lambda3*lambda3*lambda3)/delta2/delta2/delta2+(2*lambda1*lambda2*lambda3*lambda4)/(delta1*delta2*delta2)
				-(2*lambda1*lambda2*lambda3*lambda4)/(delta1*delta1*delta2)+(0.5*lambda1*lambda1*lambda4*lambda4)/(delta1*delta2*delta2)+(0.5*lambda1*lambda1*lambda4*lambda4)/(delta1*delta1*delta2)
				-(2*lambda2*lambda2*lambda4*lambda4)/delta1/delta1/delta1-(0.5*lambda3*lambda3*lambda4*lambda4)/(delta1*delta2*delta2)+(0.5*lambda3*lambda3*lambda4*lambda4)/(delta1*delta1*delta2)
				+(0.5*lambda4*lambda4*lambda4*lambda4)/delta1/delta1/delta1+(1.5*lambda1*lambda1*lambda5*lambda5)/(delta2*delta3*delta3)+(0.5*lambda1*lambda1*lambda5*lambda5)/(delta2*delta2*delta3)
				-(0.5*lambda2*lambda2*lambda5*lambda5)/(delta1*delta3*delta3)-(0.5*lambda2*lambda2*lambda5*lambda5)/(delta1*delta1*delta3)-(0.5*lambda3*lambda3*lambda5*lambda5)/(delta2*delta3*delta3)
				-(0.5*lambda3*lambda3*lambda5*lambda5)/(delta2*delta2*delta3)+(1.5*lambda4*lambda4*lambda5*lambda5)/(delta1*delta3*delta3)+(0.5*lambda4*lambda4*lambda5*lambda5)/(delta1*delta1*delta3)
				+(1.5*lambda5*lambda5*lambda5*lambda5)/delta3/delta3/delta3-(0.25*lambda2*lambda2*lambda6*lambda6)/(delta1*ddelta*ddelta)+(0.75*lambda2*lambda2*lambda6*lambda6)/(delta1*delta1*ddelta)
				+(lambda2*lambda3*lambda6*lambda7)/(delta1*delta2*ddelta);		
	double rab_amp = M_PI/2/rab_amp_th;//*(tries*2+130.0)/400.0;
			
	double rab = sqrt(abs(tgate*rab_amp));
		//	cout << tries << "\trab\t" <<rab_amp <<endl;
				



	cout << " lambda1 = "  << lambda1 << endl;
	cout << " lambda2 = "  << lambda2 << endl;
	cout << " lambda3 = "  << lambda3 << endl;
	cout << " lambda4 = "  << lambda4 << endl;
//	cout << " lambda5 = "  << lambda5 << endl;
	cout << " delta1 = "  << delta1 << endl;
	cout << " delta2 = "  << delta2 << endl;
//	cout << " delta3 = "  << delta3 << endl;
//	cout << " deriv1 = "  << lambda2/delta1 << " deriv2=  "<< lambda5/delta3 << endl;
	
	//cout << HdriftQ << endl;		
	//cout << HdriftCav << endl;	
	cout << "Hd " << Hdrift << endl;
	
//	dimQ=0;
	U_desired(dimQ+1,0)=1;//std::complex<double>(0,1/sqrt(2));
	U_desired(0,dimQ+1)=1;//std::complex<double>(0,1/sqrt(2));
	U_desired(0,0)=U_desired(dimQ+1,dimQ+1)=0;//1/sqrt(2);
	U_desired= ExpM::EigenMethod(Hdrift, -ii*tgate)*U_desired;
//	U_desired(dimQ+1,dimQ+1)=-1;
	
	//MOs::Null(HcontrolX);
	//HcontrolX(0,3) = HcontrolX(3,0)=0.01;
	
	cout << "HX " << HcontrolX << endl;	
	
	vector<double> ucontrol0i(num_time), ucontrol1i(num_time), ucontrol2i(num_time), ucontrol3i(num_time);
	
	Control u0(HcontrolX, dt, num_time, 10, &Control::ShiftedGaussian, rab_amp, 2, 2), 
			u1(HcontrolY,dt, num_time, 10, &Control::Derivative, -0.8/delta1, 1, 0), 
			u2(HcontrolZ, dt, num_time, 10, &Control::StarkShift, shift1, 1, 0), 
			u3(HcontrolZ2, dt, num_time, 10, &Control::StarkShift, 0, 1, 0);
	ofstream dataout;
	UFs::OpenFile("rab.dat",dataout, 16);
	
	cout << "rab amp = " <<  rab_amp << endl;

	sys.Setucontrol(&u0,0); //0th control
	sys.Setucontrol(&u1,1); //1st control
	sys.Setucontrol(&u2,2); //2nd control
	sys.Setucontrol(&u3,3); //2nd control
	
	sys.SetHdrift(Hdrift);
	sys.SetRhoDesired(U_desired);

	
int k=0;
//for ( k=0; k<50; k++) 
{
	double sum=0, sum2=0;
	
	double rab = sqrt(abs(tgate*rab_amp) );
	u0.ShiftedGaussian( 0, tgate, rab, tgate/2, &u0);
	//u0.SquarePulse( 0, tgate, rab_amp/tgate);
	//u0.TanhPulse( 0, tgate, M_PI, tgate/16);		
	double norm = u0.SquareNormalize(rab_amp);	//cout << "norm " << norm << endl;
	//u1.GaussianDerivative(0, tgate, -rab_amp*norm/delta1, tgate/2);
	for(size_t j = 0; j < num_time; ++j){
		//ucontrol0i[j]=USs::ShiftedGaussian(j*dt, 0, tgate-1, rab_amp, tgate/4)/2.3;
		//ucontrol0i[j]= rab_amp/tgate;
		sum+=ucontrol0i[j]*dt;
		//ucontrol1i[j]=USs::GaussianDerivative(j*dt, 0, tgate-1, rab_amp, tgate/4)/delta1;
		//	sum2+=ucontrol0i[j]*dt;
		//	ucontrol1[j]=USs::GaussianDerivative(j*dt, 0, tgate-1, M_PI, tgate/2)/20;
		//cout << ucontrol0i[j] << "\t" << ucontrol1i[j]<< endl;
		u2.u_[j]=0;
		u3.u_[j]=u0.u_[j]*u0.u_[j]*(lambda1*lambda1/(delta2)-lambda2*lambda2/delta1+lambda4*lambda4/delta1 -lambda3*lambda3/delta2);//+lambda5*lambda5*u0.u_[j]*u0.u_[j]/delta3-0.0*lambda1*lambda5*u0.u_[j]*u0.u_[j]/delta3;
		
						//cav ram
//				theControls_[0]->u_[j]-= sqrt(abs(x*x*x*x*shift2));
//				theControls_[0]->u_[j] = sqrt(abs( x * x - rab_amp2c/rab_amp_th*x*x*x*x));
	}
	//u0.ShiftedGaussian( 0, tgate, rab_amp, tgate/2);
				//cav ram
			//theControls_[1]->Null();
//			theControls_[1]->init(theControls_[0]) ;//Derivative(0, 0, -0.8/delta1, 0, theControls_[0]);			
			//theControls_[2]->Null();
//			theControls_[2]->StarkShift((tries2+90.0)/100.0*shift1, theControls_[0]);
//			theControls_[2]->init(theControls_[0]);//StarkShift(0, 0, shift1, 0, theControls_[0]);
//			theControls_[2]->Stark2Shift(shift2b, theControls_[0]);

//			theControls_[2]->SquarePulse(0,tgate_,(tries2+110.0)/100.0*theControls_[0]->u_[num_time_/2]*theControls_[0]->u_[num_time_/2]*shift1);
//			theControls_[2]->Pixelate();

	cout << sum/rab_amp << endl;
	cout << "c0 " <<  ucontrol0i[2] << " c2  " << ucontrol2i[2]  <<endl;

	//sys.UnitaryTransfer(static_cast< Fid >(&GrapeUnitaryCavity::Phi43levCav), static_cast< GradFid >( &GrapeUnitaryCavity::GradPhi43levCav) ); 
	sys.sweeptimes(argv[1], static_cast< Fid >(&GrapeUnitaryCavity::Phi43levCav), static_cast< GradFid >( &GrapeUnitaryCavity::GradPhi43levCav));
	
	//sys.UnitaryTransfer(static_cast< Fid >(&GrapeUnitaryCavity::Phi4TraceCav), static_cast< GradFid >( &GrapeUnitaryCavity::GradPhi4TraceCav) ); 
	
	//dataout << k << "\t" << sys.Phi43levCav() << "\n" ;
	
	//for (int k=0; k<num_time; k++)
	//{	//dataout << k << "\t" <<  abs(sys.rho_[k](0,0))<< "\t" <<  abs(sys.rho_[k](2,2))<< "\t" <<  abs(sys.rho_[k](3,3))<< "\t" <<  abs(sys.rho_[k](0,0)) <<std::endl;
	//	cout << k << "\t" <<  abs(sys.rho_[k](0,4)) <<std::endl;
	//}	
	
}	
	dataout.close();
	
	return 0;
	
//Pick desired Gate
	//PI pulse
	
//	U_desired(2,2)=std::complex<double>(1,0);
	
	//PI/2 pulse	
	//U_desired(0,0)=1/sqrt(2);
	//U_desired(1,0)=-std::complex<double>(0,1/sqrt(2));
	//U_desired(1,1)=1/sqrt(2);
	//U_desired(0,1)=-std::complex<double>(0,1/sqrt(2));
	//U_desired(2,2)=1;
	
	//2PI pulse
	//MOs::Identity(U_desired);

	//Hdrift = (delta-0.5*Delta)*n+ 0.5*Delta*n*n;	
///	Hdrift(2,2) =  Delta;	
///	Hdrift(3,3) = 3.0*Delta;	
	
	
	
	
	//1st control
	
	
//	sys.Setucontrol(ucontrol2,2);
//	sys.SetHcontrol(0.0*HcontrolZ,2);				
		
		
	
	//run grape
	double to_ms =1000.0/CLOCKS_PER_SEC;
	int clo ;
	clo = clock();	
	sys.Normalizeucontrol(0, 3.14159265358*1);
//	sys.UnitaryTransfer(&GrapeUnitaryCavity::Phi4Sub2, &GrapeUnitaryCavity::GradPhi4Sub2); 
	
	double* newcontrols0 = sys.Getucontrol(0);
	double* newcontrols1 = sys.Getucontrol(1);
//	double* newcontrols2 = sys.Getucontrol(2);
//	cout << "con3 " << newcontrols2[3] << endl;
	
	clo = to_ms*(clock() - clo);
	cout << "time " << clo << endl;


	return 0;
}