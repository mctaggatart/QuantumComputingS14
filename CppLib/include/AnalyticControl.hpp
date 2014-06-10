/*
 *  Control.hpp
 *  
 *
 *  Created by Felix on 30/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */
#include "MatrixExponential.hpp"
#include "QuantumOperations.hpp"
#include <fstream>
#include "Control.hpp"
#include "OptimizeEvolution.hpp"

class AnalyticControl;
typedef void (AnalyticControl::* ptrSetTransition)(matrix<complex<double> >* H_drift, matrix<complex<double> >* H_drive);

class AnalyticControl: public Control {
	protected:
	
	complex<double> *sysdata[4];
	ptrSetTransition pSetTransition;
	void (AnalyticControl::* pInit_)(const double amp, const double sigma, complex<double>** sysparams);
	AnalyticControl* refcontrol;	
	
	public:
	size_t bNorm_;
	double amp_, alpha_;
	double scaling;

	inline AnalyticControl(const matrix<complex<double> > &Hcontrol,const double dt,const size_t num_time, const size_t nsubpixels, Evolution* evol, char* filename=NULL) 
		: Control(Hcontrol,dt,num_time, nsubpixels, evol, filename), scaling(0) 
	{	
		pSetTransition = &AnalyticControl::setdriving01;
		pInit_= &AnalyticControl::Null;
		refcontrol=NULL;
		if(evol) (this->*pSetTransition)(evol->GetHdrift(), &(evol->GetuControl(0)->Hcontrol_));	
	}
	
	inline AnalyticControl(const AnalyticControl& control, const size_t nsubpixels )
		: Control(control, nsubpixels), amp_(control.amp_), bNorm_(control.bNorm_), alpha_(control.alpha_), scaling(0)
	{	pSetTransition = control.pSetTransition;
		pInit_ = control.pInit_;
		refcontrol=control.refcontrol;
		//	(this->*pSetTransition)(evol.GetHdrift());
		if(sysdata!=NULL)
			for(size_t j=0;j<4;j++)
				sysdata[j] = control.sysdata[j];
	}
	
	inline ~AnalyticControl(){	}
	
	inline void setparams(void (AnalyticControl::* pInit)(const double amp, const double alpha, complex<double>** sysparams), 
			double amp, double alpha, size_t bNorm, ptrSetTransition pSetTransition=NULL , AnalyticControl* refcontrol=NULL)
	{
		 this->pInit_=pInit;
		 this->amp_=amp;
		 this->bNorm_=bNorm;
		 this->alpha_=alpha;
		 this->refcontrol = refcontrol;
		 if(pSetTransition!=NULL) this->pSetTransition = pSetTransition;
		 
		 if(evol!=NULL) (this->*this->pSetTransition)(evol->GetHdrift(),  &(evol->GetuControl(0)->Hcontrol_));	
	}
	
	inline void init()
	{

		Null(); 
		if(refcontrol!=NULL) //&&npixels_!=refcontrol->npixels_)
		{
			refcontrol->settotaltime(npixels_*dt_); 
			refcontrol->init();
		//	((OptimizeEvolution*)(refcontrol->evol))->UnitaryTransfer();
		}		
		

		if(bNorm_==2) { 
		  double rtAmpSq = sqrt(abs(dt_*npixels_*amp_));
		  (this->*pInit_)(rtAmpSq, alpha_, sysdata);
		  
		} else { 
		  (this->*pInit_)(amp_, alpha_, sysdata);
		}
				
		writefilewithsuffix("_0");			
		
		Interpolate();
		
		if(bNorm_==2) SquareNormalize(amp_);
		else if(bNorm_==1) Normalize(amp_);
		
		writefilewithsuffix(filesuffix);
			
	}	 
	
	inline void setdriving01(matrix<complex<double> >* H_drift, matrix<complex<double> >* H_drive)
	{
		if(H_drift==NULL || H_drift->GetLD()<3 || H_drive->GetLD()<3){ return; }
		static complex<double>empty=0;
		sysdata[0] = &((*H_drift)(2,2));
		sysdata[1] =  &((*H_drive)(1,2));
		sysdata[2] =  &((*H_drift)(0,0));
		sysdata[3] =  &empty;	
	}
	
	inline void setdriving12(matrix<complex<double> >* H_drift, matrix<complex<double> >* H_drive)
	{
		if(H_drift==NULL) return;
		sysdata[0] = &((*H_drift)(4,4));
		sysdata[1] =  &((*H_drive)(2,3));
		sysdata[2] =  &((*H_drift)(0,0));
		sysdata[3] =  &((*H_drive)(0,1));	
	}
	
	using Control::Null;	
	inline void Null(const double amp, const double sigma, complex<double>** sys_params)
	{	Null();
	}
	
	using Control::Replicate;
	inline void Replicate(const double amp, const double sigma, complex<double>** sys_params)
	{	
		for(int j=0; j<npixels_; j++)
			u_[j] = refcontrol->u_[(size_t)round(refcontrol->npixels_*(double)j/npixels_)];
	}
	
	using Control::RandomControl;
	inline void RandomControl(const double min, const double max, complex<double>** sys_params)
	{	RandomControl(min, max);
	}	
	
	using Control::Sine;
	inline void Sine(const double amp, const double freq, complex<double>** sys_params)
	{	Sine(amp, freq);
	}
	
	
	inline void scanamp(const double min, const double max, complex<double>** sys_params)
	{
		for(size_t j =0; j < npixels_; j++)
		 {
			u_[j]= (min + scaling * (max-min)) * refcontrol->u_[j];
		 }	
	}
	
	inline void RandomMatrixControl(const double min, const double max, complex<double>** sys_params)
	{	
		MOs::RandomHermitian(*(evol->GetHdrift()));
		*(evol->GetHdrift()) = *(evol->GetHdrift()) * 5.0;
		MOs::RandomHermitian(Hcontrol_);
		RandomControl(min, max);
	}		
	
	using Control::SquarePulse;
	inline void SquarePulse(const double amp, const double sigma, complex<double>** sys_params)
	{	SquarePulse(amp);
	}
	
	inline void Gaussian(const double amp, const double alpha, complex<double>** sys_params){
		//normalized so that the area is amp
		 double sigma = dt_*npixels_/alpha;
		 double td = dt_*npixels_;
		 double t=dt_*0.5;
		 for(size_t j =0; j < npixels_; j++)
		 {
			u_[j]= amp*exp(-0.5*pow((t-td*0.5),2)/pow(sigma,2)) / (sigma*M_SQRT2*sqrt(M_PI));
			t+=dt_;
		 }
	}	
	inline void ShiftedGaussian(const double amp, const double alpha, complex<double>** sys_params){
		//normalized so that the area is ~amp
		 double sigma = dt_*npixels_/alpha;
		 double td = dt_*npixels_;
		 double  t=dt_*0.5;
		 for(size_t j =0; j < npixels_; j++)
		 {
			u_[j]= M_PI*amp*(exp(-0.5*pow((t-td*0.5)/sigma,2)) - exp(-0.5*pow((-td*0.5)/sigma,2)) ) / (sigma*M_SQRT2*sqrt(M_PI))/erf(td*0.5*M_SQRT1_2/sigma);
			t+=dt_;
		 }
	}

	
	
	using Control::StarkShift;
	inline void StarkShift(const double amp, const double power, complex<double>** sys_params)  
	{	StarkShift(amp, power, refcontrol);
	}
	
	using Control::Product;
	inline void Product(const double amp, const double power, complex<double>** sys_params)  
	{	Product(amp, power, refcontrol);
	}
	
	inline void DerivativeDRAGB2(const double amp, const double sigma, complex<double>** sys_params){  
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		double Delta2 = real(*sys_params[2]);
		double lamb2 = 2*real(*sys_params[3]);
		if(!Delta || ! Delta2) 
			UFs::MyError("Divide by zero in DerivativeDRAGB2\n");		
		Derivative(-(lambda*lambda/Delta-lamb2*lamb2/Delta2)/4, refcontrol);
	}

	inline void Derivative7DRAGC(const double amp, const double sigma, complex<double>** sys_params){  
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		
		StarkShift( -1/Delta, 1, refcontrol);
		StarkShift( -11*(-2+lambda*lambda)/2/pow(Delta,3)/3/pow(2.0,2), 3, refcontrol);
		StarkShift( -(2924-2438*pow(lambda,2)+587*pow(lambda,4) )/24/pow(Delta,5)/5/pow(2.0,4), 5, refcontrol);
		Derivative(1, this);
		
		Control temp(*refcontrol, refcontrol->nsubpixels_);
		temp.Null();
		temp.StarkShift(1,1,refcontrol);
		temp.Derivative(1,&temp);
		StarkShift(-(-47 +22*pow(lambda,2))/(1*pow(Delta,5))/pow(2.0,2),3,&temp);
		
		Control temp2(*refcontrol, refcontrol->nsubpixels_);
		temp2.Null();
		temp2.StarkShift(1,1,&temp);
		temp2.Derivative(1,&temp2);
		
		temp.Product(1,1,&temp2);
		temp.Product( -2*(-71 +33*pow(lambda,2))/(1*pow(Delta,5))/pow(2.0,2),1,refcontrol);
		StarkShift(1,1,&temp);
		
		temp2.Derivative(1,&temp2);
		temp2.Product( -(-24 +11*pow(lambda,2))/(1*pow(Delta,5))/pow(2.0,2),2,refcontrol);
		StarkShift(1,1,&temp2);		
		
	}
	
	using Control::Derivative;
	inline void Derivative(const double amp, const double sigma, complex<double>** sys_params){  
	  std::cout<<"Derivative";
		Derivative(amp, refcontrol);
	}	
	
	inline void DerivativeTrig(const double amp, const double sigma, complex<double>** sys_params){  
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		
		for(size_t j =0; j < npixels_; j++)
		{
			double A = refcontrol->u_[j]/Delta/2;
			double sum;
			for (int m=1; m<6;m++)
			{
				sum = refcontrol->u_[j]/Delta/2;
				double fac = 1;
				for (int n=1; n<m+1;n++)
				{
					fac *= (2*n-1)*(2*n);
					sum += (pow(-1.0,n)/fac)*pow(sqrt(3)*A,2*n)*( refcontrol->u_[j]*pow(3.0,2*n) +2*Delta*A*(1-pow(2.0,2*n+2))/(2*n+1) )/Delta/6;
				}
				A=sum;				
			}
			u_[j]=-2*sum;			
		}		
		Derivative(1,this);
	}	
	
	
	inline void StarkShiftDRAGAFreq(const double amp, const double sigma, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		Null();
		StarkShift( ((lambda*lambda)/Delta)/4/4, 2, refcontrol);
		average();
		for(int j=1; j<npixels_; j++)
			u_[j] += u_[j-1];	
		for(int j=0; j<npixels_; j++)
			u_[j]/=(j+0.5);
	//		for(int j=0; j<npixels_; j++)
	//		u_[j]-=2*avg;
		for(int j=0; j<npixels_; j++)
			u_[j]=avg;
	}
	
	inline void StarkShiftDRAGA(const double amp=0, const double sigma=0, complex<double>** sys_params=NULL){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		cout << Delta << " " << lambda << endl;
		Null();
		StarkShift( ((lambda*lambda)/Delta)/4, 2, refcontrol);
		average();
	//	for(int j=0; j<npixels_; j++)
	//		u_[j] = avg;
	}

	inline void StarkShiftDRAGA2(const double amp, const double sigma, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		double Delta2 = real(*sys_params[2]);
		double lamb2 = 2*real(*sys_params[3]);
		if(!Delta || ! Delta2) 
			UFs::MyError("Divide by zero in StarkShiftDRAGA2\n");
		StarkShift( ( -lamb2*lamb2/Delta2 +(lambda*lambda)/Delta)/4, 2, refcontrol);
	}

	inline void Twophoton(const double amp, const double sigma, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		StarkShift( -lambda/8/Delta/Delta, 2, refcontrol);
		Derivative(1,this);
	}
	
	inline void TwophotonTrig(const double amp, const double sigma, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		
		Control temp(*refcontrol, refcontrol->nsubpixels_);
		temp.Null();
		temp.StarkShift(1,1,refcontrol);
		temp.Derivative(1,&temp);
		for(size_t j =0; j < npixels_; j++)
		 {
//			double A = ctrl->u_[j]/Delta/2 - pow(ctrl->u_[j], 3) / (4*pow(Delta,3));// +((231* pow(ctrl->u_[j]/Delta, 5))/(640)) -((4437* pow(ctrl->u_[j]/Delta, 7))/(7168))+((39959* pow(ctrl->u_[j]/Delta, 9))/(32768));
			double A = refcontrol->u_[j]/Delta/2;
			double sum;
			for (int m=1; m<6;m++)
			{
				sum = refcontrol->u_[j]/Delta/2;
				double fac = 1;
				for (int n=1; n<m+1;n++)
				{
					fac *= (2*n-1)*(2*n);
					sum += (pow(-1.0,n)/fac)*pow(sqrt(3)*A,2*n)*( refcontrol->u_[j]*pow(3.0,2*n) +2*Delta*A*(1-pow(2.0,2*n+2))/(2*n+1) )/Delta/6;
				}
				A=sum;
			}		
			u_[j]= -((4*sqrt(2)*pow(sin((sqrt(3)*A)/2),2)*(6*Delta + 2*Delta*cos(sqrt(3)*A) + Delta*cos(2*sqrt(3)*A) + 3*sqrt(3)*refcontrol->u_[j]*sin(sqrt(3)*A)))
			/(9*(8 - 3*cos(sqrt(3)*A) + cos(3*sqrt(3)*A))));
		}
	}
	
	
	inline void StarkShiftDRAGC(const double amp, const double sigma, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		StarkShift( (lambda*lambda - 4)/Delta/4, 2, refcontrol);
		average();
		for(int j=0; j<npixels_; j++)
			u_[j] = avg;
	}
	
	inline void StarkShift4DRAGC(const double amp, const double sigma, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		StarkShift( (lambda*lambda - 4)/Delta/4, 2, refcontrol);
		StarkShift( -(12 - 7*lambda*lambda + lambda*lambda*lambda*lambda)/Delta/Delta/Delta/pow(2.0,6), 4, refcontrol);
		StarkShift( (-206 + 170*lambda*lambda - 47*pow(lambda,4) +3*pow(lambda,6))/3/pow(Delta,5)/pow(2.0,6), 6, refcontrol);
		
		Control temp(*refcontrol, refcontrol->nsubpixels_);
		temp.Null();
		temp.StarkShift(1,1,refcontrol);
		temp.Derivative(1,&temp);
		temp.Derivative(1,&temp);
		temp.Product( 12*(24 -11*pow(lambda,2))/(3*pow(Delta,5))/pow(2.0,3),3,refcontrol);
		//Null();
		StarkShift(1,1,&temp);
		
	}
	
	inline void StarkShiftSelect(const double amp, const double sigma,  complex<double>** sys_params){
		Null();
		StarkShift( 1, 1, refcontrol);
		StarkShift( amp, 3, refcontrol);		
	}
	
	inline void StarkShiftTrig(const double amp, const double sigma, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		
		Control temp(*refcontrol, refcontrol->nsubpixels_);
		temp.Null();
		temp.StarkShift(1,1,refcontrol);
		temp.Derivative(1,&temp);
		for(size_t j =0; j < npixels_; j++)
		 {
//			double A = refcontrol->u_[j]/Delta/2 - pow(refcontrol->u_[j], 3) / (4*pow(Delta,3));// +((231*pow(refcontrol->u_[j]/Delta, 5))/(640)) -((4437* pow(refcontrol->u_[j]/Delta, 7))/(7168))+((39959* pow(refcontrol->u_[j]/Delta, 9))/(32768));
			double A = refcontrol->u_[j]/Delta/2;
			double sum;
			for (int m=1; m<6;m++)
			{
				sum = refcontrol->u_[j]/Delta/2;
				double fac = 1;
				for (int n=1; n<m+1;n++)
				{
					fac *= (2*n-1)*(2*n);
					sum += (pow(-1.0,n)/fac)*pow(sqrt(3)*A,2*n)*( refcontrol->u_[j]*pow(3.0,2*n) +2*Delta*A*(1-pow(2.0,2*n+2))/(2*n+1) )/Delta/6;
				}
				A=sum;
			}
			
			u_[j]= (4.0/9)*(Delta + (3*(-6*Delta + 4*Delta*cos(sqrt(3)*A) + sqrt(3)*refcontrol->u_[j]*(-2*sin(sqrt(3)*A) + sin(2*sqrt(3)*A))))/(8 - 3*cos(sqrt(3)*A) + cos(3*sqrt(3)*A)));
			//u_[j]= (4/9)*(Delta + (6*(-3*Delta + 2*Delta*cos(sqrt(3)*A) + sqrt(3)*ctrl->u_[j]/2*(-2*sin(sqrt(3)*A) + sin(2*sqrt(3)*A))))/(8 - 3 *cos(sqrt(3)*A) + cos(3*sqrt(3)*A)));
		//if(j==39) cout << j << " " << sum << endl;					
			
		}
	}
	
	inline void StarkShiftDRAGC2(const double amp, const double sigma, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		double Delta2 = real(*sys_params[2]);
		double lamb2 = 2*real(*sys_params[3]);
		if(!Delta || ! Delta2) 
			UFs::MyError("Divide by zero in StarkShiftDRAGC2\n");
		StarkShift( ( -lamb2*lamb2/Delta2 +(lambda*lambda - 4)/Delta)/4, 2, refcontrol);
	}
	
	
	inline void TanhPulse(const double amp, const double alpha, complex<double>** sys_params){
		// height is amp and duration is td area is approximately amp*td
		//n pi/2 pulse area = n pi/2 -> amp  = n pi/2td
			
		 double sigma = dt_*npixels_/alpha;
		 double td=dt_*npixels_;
		 double t=dt_*0.5;
		 double one_on_sigma=1.0/sigma;
		 for(size_t j =0; j < npixels_; j++)
		 {
			u_[j]=  0.5*amp*(tanh(one_on_sigma*(t)) + tanh(-one_on_sigma*(t-td)) );
			u_[j]-= 0.5*amp*(tanh(-one_on_sigma*(-td)) );
			t+=dt_;
		}		
	}

	inline void GaussianDRAGA(const double amp, const double alpha, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		ShiftedGaussian(amp, alpha, sys_params);
		Normalize(amp);		
		StarkShift( (lambda*lambda/Delta/Delta)/8, 3, this);
	}
	inline void GaussianDRAGA2(const double amp, const double alpha, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		double Delta2 = real(*sys_params[2]);
		double lamb2 = 2*real(*sys_params[3]);
		if(!Delta || ! Delta2) 
			UFs::MyError("Divide by zero in GaussianDRAGA2\n");		
		ShiftedGaussian(amp, alpha, sys_params);
		Normalize(amp);		
		StarkShift( (lambda*lambda/Delta/Delta+lamb2*lamb2/Delta2/Delta2)/8, 3, this);
	}
	inline void GaussianDRAGB(const double amp, const double alpha, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		if(!Delta) 
			UFs::MyError("Divide by zero in GaussianDRAGB\n");		
		ShiftedGaussian(amp, alpha, sys_params);
		Normalize(amp);		
		StarkShift( -0, 3, this);
	}
	
	inline void GaussianDRAGB2(const double amp, const double alpha, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		double Delta2 = real(*sys_params[2]);
		double lamb2 = 2*real(*sys_params[3]);
		if(!Delta || ! Delta2) 
			UFs::MyError("Divide by zero in GaussianDRAGB2\n");		
		ShiftedGaussian(amp, alpha, sys_params);
		Normalize(amp);		
		StarkShift( ( lamb2*lamb2/Delta2/Delta2 )/4, 3, this);
	}
	inline void GaussianDRAGC2(const double amp, const double alpha, complex<double>** sys_params){		
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		double Delta2 = real(*sys_params[2]);
		double lamb2 = 2*real(*sys_params[3]);
		if(!Delta || ! Delta2) 
			UFs::MyError("Divide by zero in GaussianDRAGC2\n");
		ShiftedGaussian(amp, alpha, sys_params);
		Normalize(amp);		
		StarkShift( -( (4-lambda*lambda)/Delta/Delta - lamb2*lamb2/Delta2/Delta2 )/8, 3, this);
	}
	inline void GaussianDRAGC(const double amp, const double alpha, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		ShiftedGaussian(amp, alpha, sys_params);
		Normalize(amp);		
		StarkShift( ( (lambda*lambda-4)/Delta/Delta )/8, 3, this);
	}
	inline void Gaussian5DRAGC(const double amp, const double alpha, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		ShiftedGaussian(amp, alpha, sys_params);
		Normalize(amp);		
		StarkShift( ( (lambda*lambda-4)/Delta/Delta )/8, 3, this);
		StarkShift( -(112  - 76*lambda*lambda + 13*pow(lambda,4))/(8*pow(Delta,4))/pow(2.0,4), 5, refcontrol);
		StarkShift( (-5504  + 5768*lambda*lambda - 2120*pow(lambda,4) +255*pow(lambda,6))/(48*pow(Delta,6))/pow(2.0,6), 7, refcontrol);		
		
		Control temp(*refcontrol, refcontrol->nsubpixels_);
		temp.Null();
		temp.StarkShift(1,1,refcontrol);
		temp.Derivative(1,&temp);
		temp.Derivative(1,&temp);
		temp.Product( (576 -286*pow(lambda,2) + 11*pow(lambda,4) )/(6*pow(Delta,6))/pow(2.0,4),4,refcontrol);
		//Null();
		StarkShift(1,1,&temp);
		
	}
	
	inline void GaussianTrig(const double amp, const double alpha, complex<double>** sys_params){
		double Delta = real(*sys_params[0]);
		double lambda = 2*real(*sys_params[1]);
		
		Control temp(*refcontrol, refcontrol->nsubpixels_);
		temp.Null();
		temp.StarkShift(1,1,refcontrol);
		temp.Derivative(1,&temp);
		for(size_t j =0; j < npixels_; j++)
		 {
//			double A = refcontrol->u_[j]/Delta/2 - pow(ctrl->u_[j], 3) / (4*pow(Delta,3));// +((231* pow(ctrl->u_[j]/Delta, 5))/(640)) -((4437* pow(ctrl->u_[j]/Delta, 7))/(7168))+((39959* pow(ctrl->u_[j]/Delta, 9))/(32768));
			double A = refcontrol->u_[j]/Delta/2;
			double sum;
			for (int m=1; m<6;m++)
			{
		//		if(m==6) cout << j << "j\n";
				sum = refcontrol->u_[j]/Delta/2;
				double fac = 1;
				for (int n=1; n<m+1;n++)
				{
					fac *= (2*n-1)*(2*n);
					sum += (pow(-1.0,n)/fac)*pow(sqrt(3)*A,2*n)*( refcontrol->u_[j]*pow(3.0,2*n) +2*Delta*A*(1-pow(2.0,2*n+2))/(2*n+1) )/Delta/6;
				}
				A=sum;
			}			
			u_[j]= (2*(9*refcontrol->u_[j]*cos(sqrt(3)*A) - 4*sqrt(3)*Delta*pow(sin(sqrt(3)*A),3)))/(3*(8 - 3*cos(sqrt(3)*A) + cos(3*sqrt(3)*A)));
			//if(j==39) cout << j << " " << sum << endl;					
			//u_[j]= 2*(2* (9 * ctrl->u_[j]/2 * cos(sqrt(3)* A) - 2 * sqrt(3) * Delta * pow(sin(sqrt(3) *A),3) ))/(3* (8 - 3*cos(sqrt(3)* A) + cos(3* sqrt(3)* A)));
		 }
	}
	
	inline void GaussianDerivative(const double amp, const double alpha, complex<double>** sys_params){
		 double sigma = dt_*npixels_/alpha;
		double td = dt_*npixels_;
		 double t=dt_*0.5;
		 for(size_t j =0; j < npixels_; j++)
		 {
			u_[j]= -2*0.5*(t-td*0.5)/pow(sigma,2) * M_PI*amp*exp(-0.5*pow((t-td*0.5)/sigma,2)) / (sigma*M_SQRT2*sqrt(M_PI))/erf(td*0.5*M_SQRT1_2/sigma);
			t+=dt_;
		 }
	}
	
	inline void TruncatedGaussian(const double amp, const double alpha, complex<double>** sys_params){
		double sigma = dt_*npixels_/alpha;
		double td = dt_*npixels_;
		 double t=dt_*0.5;
		 for(size_t j =0; j < npixels_; j++)
		 {
			u_[j]= amp*exp(-0.5*pow((t-td*0.5),2)/pow(sigma,2))/(sigma*M_SQRT2*sqrt(M_PI))/erf(td*0.5*M_SQRT1_2/sigma);
			t+=dt_;
		 } 	
	}
	
	

	
	
};
