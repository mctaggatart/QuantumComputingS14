/*
 *  Control.hpp
 *  
 *
 *  Created by Felix on 30/09/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef Control_h
#define Control_h 
 
#include "MatrixExponential.hpp"
#include "QuantumOperations.hpp"
#include <fstream>

class Control;

using namespace std;
typedef void (Control::* Hcontroltdep)(size_t j, matrix<complex<double> >* Hreturn);
typedef void (Control::* ptrSubGradients)(size_t j, double gradient, double* gradients);
typedef void (Control::* ptrInterpolate)();


class Control {
	public:
	
	size_t npixels_;//  the number of time points	
	size_t nsubpixels_;//  the number of time points
	double dt_;	
	size_t dim_;
	
	Hcontroltdep pSetHcontrol;				//pixel Hamiltonian element
	Hcontroltdep pSetHgradient;				//pixel Hamiltonian element gradient
	
	ptrInterpolate pInterpolate;			//pixel dependence on control
	ptrSubGradients pSubGradients;			//pixel gradient dependence on control
	
	//The physical system, Hamiltonian, and eigenvectors/values
	matrix<complex<double> > Hcontrol_;
	Evolution* evol;
	matrix<std::complex<double> > Z;   //transormation matrices for each control/drift
	matrix<std::complex<double> > ZeH;   //transormation matrices for each control/drift
	double *W;						//eigenvalues for each control/drift
	
	//memory allocated to controls
	double *u_; //control array of size npixels	
	double *u_filt; //filtered control array of size npixels	
	double *u_best;		
	complex<double> *w_; //FT of control array of size npixels		
	size_t buffer; // number of pixels on either side to set to 0

	//physics variables
	double avg;
	Control* detctrl;
	Control* ampctrl;
	double *freqs_;
	double framefreq_, drivefreq_;
	double relphase_;
	double* config_;
	size_t* nconfigs_;
	double phase_;

	//file variables
	char filename[127];
	char filesuffix[20];

	public:
	inline Control(const matrix<complex<double> > &Hcontrol,const double dt,const size_t num_time, const size_t nsubpixels, Evolution* evol, char* filename) 
		: dt_(dt), npixels_(num_time), nsubpixels_(nsubpixels), evol(evol), config_(NULL), nconfigs_(NULL), relphase_(0), phase_(0)
	{
		u_ = new double[npixels_];
		u_filt = u_;
		u_best = new double[npixels_];
		w_ = new complex<double>[npixels_];
		Hcontrol_ = Hcontrol;
		Hcontrol_.SetOutputStyle(Matrix);	
		dim_ = Hcontrol.GetLD();
		Z.initialize(dim_,dim_);
		Z.SetOutputStyle(Matrix);
		ZeH.initialize(dim_,dim_);
		W = new double[dim_];
	//	ExpM::EigenMethod(Hcontrol_,1.0, &(Z), W);
		Null();
		avg=0;
		buffer=0;
		pSetHcontrol = &Control::linearControl;		
		pSetHgradient = &Control::linearGradient;			
		pInterpolate  = &Control::Pixelate;					
		pSubGradients = &Control::setPixelsGradients;
		if(evol!=NULL)	evol->Setucontrol(this);
		if(filename) strcpy(this->filename, filename);
		strcpy(filesuffix,"");
		if(filename) cout << filename;
		cout << " hamiltonian:\n" <<Hcontrol_ << endl;
				
	}
	
	inline Control(const Control& control, const size_t nsubpixels )
		: nsubpixels_(nsubpixels), config_(NULL), nconfigs_(NULL), relphase_(0)
		{
		npixels_ = control.npixels_*nsubpixels/control.nsubpixels_;
		Hcontrol_ = control.Hcontrol_;
		dim_ = control.dim_;
		Z.initialize(dim_,dim_);
		W = new double[dim_];
		ExpM::EigenMethod(Hcontrol_,1.0, &(Z), W);
		relphase_ = control.relphase_;
		framefreq_ = control.framefreq_;
		phase_ = control.phase_;
		evol = control.evol;
		drivefreq_ = control.drivefreq_;
		strcpy(filename,control.filename);
		u_ = new double[npixels_];
		u_filt = u_;
		u_best = new double[npixels_];
		w_ = new complex<double>[npixels_];
		dt_ = control.dt_*control.nsubpixels_/nsubpixels;
		size_t k=0;	
		for(size_t j =0; j < npixels_; j+=nsubpixels_, k+=control.nsubpixels_)
		{	for(size_t l=0; l<nsubpixels_; l++)	
			{	u_[j+l]=control.u_[k];
				w_[j+l]=control.w_[k];		
			}			
		}
		pSetHcontrol = control.pSetHcontrol;	
		pSetHgradient = control.pSetHgradient;
		pInterpolate = control.pInterpolate;
		pSubGradients = control.pSubGradients;
		avg=control.avg;
		strcpy(filesuffix,"");
		
	}
	
	inline ~Control(){
		if(u_filt!=u_)
			delete [] u_filt;	
		if(u_!=NULL)
			delete [] u_;
		if(u_best!=NULL)
			delete [] u_best;
		if(w_!=NULL)
			delete [] w_;	
		delete [] W;	
	}	
	
	inline void settotaltime(double totaltime)
	{
		settotalpixels(totaltime/dt_);	
	}
	
	
	inline void settotalpixels(size_t newnpixels)
	{
		size_t oldnpixels = npixels_;
		npixels_ = newnpixels;
		if(nsubpixels_==oldnpixels) //square pulses
			nsubpixels_=newnpixels;
		if(newnpixels<=oldnpixels) return;
		double* newcontrols = new double[newnpixels];
		u_best = new double[newnpixels];
		for(int j=oldnpixels; j<newnpixels; j++) newcontrols[j] = 0;
		for(int j=0; j<oldnpixels; j++) newcontrols[j] = u_[j];
		if(u_filt!=u_)
			delete [] u_filt;
		delete [] u_, w_, u_best;
		u_ = u_filt = newcontrols;
		w_ = new complex<double>[newnpixels];
		ExpM::EigenMethod(Hcontrol_,1.0, &(Z), W);	
											
	}

		//pre: npixels mod nsubpixels == 0
	inline void setsubpixels(size_t newnsubpixels)
	{
		npixels_=npixels_/nsubpixels_*newnsubpixels;
		nsubpixels_ = newnsubpixels;
		//update redund subpixels 
	}
				
	inline void setasbest()
	{
		for(size_t j =0; j < npixels_; j++)
			u_best[j] = u_[j];
	}			
				
	inline void writefile(char* outfile=NULL)
	{
		ofstream dataout;

		if(outfile==NULL) outfile = filename;
		
		if(outfile!=NULL)
			UFs::OpenFile(outfile,dataout, 16);
		else return;

		for(size_t j =0; j < npixels_; j++){
			dataout <<  dt_*(j+0.5) <<"\t" << u_filt[j] << "\n";
		}
		dataout.close();
	}
	
	inline void writefilewithsuffix(char* suffix)
	{
		char tempfilename[80];
		if(filename!=NULL && suffix!=NULL)
			writefile(strcat(strcat(strcpy(tempfilename,filename), suffix), ".dat"));
	}
	
	inline void writefilewithprefixandsuffix(char* prefix, char* suffix)
	{
		char tempfilename[80];
		if(filename!=NULL && suffix!=NULL)
			writefile(strcat(strcat(strcat(strcat(strcpy(tempfilename,prefix),"/"), filename), suffix), ".dat"));
	}
	
	inline void readfile(char* infile=NULL)
	{
		ifstream datain;
		
		if(infile==NULL) infile = filename;
		
		if(infile!=NULL)
			datain.open(infile, std::ios::in);
		else return;
		
		char c;
		double d;
		for(size_t j =0; j < 20150-buffer*nsubpixels_; j++){
			datain >> c; datain >> c; datain >> c;
			datain >> d;
			datain >> c; 
			datain >> d;			
		}
		for(size_t j =0; j < npixels_; j++){
			datain >> c; datain >> c; datain >> c;
			datain >> u_[j];
			datain >> c; 
			datain >> u_[j];			
		}
		datain.close();
	}
	
	inline double Normalize(double value)
	{
		double sum=0;
		for(size_t j = 0; j < npixels_; ++j){
			sum+=u_[j];
		}
		if(sum)sum=value/sum/dt_;
		else sum=1;
		for(size_t j = 0; j < npixels_; ++j){
			u_[j]*=sum;
		}
		return sum;
	}
	
	inline void rescale(double value)
	{
		for(size_t j = 0; j < npixels_; ++j){
			u_[j]*=value;
		}
	}

	inline double SquareNormalize(double value)
	{
		double sum=0;
		for(size_t j = 0; j < npixels_; ++j){
			sum+=u_[j]*u_[j];
		}
		//cout << value/sum/dt_ << endl;
		if(sum)sum=sqrt(abs(value/sum/dt_));
		else sum=1;
		for(size_t j = 0; j < npixels_; ++j){
			u_[j]*=sum;
		}
		return sum;
	}

	inline double average()
	{
		double sum=0;
		for(size_t j = 0; j < npixels_; ++j){
			sum+=u_[j];
		}
		return avg=sum/npixels_;
	}
	
	inline void Null()
	{
		for(size_t j =0; j < npixels_; j++)
			u_[j]= 0;
	}
	
	inline void SquarePulse(const double amp){
		// height is amp and duration is td .. thus area is amp*td
		//n pi/2 pulse area = n pi/2 -> amp  = n pi/2td
	 	 double ton=dt_*buffer*nsubpixels_;
		 double td=dt_*(npixels_-2*buffer*nsubpixels_);
		 double t=dt_*0.5;
		 for(size_t j =0; j < npixels_; j+=nsubpixels_/2)
		 {
			u_[j]= (t>ton && t < td+ton) ? amp : 0.0;
			t+=dt_*nsubpixels_/2;
		 }
		 
		 //horns
	//	 u_[buffer*nsubpixels_+nsubpixels_/2] = u_[npixels_-buffer*nsubpixels_-nsubpixels_/2] = 1.2* amp; 
		 Pixelate();
	}
	
	inline void StarkShift(const double amp, const double power, const Control* ctrl){  
		for(size_t j=0; j < npixels_; j++)
{			u_[j] += amp*pow(ctrl->u_[j],power);
//			cout << j << " " << u_[j] << " " << ctrl->u_[j] << endl;
}
	}
	
	inline void Product(const double amp, const double power, const Control* ctrl){  
		for(size_t j=0; j < npixels_; j++)
			u_[j] *= amp*pow(ctrl->u_[j],power);
	}
	
	
	inline void Derivative(const double amp, const Control* ctrl){  
		double last, curr;
		last = ctrl->u_[0];
		u_[0] = amp*(ctrl->u_[1]-ctrl->u_[0])/dt_;
		 for(size_t j =1; j < npixels_-1; j++)
		{
			curr = ctrl->u_[j];
			u_[j] = amp*0.5*(ctrl->u_[j+1]-last)/dt_;
			last = curr;
		}
		 u_[npixels_-1] =amp*(ctrl->u_[npixels_-1]-last)/dt_;
	}	
	
	
	inline void Replicate(Control* source)
	{	
		if(npixels_==source->npixels_)
			for(int j=0; j<npixels_; j++)
				u_[j] = source->u_[j];
		else
			for(int j=0; j<npixels_; j++)
			{	u_[j] = source->u_[(size_t)floor(source->npixels_*(double)j/npixels_)];
	//	cout<<  j << " " << u_[j] << " " << (size_t)round(source->npixels_*(double)j/npixels_) << " " << source->u_[(size_t)round(source->npixels_*(double)j/npixels_)] << endl;
				
				}
	//	for(int j=0; j<source->npixels_; j++)		
	//		cout<<  j << " " << source->u_[j] << endl;
						
	}	
			
					
	inline void Sine(const double amp, const double freq)
	{	 
		 for(size_t j =0; j < npixels_; j++)
		 {
			 u_[j]= amp*sin(j*dt_*freq);
		 } 
		// Pixelate();
	}		
					
	inline void RandomControl(const double min, const double max)
	{	 
		 for(size_t j =0; j < npixels_; j+=nsubpixels_)
		 {
			 u_[j]= min + (rand()/(double)RAND_MAX)*(max - min);
			 for(size_t k =1; k < nsubpixels_; k++)
				u_[j+k]=u_[j];
		 } 
	}
		
	inline void DFT()
	{	 std::complex<double> i(0.0,1.0);
		 for(size_t k =0; k < npixels_; k++)
		{	w_[k] = 0;
			for(size_t j =0; j < npixels_; j++)
			{
				w_[k] = w_[k]+ dt_*u_[j]*(cos(100.0*dt_*j * k * 2.0 * M_PI / npixels_) + i * sin(100.0*dt_*j * k * 2.0 * M_PI / npixels_));
			}	
		 } 
		 for(size_t j =0; j < npixels_; j++)
			u_[j] = abs(w_[j]);
	}

	inline void iDFT()
	{	 std::complex<double> i(0.0,1.0);
		std::complex<double> cnum;
		 for(size_t j =0; j < npixels_; j++)
		 {	cnum = 0;
			for(size_t k =0; k < npixels_; k++)
			{
				cnum = cnum + (w_[k]*(cos(j * k * 2.0 * M_PI / npixels_) - i * sin(j * k * 2.0 * M_PI / npixels_)));
			}	
			u_[j] = real(cnum);
		 } 
	}
	
	
	////////////////////////////////////////////////////////
	//Operations for matrix element / sub-pixel  time-dependence and gradients
	////////////////////////////////////////////////////////
	
	inline void getMatrixControl(size_t j, matrix <complex<double> >* Hreturn)
	{
		(this->*pSetHcontrol)(j,Hreturn);		
	}
	
	inline void getMatrixGradient(size_t j, matrix <complex<double> >* Hreturn)
	{
		(this->*pSetHgradient)(j,Hreturn);		
	}
	
	void linearControl(size_t j, matrix <complex<double> >* Hreturn)
	{
		*Hreturn = u_filt[j] * Hcontrol_; 
	}
	
	void linearGradient(size_t j, matrix <complex<double> >* Hreturn)
	{
		*Hreturn = Hcontrol_; 
	}
	
	void oscillatoryControl(size_t j, matrix <complex<double> >* Hreturn)
	{
		double phase;
		if(config_==NULL || nconfigs_==NULL || nconfigs_==0)
			phase = 0;
		else phase= (*config_)*2*M_PI/(*nconfigs_);
			
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				cout << Hcontrol_(p,q)*exp(-complex<double>(0,-1)*freqs_[q]*(j*dt_ + phase)) << endl;
				(*Hreturn)(p,q)= exp(-std::complex<double>(0.0,-1.0*freqs_[q]*(j*dt_ + phase)))*u_[j]*Hcontrol_(p,q)*exp(std::complex<double>(0.0,-1.0*freqs_[p]*(j*dt_ + phase)));
			}
		}
		
	}
	
	void oscillatoryControlGradient(size_t j, matrix <complex<double> >* Hreturn)
	{
		double phase;
		if(config_==NULL || nconfigs_==NULL || nconfigs_==0)
			phase = 0;
		else phase= (*config_)*2*M_PI/(*nconfigs_);
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				(*Hreturn)(p,q)= exp(-std::complex<double>(0.0,-1.0*freqs_[q]*(j*dt_ + phase)))*Hcontrol_(p,q)*exp(std::complex<double>(0.0,-1.0*freqs_[p]*(j*dt_ + phase)));
			}
		}
	}
	
	void drivenControl(size_t j, matrix <complex<double> >* Hreturn)
	{
	/*	double phase;
		if(config_==NULL || nconfigs_==NULL)
			phase = 0;
		else phase= (*config_)*2*M_PI/(*nconfigs_);*/
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				(*Hreturn)(p,q)= exp(-std::complex<double>(0.0,framefreq_*q*(double)(j+0.5)*dt_ + q*phase_))*u_[j]
				*Hcontrol_(p,q)*cos(-relphase_+drivefreq_*(j+0.5)*dt_ + phase_)*exp(std::complex<double>(0.0,framefreq_*p*(double)(j+0.5)*dt_ + p*phase_));
			}
		}
	}

	void drivenControlGradient(size_t j, matrix <complex<double> >* Hreturn)
	{
	/*	double phase;
		if(config_==NULL || nconfigs_==NULL)
			phase = 0;
		else phase= (*config_)*2*M_PI/(*nconfigs_);*/
		
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				(*Hreturn)(p,q)= exp(-std::complex<double>(0.0,framefreq_*q*(double)(j+0.5)*dt_ + q*phase_))
				*Hcontrol_(p,q)*cos(-relphase_+drivefreq_*(j+0.5)*dt_ + phase_)*exp(std::complex<double>(0.0,framefreq_*p*(double)(j+0.5)*dt_ + p*phase_));
			}
		}
	}
	
	//void phaseramp()
				//theControls_[0]->u_[j]=x*cos(integ) - theControls_[1]->u_[j]*sin(integ);
				//theConwtrols_[1]->u_[j]=theControls_[1]->u_[j]*cos(integ) + x*sin(integ);
/*	void fmdrivenControl(size_t j, matrix <complex<double> >* Hreturn)
	{	
		if(detctrl==NULL)
		{	detctrl=this;
			avg=0;
		}	
		double phase;
		if(config_==NULL || nconfigs_==NULL)
			phase = 0;
		else phase= (*config_)*2*M_PI/(*nconfigs_);
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
			//	(*Hreturn)(p,q)= exp(-std::complex<double>(0.0,(framefreq_+detctrl->u_[npixels_-1]/(npixels_))*q*(double)(j+0.5)*dt_))*u_[j]*Hcontrol_(p,q)*cos(-relphase_+(drivefreq_+detctrl->u_[npixels_-1]/(npixels_))*(j+0.5)*dt_)*exp(std::complex<double>(0.0,(framefreq_+detctrl->u_[npixels_-1]/(npixels_))*p*(double)(j+0.5)*dt_));
				(*Hreturn)(p,q)= exp(std::complex<double>(0.0,(framefreq_-detctrl->avg)*p*(j+0.5)*dt_ + p*phase))*
							u_[j]*Hcontrol_(p,q)*cos(-relphase_+(drivefreq_*dt_-detctrl->u_[j]*M_PI/npixels_)*(j+0.5) + phase)*
							exp(-std::complex<double>(0.0,(framefreq_-detctrl->avg)*q*(j+0.5)*dt_ + q*phase))
							+ 0.5*detctrl->avg*q*(q==p);
		//		if(j>0)	
		//			(*Hreturn)(p,q)+= 0.5* (detctrl->u_[j]*(j+0.5)-detctrl->u_[j-1]*(j-0.5))*q*(q==p); 
			}
		}
	}
	
	void fmdrivenControlGradient(size_t j, matrix <complex<double> >* Hreturn)
	{
		double phase;
		if(config_==NULL || nconfigs_==NULL)
			phase = 0;
		else phase= (*config_)*2*M_PI/(*nconfigs_);
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				(*Hreturn)(p,q)= exp(std::complex<double>(0.0,(framefreq_-detctrl->avg)*p*(j+0.5)*dt_ + p*phase))*
							Hcontrol_(p,q)*cos(-relphase_+(drivefreq_*dt_-detctrl->u_[j]*M_PI/npixels_)*(j+0.5) + phase)
							*exp(-std::complex<double>(0.0,(framefreq_-detctrl->avg)*q*(j+0.5)*dt_ + q*phase));
			}
		}
	}
	
	void drivefreqControl(size_t j, matrix <complex<double> >* Hreturn)
	{
		*Hreturn =  0.0*Hcontrol_;
	}

	void drivefreqControlGradient(size_t j, matrix <complex<double> >* Hreturn) //todo 'avg' assumed to be 0 here
	{
		double phase;
		if(config_==NULL || nconfigs_==NULL)
			phase = 0;
		else phase= (*config_)*2*M_PI/(*nconfigs_);
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				(*Hreturn)(p,q)= +(j+0.5)*M_PI/npixels_*exp(std::complex<double>(0.0,(ampctrl->framefreq_-avg)*p*(j+0.5)*dt_ + p*phase))*ampctrl->u_[j]
							*ampctrl->Hcontrol_(p,q)*sin(-ampctrl->relphase_+(ampctrl->drivefreq_*dt_-u_[j]*M_PI/npixels_)*(j+0.5) + phase)
							*exp(-std::complex<double>(0.0,(ampctrl->framefreq_-avg)*q*(j+0.5)*dt_ + q*phase));
			}
		}
	}	
	
	/*/
	void fmdrivenControl(size_t j, matrix <complex<double> >* Hreturn)
	{	
		if(detctrl==NULL)
		{	detctrl=this;
			avg=0;
		}	
		double phase;
		if(config_==NULL || nconfigs_==NULL)
			phase = 0;
		else phase= (*config_)*2*M_PI/(*nconfigs_);
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
			//	(*Hreturn)(p,q)= exp(-std::complex<double>(0.0,(framefreq_+detctrl->u_[npixels_-1]/(npixels_))*q*(double)(j+0.5)*dt_))*u_[j]*Hcontrol_(p,q)*cos(-relphase_+(drivefreq_+detctrl->u_[npixels_-1]/(npixels_))*(j+0.5)*dt_)*exp(std::complex<double>(0.0,(framefreq_+detctrl->u_[npixels_-1]/(npixels_))*p*(double)(j+0.5)*dt_));
				(*Hreturn)(p,q)= exp(std::complex<double>(0.0,(framefreq_-detctrl->u_[j]*M_PI/npixels_-detctrl->avg)*p*(j+0.5)*dt_ + p*phase))*
							u_[j]*Hcontrol_(p,q)*cos(-relphase_+(drivefreq_*dt_)*(j+0.5) + phase)*
							exp(-std::complex<double>(0.0,(framefreq_-detctrl->u_[j]*M_PI/npixels_-detctrl->avg)*q*(j+0.5)*dt_ + q*phase))
							+ 0.5*detctrl->avg*q*(q==p);
		//		if(j>0)	
		//			(*Hreturn)(p,q)+= 0.5* (detctrl->u_[j]*(j+0.5)-detctrl->u_[j-1]*(j-0.5))*q*(q==p); 
			}
		}
	}
	
	void fmdrivenControlGradient(size_t j, matrix <complex<double> >* Hreturn)
	{
		double phase;
		if(config_==NULL || nconfigs_==NULL)
			phase = 0;
		else phase= (*config_)*2*M_PI/(*nconfigs_);
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				(*Hreturn)(p,q)= exp(std::complex<double>(0.0,(framefreq_-detctrl->u_[j]*M_PI/npixels_-detctrl->avg)*p*(j+0.5)*dt_ + p*phase))*
							Hcontrol_(p,q)*cos(-relphase_+(drivefreq_*dt_)*(j+0.5) + phase)
							*exp(-std::complex<double>(0.0,(framefreq_-detctrl->u_[j]*M_PI/npixels_-detctrl->avg)*q*(j+0.5)*dt_ + q*phase));
			}
		}
	}
	
	void drivefreqControl(size_t j, matrix <complex<double> >* Hreturn)
	{
		*Hreturn =  0.0*Hcontrol_;
	}

	void drivefreqControlGradient(size_t j, matrix <complex<double> >* Hreturn) //todo 'avg' assumed to be 0 here
	{
		if(detctrl==NULL)
		{	detctrl=this;
			avg=0;
		}	
		double phase;
		if(config_==NULL || nconfigs_==NULL)
			phase = 0;
		else phase= (*config_)*2*M_PI/(*nconfigs_);
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
			//	(*Hreturn)(p,q)= exp(-std::complex<double>(0.0,(framefreq_+detctrl->u_[npixels_-1]/(npixels_))*q*(double)(j+0.5)*dt_))*u_[j]*Hcontrol_(p,q)*cos(-relphase_+(drivefreq_+detctrl->u_[npixels_-1]/(npixels_))*(j+0.5)*dt_)*exp(std::complex<double>(0.0,(framefreq_+detctrl->u_[npixels_-1]/(npixels_))*p*(double)(j+0.5)*dt_));
	//			(*Hreturn)(p,q)=               exp(std::complex<double>(0.0,(framefreq_-detctrl->u_[j]*M_PI/npixels_-detctrl->avg)*p*(j+0.5)*dt_ + p*phase))*
	///						u_[j]*Hcontrol_(p,q)*cos(-relphase_+(drivefreq_*dt_)*(j+0.5) + phase)*
	///					exp(-std::complex<double>(0.0,(framefreq_-detctrl->u_[j]*M_PI/npixels_-detctrl->avg)*q*(j+0.5)*dt_ + q*phase));
		//		if(j>0)	
		//			(*Hreturn)(p,q)+= 0.5* (detctrl->u_[j]*(j+0.5)-detctrl->u_[j-1]*(j-0.5))*q*(q==p); 
		
		
				(*Hreturn)(p,q)= 0.0*( std::complex<double>(0.0,(-u_[j]*M_PI/npixels_)*p*(j+0.5)*dt_)-std::complex<double>(0.0,(-u_[j]*M_PI/npixels_)*q*(j+0.5)*dt_) ) 
							* exp(std::complex<double>(0.0,(framefreq_-u_[j]*M_PI/npixels_)*p*(j+0.5)*dt_ + p*phase))*
							ampctrl->u_[j]*Hcontrol_(p,q)*cos(-relphase_+(drivefreq_*dt_)*(j+0.5) + phase)*
							exp(-std::complex<double>(0.0,(framefreq_-u_[j]*M_PI/npixels_)*q*(j+0.5)*dt_ + q*phase));
			}
		}
	}	
	

	
	inline void Interpolate()
	{
		(this->*pInterpolate)();
	}
	
	inline void Void(){
		strcpy(filesuffix,"_nointer");
	}
	
	inline void Pixelate()
	{
		strcpy(filesuffix,"_pixels");
		for(size_t j =0; j < npixels_; j++)
		{	u_filt[j]=u_[j]=u_[(j/nsubpixels_)*nsubpixels_+nsubpixels_/2];
		}
	}
	
	inline void LinearInterpolate()
	{
		strcpy(filesuffix,"_lininter");
		double alpha;
		double temp[npixels_];
		for(size_t j = nsubpixels_/2; j < npixels_; j+=nsubpixels_/2) temp[j] = u_[j];
		for(size_t j = 0; j < npixels_; ++j){
			alpha = ((double)((j+nsubpixels_/2)%nsubpixels_))/nsubpixels_;
			double val = 0, val1=0, val2=0;
			if(j>=npixels_-nsubpixels_/2) alpha*=2;
			if(j>=nsubpixels_/2)
				val1 +=   (1-alpha)*temp[((j+nsubpixels_/2)/nsubpixels_ - 1)*nsubpixels_+nsubpixels_/2];
			else alpha=(alpha-0.5)*2;
			if(j<npixels_-nsubpixels_/2)
				val2 +=   (alpha)*temp[((j+nsubpixels_/2)/nsubpixels_)*nsubpixels_+nsubpixels_/2];
			val = val1+val2;
			u_filt[j]=u_[j] = val;
		}
	}
	
	inline void CubicInterpolate()
	{
		if(u_filt==u_)
		{	u_filt=new double[npixels_];			
		}
		Pixelate();	
		strcpy(filesuffix,"_cubinterpol");
		double y0, y1, y2, y3;
		double h, h3;
		for(size_t jc = 0; jc < npixels_; ++jc){
			u_filt[jc]=0;
			for(size_t offset=0; offset<2; offset++)
			{
				size_t j=jc+1-offset+offset*(nsubpixels_%2);
				//if( !((j-nsubpixels_/2)%nsubpixels_)) continue;
				h = nsubpixels_*dt_;
			//	if(j<nsubpixels_/2 || j>=npixels_-nsubpixels_/2) h/=2;
				h3 = pow(h,3);
				double lx = ((j+nsubpixels_/2)%nsubpixels_)*dt_;
				double rx = lx-h;
				double lx3 = pow(lx,3);
				double rx3 = pow(rx,3);
				double lrx2 = lx*pow(rx,2);
							
				y0 = (jc>=3*nsubpixels_/2) ? u_[((j+nsubpixels_/2)/nsubpixels_ - 2)*nsubpixels_+nsubpixels_/2] : 0;//-u_[nsubpixels_/2] ;//
				y1 = (jc>=nsubpixels_/2) ? u_[((j+nsubpixels_/2)/nsubpixels_ - 1)*nsubpixels_+nsubpixels_/2] : 0;// -u_[nsubpixels_/2] ;//0;//-u_[nsubpixels_/2] ;//0;//
				y2 = (jc<npixels_-nsubpixels_/2) ? u_[((j+nsubpixels_/2)/nsubpixels_)*nsubpixels_+nsubpixels_/2] : 0;//- u_[npixels_-nsubpixels_/2] ;//0;//
				y3 = (jc<npixels_-3*nsubpixels_/2) ? u_[((j+nsubpixels_/2)/nsubpixels_ + 1)*nsubpixels_+nsubpixels_/2] : 0;// -u_[npixels_-nsubpixels_/2] ;//0;//
				
				//u_[j] = y2*lx3/h3 + (6*y2-y1-y3)*rx3/2/h3 + (2*y3 - 11*y2 + 4*y1 - y0)*lrx2/2/h3 + (y3 - y1 -6*y2)*rx/2/h;			
				u_filt[jc] += (-lrx2/2/h3)*y0 + (-rx3/2/h3+4*lrx2/2/h3-rx/2/h)*y1 + (lx3/h3+6*rx3/2/h3-11*lrx2/2/h3-6*rx/2/h)*y2 + (-rx3/2/h3+2*lrx2/2/h3+rx/2/h)*y3;
			}
			u_filt[jc]/=2.0;
			//cout << u_[jc] << " " << u_filt[jc] << endl;
		}
	}
/*inline void CubicInterpolate()
	{
		if(u_filt==u_)
		{	u_filt=new double[npixels_];			
		}
		Pixelate();	
		strcpy(filesuffix,"_cubinterpol");
		double y0, y1, y2, y3;
		double h, h3;
		for(size_t j = 0; j < npixels_; ++j){
			//if( !((j-nsubpixels_/2)%nsubpixels_)) continue;
			h = nsubpixels_*dt_;
		//	if(j<nsubpixels_/2 || j>=npixels_-nsubpixels_/2) h/=2;
			h3 = pow(h,3);
			double lx = ((j+nsubpixels_/2)%nsubpixels_)*dt_;
			double rx = lx-h;
			double lx3 = pow(lx,3);
			double rx3 = pow(rx,3);
			double lrx2 = lx*pow(rx,2);
						
			y0 = (j>=3*nsubpixels_/2) ? u_[((j+nsubpixels_/2)/nsubpixels_ - 2)*nsubpixels_+nsubpixels_/2] : 0;//-u_[nsubpixels_/2] ;//
			y1 = (j>=nsubpixels_/2) ? u_[((j+nsubpixels_/2)/nsubpixels_ - 1)*nsubpixels_+nsubpixels_/2] : 0;// -u_[nsubpixels_/2] ;//0;//-u_[nsubpixels_/2] ;//0;//
			y2 = (j<npixels_-nsubpixels_/2) ? u_[((j+nsubpixels_/2)/nsubpixels_)*nsubpixels_+nsubpixels_/2] : 0;//- u_[npixels_-nsubpixels_/2] ;//0;//
			y3 = (j<npixels_-3*nsubpixels_/2) ? u_[((j+nsubpixels_/2)/nsubpixels_ + 1)*nsubpixels_+nsubpixels_/2] : 0;// -u_[npixels_-nsubpixels_/2] ;//0;//
			
			//u_[j] = y2*lx3/h3 + (6*y2-y1-y3)*rx3/2/h3 + (2*y3 - 11*y2 + 4*y1 - y0)*lrx2/2/h3 + (y3 - y1 -6*y2)*rx/2/h;			
			u_filt[j] = (-lrx2/2/h3)*y0 + (-rx3/2/h3+4*lrx2/2/h3-rx/2/h)*y1 + (lx3/h3+6*rx3/2/h3-11*lrx2/2/h3-6*rx/2/h)*y2 + (-rx3/2/h3+2*lrx2/2/h3+rx/2/h)*y3;
			cout << u_[j] << " " << u_filt[j] << endl;
		}
	}*/
	inline void CubicInterpolate2()
	{
		strcpy(filesuffix,"_cubinter2");
		double y0, y1, y2, y3;
		double h, h3, A, B;
			
		for(size_t j = 0; j < npixels_; ++j){
			//if( !((j-nsubpixels_/2)%nsubpixels_)) continue;
			h = nsubpixels_*dt_;
			double lx = ((j+nsubpixels_/2)%nsubpixels_)*dt_;
			double rx = lx-h;
			
			y0 = (j>=3*nsubpixels_/2) ? u_[((j+nsubpixels_/2)/nsubpixels_ - 2)*nsubpixels_+nsubpixels_/2] :0;// -u_[nsubpixels_/2] ;//
			y1 = (j>=nsubpixels_/2) ? u_[((j+nsubpixels_/2)/nsubpixels_ - 1)*nsubpixels_+nsubpixels_/2] :0;// -u_[nsubpixels_/2] ;//
			y2 = (j<npixels_-nsubpixels_/2) ? u_[((j+nsubpixels_/2)/nsubpixels_)*nsubpixels_+nsubpixels_/2] :0;// - u_[npixels_-nsubpixels_/2] ;
			y3 = (j<npixels_-3*nsubpixels_/2) ? u_[((j+nsubpixels_/2)/nsubpixels_ + 1)*nsubpixels_+nsubpixels_/2] :0;// -u_[npixels_-nsubpixels_/2] ;

			A = -rx/h;
			B = lx/h;

			u_[j] = (1.0/6)*(A*A*A - A)*h*h*(y2+y0-2*y1)/h/h  + (1.0/6)*(B*B*B - B)*h*h*(y3+y1 - 2*y2)/h/h + A*y1 + B*y2;
			//u_[j] = -y0*lrx2/2/h3 + y1*( -rx/2/h + 3 *rx3/2/h3 - 4 * lrx2/h3) + y2*( lx3/h3 + 3*rx3/h3 + lrx2/h3 - 3*rx/h) + y3*(rx/2/h -rx3/2/h3 + lrx2/h3);
			//u_[j] = y2*lx3/h3 + (6*y2-y1-y3)*rx3/2/h3 + (2*y3 - 11*y2 + 4*y1 - y0)*lrx2/2/h3 + (y3 - y1 -6*y2)*rx/2/h;
		}
	}

	inline void GaussianFilter()
	{
		if(u_filt==u_)
		{	u_filt=new double[npixels_];			
		}	
		Pixelate();	
		strcpy(filesuffix,"_gaussfilt");
		double area=0;
		int left=-(2*nsubpixels_), right=2*nsubpixels_;
		double amps[right-left];
		for(int k = left; k < right; ++k)
			area+=amps[k-left]=exp(-(double)k*k*dt_*dt_/0.5) ;///0.03125);	//0,529 experimental value	
//		double* res = new double[npixels_];

		for(int j = 0; j < npixels_; ++j){
			u_filt[j]=0.0;
			for(int k = left; k < right; ++k)		
				if(j+k>=0)
					if(j+k<npixels_)	u_filt[j] +=  amps[k-left]*u_[j+k];
//					else u_filt[j] +=  amps[k-left]*u_[npixels_-1];
//				else  u_filt[j] +=  amps[k-left]*u_[0];
			if(area) u_filt[j]/=area;	
		}
//		delete u_;
//		u_ = res;
	}

	inline void setSubGradients(size_t j, double gradient, double* gradients)
  {
   
		(this->*pSubGradients)(j, gradient, gradients);		
	}

	inline void setNoSubPixelGradients(size_t j, double gradient, double* gradients){	
	  gradients[j] = gradient;				
	}
	
	inline void setAMGradient(size_t j, double gradient, double* gradients){
		gradients[0] += gradient;
		if(!j) {
			gradients[0]/=npixels_/u_[npixels_/2];
			for(int l = npixels_-1; l >=0 ; --l)
				gradients[l] = gradients[0]*u_[l];
			
		}	 				
	}

	inline void setPixelsGradients(size_t j, double gradient, double* gradients)
  {
    //std::cout<<"GRADIENT TEST, j is "<<j<<"subpixels min equ "<<buffer*nsubpixels_<<"j is less than "<< npixels_-buffer*nsubpixels_;
    if(j>=buffer*nsubpixels_ && j<npixels_-buffer*nsubpixels_){
      //  std::cout<<"j is"<< j;
      gradients[(j/nsubpixels_)*nsubpixels_+nsubpixels_/2] += gradient;}
	}

	inline void setLinearSplineGradients(size_t j, double gradient, double* gradients)
	{
		double alpha = ((double)((j+nsubpixels_/2) % nsubpixels_)) /nsubpixels_;
		if(j>=npixels_-nsubpixels_/2) alpha*=2;
		if(j>=nsubpixels_/2+buffer*nsubpixels_)
				gradients[((j+nsubpixels_/2)/nsubpixels_ - 1)*nsubpixels_+nsubpixels_/2] += (1-alpha)*gradient;	
		else alpha=(alpha-0.5)*2;
		if(j<npixels_-nsubpixels_/2-buffer*nsubpixels_)
				gradients[((j+nsubpixels_/2)/nsubpixels_)*nsubpixels_+nsubpixels_/2] += alpha*gradient;		
	}

	inline void GaussianFilterGradient(size_t j, double gradient, double* gradients)
	{
		double area=0;
		int left=-(2*nsubpixels_), right=2*nsubpixels_;
		double amps[right-left];
		for(int k = left; k < right; ++k)
			area+=amps[k-left]=exp(-(double)k*k*dt_*dt_/0.529);		
				
		for(int k = left; k < right; ++k)				
			if(j+k>=buffer*nsubpixels_ && j+k<npixels_-buffer*nsubpixels_)	
				gradients[((j+k)/nsubpixels_)*nsubpixels_+nsubpixels_/2] +=  amps[k-left]*gradient;				
				
	}	
	
	inline void setCubicSplineGradients(size_t j, double gradient, double* gradients)
	{
		size_t x0, x1, x2, x3;
		double h, h3;
		
		for(size_t offset=0; offset<2; offset++)
		{
			size_t jc=j+1-offset+offset*(nsubpixels_%2);
			//if( !((jc-nsubpixels_/2)%nsubpixels_)) continue;
			h = nsubpixels_*dt_;
			//	if(jc<nsubpixels_/2 || j>=npixels_-nsubpixels_/2) h/=2;
			h3 = pow(h,3);
			double lx = ((jc+nsubpixels_/2)%nsubpixels_)*dt_;
			double rx = lx-h;
			double lx3 = pow(lx,3);
			double rx3 = pow(rx,3);
			double lrx2 = lx*pow(rx,2);
			
			if(jc>=3*nsubpixels_/2) 
				gradients[((jc+nsubpixels_/2)/nsubpixels_ - 2)*nsubpixels_+nsubpixels_/2]+= -lrx2/2/h3*gradient;
			if(jc>=nsubpixels_/2) 
				gradients[((jc+nsubpixels_/2)/nsubpixels_ - 1)*nsubpixels_+nsubpixels_/2]+= (-rx3/2/h3+4*lrx2/2/h3-rx/2/h)*gradient;
			if(jc<npixels_-nsubpixels_/2) 
				gradients[((jc+nsubpixels_/2)/nsubpixels_)*nsubpixels_+nsubpixels_/2]+= (lx3/h3+6*rx3/2/h3-11*lrx2/2/h3-6*rx/2/h)*gradient;
			if(jc<npixels_-3*nsubpixels_/2) 
				gradients[((jc+nsubpixels_/2)/nsubpixels_ + 1)*nsubpixels_+nsubpixels_/2]+= (-rx3/2/h3+2*lrx2/2/h3+rx/2/h)*gradient;
				
		}
	}
	
	
	
	
};


#endif
