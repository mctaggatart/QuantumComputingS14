
/*
Name: Evolution.h
Author: Jay M. Gambetta and Felix Motzoi, Anastasia McTaggart
This does the actual computation of the gradient			
Dependences: Controls, but relies on
include "MatrixExponential.hpp"
#include "MatrixOperations.hpp"
#include "MatrixOperations.hpp"
#include "QuantumOperations.hpp"	
#include "QuantumOperations.hpp"
*/

#ifndef Evolution_h
#define Evolution_h
	
class Control;
class Evolution;
	
typedef void (Evolution::* ptrPropagate)();
#include "MatrixExponential.hpp"
#include "MatrixOperations.hpp"
#include "MatrixOperations.hpp"
#include "QuantumOperations.hpp"	
#include "QuantumOperations.hpp"
class Evolution {
	public:
		Evolution(size_t dim, size_t num_controls, size_t num_time, const double h);// number of controls, the number of time points
		virtual ~Evolution();
	
		void SetNumTimes(const size_t newnumtimes);
		
		matrix<std::complex<double> >* GetHdrift();
		void SetHdrift( matrix<std::complex<double> >& H_drift);
		
		matrix<std::complex<double> > GetRhoInitial();
  void SetRhoInitial(const matrix<std::complex<double> >& rho_initial); //initial state
		
		Control* GetuControl(const size_t k);
		void Setucontrol(Control* control, const int index=-1);
		
  matrix<std::complex<double> > GetUDesired() { return U_desired_; } //desired Unitiary
  void SetUDesired(const matrix<std::complex<double> >& U_desired);
  void SetTrueRhoDesired(const matrix<std::complex<double> >& Udes);//sets desired Rho

		void writepopulations(char* outfile);
						
		ptrPropagate pPropagate;


		void forwardpropagate();
                void forwardpropagate_Density();
                void forwardpropagate_Taylor();
		void forwardpropagate_Commute();
		
                matrix<std::complex<double> >* U_; //an array of the forward evolutions (rho in paper
                matrix<std::complex<double> >* Rho1_;
                matrix<std::complex<double> >** UC_; //an array of the forward evolutions (rho in paper) splitting by control
  double* freqs_;
						
	public:
		size_t num_time_;//  the number of time points	
		double h_; //the time step
		double tgate_;
		size_t dim_;//dimensions of the Hilbert space
		double one_on_dim_;//used alot so I have defined it
		
		matrix<std::complex<double> > *H_drift_; //drift hamiltonian
		matrix<std::complex<double> > H_drift_exp; //exponential of drift hamiltonian
		matrix<std::complex<double> > U_desired_; //the desired rho or unitary			//rename as U_desired_
  matrix<std::complex<double> > true_rho_desired_;
		matrix<std::complex<double> > rho_initial_; //the initial rho or unitary
		
		size_t num_controls_; //number of controls,
		Control** theControls_; //list of pointers to controls		
		size_t *controlsetflag_;//a array which when all 1 all controls and drift are set
		
		matrix<std::complex<double> >* Unitaries_;//an array of the unitaries to be applied at each time step
		matrix<std::complex<double> >** UnitariesC_;//an array of the unitaries to be applied at each time step for each control
		std::complex<double> i; // the complex number i
		size_t lastcontrol;  //index of the last control to be set
		
		long times0, times1;			//timing the fidelity and gradient steps of the algorithm

		
		matrix<std::complex<double> > *Z;   //transormation matrices for each time step
		double **W;							//eigenvalues for each time step
};

#include "Control.hpp"

inline Evolution::Evolution(size_t dim, size_t num_controls, size_t num_time, const double h) 
	: dim_(dim), num_controls_(num_controls), num_time_(num_time), h_(h), lastcontrol(0)
{	
  freqs_=NULL;
	Z = new matrix<std::complex<double> >[num_time_];	
	W = new double *[num_time_];
	for(size_t j = 0; j < num_time_; ++j)
	{ Z[j].initialize(dim_,dim_);
	  W[j] = new double[dim_];	
	}
	//SET ppropogate to the one you want to use!
	pPropagate = &Evolution::forwardpropagate_Density;
	theControls_ = new Control*[num_controls_];
	controlsetflag_=new size_t[num_controls_+3]; // one more for Hdrift and one more for rho_desired and rho_initial
	UC_ = new  matrix<std::complex<double> >*[num_controls_];
	UnitariesC_ = new  matrix<std::complex<double> >*[num_controls_];
	for(size_t k=0; k < num_controls_; ++k){
	      
		UC_[k] = new matrix<std::complex<double> >[num_time_];
		UnitariesC_[k] = new matrix<std::complex<double> >[num_time_];
		for(size_t j=0; j < num_time_; ++j){
			UC_[k][j].initialize(dim_,dim_);
			UnitariesC_[k][j].initialize(dim_,dim_);
		}
		controlsetflag_[k]=1;//set for the Hcontros
	}
	//config_=0;
	controlsetflag_[num_controls_]=1; //set for the Hdrift
	controlsetflag_[num_controls_+1]=1; //set for the rho_desired
	controlsetflag_[num_controls_+2]=1; //set for the rho_initial
	lastcontrol = 0;
	U_=new  matrix<std::complex<double> >[num_time_];
	Rho1_= new matrix<std::complex<double> >[num_time_];
	Unitaries_=new  matrix<std::complex<double> >[num_time_];
	for(size_t j=0; j < num_time_; ++j){
		U_[j].initialize(dim_,dim_);
		Rho1_[j].initialize(dim_, dim_);
		Unitaries_[j].initialize(dim_,dim_);
	}
	
	i=std::complex<double>(0.0,1.0);
	one_on_dim_=1/double(dim_);	
	H_drift_exp.initialize(dim_,dim_);
	rho_initial_.initialize(dim_,dim_);
	rho_initial_(0,0)=1;
}

Evolution::~Evolution(){
	for( size_t k = 0; k < num_controls_; ++k)
	{	delete [] UC_[k];
		delete [] UnitariesC_[k];
	}
	for(size_t j=0; j < num_time_; ++j){
		
		delete [] W[j];
	}
	delete [] UC_;
	delete [] UnitariesC_;
	delete [] theControls_;
	delete [] controlsetflag_;
	delete [] U_;
	delete [] Rho1_;
	delete [] Unitaries_;
	delete []Z;
	delete []W;
}

inline void Evolution::SetNumTimes(const size_t newnumtimes)
{
	for(size_t j=0; j < num_time_; ++j)
		delete [] W[j];
	
	num_time_=newnumtimes;
	tgate_=h_*num_time_;
	cout << h_ << " " << num_time_ << " " << tgate_ << endl;
	delete [] U_;
	delete [] Rho1_;
	delete [] Unitaries_;
	delete []Z;
	delete []W;
	
	U_=new matrix<std::complex<double> >[num_time_];
	Rho1_=new matrix<std::complex<double> >[num_time_];
	Unitaries_=new  matrix<std::complex<double> >[num_time_];
	for(size_t j=0; j < num_time_; ++j){
		U_[j].initialize(dim_,dim_);
		Rho1_[j].initialize(dim_, dim_);
		Unitaries_[j].initialize(dim_,dim_);
	}

	for(size_t k=0; k < num_controls_; ++k){
		theControls_[k]->settotalpixels(num_time_);
		delete [] UC_[k];
		delete [] UnitariesC_[k];
		UC_[k] = new matrix<std::complex<double> >[num_time_];
		UnitariesC_[k] = new matrix<std::complex<double> >[num_time_];
		for(size_t j=0; j < num_time_; ++j)
		{	UC_[k][j].initialize(dim_,dim_);
			UnitariesC_[k][j].initialize(dim_,dim_);
		}
	}
	Z = new matrix<std::complex<double> >[num_time_];	
	W = new double *[num_time_];
	for(size_t j = 0; j < num_time_; ++j)
	{  Z[j].initialize(dim_,dim_);
	   W[j] = new double[dim_];	
	}
	H_drift_exp = ExpM::EigenMethod(*H_drift_,-i*(h_/num_controls_));
	for(int k=0; k<num_controls_; k++)
		if(!controlsetflag_[k])
			theControls_[k]->ZeH= MOs::Dagger(theControls_[k]->Z)*H_drift_exp;
		
}

//Set the drift Hamiltonian
inline void Evolution::SetHdrift( matrix<std::complex<double> >& H_drift){
	H_drift_=&H_drift;
	H_drift_->SetOutputStyle(Matrix);
	H_drift_exp.SetOutputStyle(Matrix);
	controlsetflag_[num_controls_]=0;
	
	H_drift_exp = ExpM::EigenMethod(H_drift,-i*(h_/num_controls_));
	for(int k=0; k<num_controls_; k++)
		if(!controlsetflag_[k])
			theControls_[k]->ZeH= MOs::Dagger(theControls_[k]->Z)*H_drift_exp;
	if(verbose==yes)
	{
		std::cout << "--------------------Hdrift is set---------------------------" << std::endl;
		std::cout << *H_drift_ << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}
//sets controls for unitaries
inline void Evolution::Setucontrol(Control* control, const int k){
	if(k>=0) lastcontrol=k;
	
	theControls_[lastcontrol]=control;
	
//	if(!controlsetflag_[num_controls_])
//		theControls_[lastcontrol]->ZeH= MOs::Dagger(theControls_[lastcontrol]->Z)*H_drift_exp;
	
	controlsetflag_[lastcontrol++]=0;
}

		
//Set the target state for transfer
inline void Evolution::SetUDesired(const matrix<std::complex<double> >& U_desired){
	U_desired_=U_desired;
	controlsetflag_[num_controls_+1]=0;
	if(verbose==yes)
	{
		U_desired_.SetOutputStyle(Matrix);
	 	std::cout << "--------------------rho desired is set---------------------------" << std::endl;
		std::cout << U_desired_ << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}
//sets actual rho desired
inline void Evolution::SetTrueRhoDesired(const matrix<std::complex<double> >& Udes){
  true_rho_desired_=Udes*rho_initial_*MOs::Dagger(Udes);
	controlsetflag_[num_controls_+1]=0;
	if(verbose==yes)
	{
	  true_rho_desired_.SetOutputStyle(Matrix);
	 	std::cout << "--------------------True rho desired is set---------------------------" << std::endl;
		std::cout << true_rho_desired_ << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}


//Set the initial state for transfer
inline void Evolution::SetRhoInitial(const matrix<std::complex<double> >& rho_initial){
	rho_initial_=rho_initial;
	controlsetflag_[num_controls_+2]=0;
	if(verbose==yes)
	{
		rho_initial_.SetOutputStyle(Matrix);
	 	std::cout << "--------------------rho initial is set---------------------------" << std::endl;
		std::cout << rho_initial_ << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}



//Return system values
inline matrix<std::complex<double> >* Evolution::GetHdrift(){
	return H_drift_;
}

inline matrix<std::complex<double> > Evolution::GetRhoInitial(){
	return rho_initial_;
}

inline Control* Evolution::GetuControl(const size_t k){
	return theControls_[k];
}


//writes populations	
inline void Evolution::writepopulations(char* outfile)
{	
	ofstream popout;
	UFs::OpenFile(outfile,popout, 16);
	
	matrix<std::complex<double> > pop = rho_initial_;			
							
	for(size_t j =0; j < num_time_; j++){
		popout << h_* j;
		//modified to use Rho_[j] instead of U_[j], as rho is the populations
		//pop = U_[j]*rho_initial_*MOs::Dagger(U_[j]);
		pop= Rho1_[j];
		
		for(int d=0; d<dim_; d++)
			popout << '\t' << pop(d,d);
		popout << endl;
	}
	popout.close();
}

//forward propogate for unitaries
   void Evolution::forwardpropagate()
   {
   	matrix<std::complex<double> > ident(dim_,dim_);  MOs::Identity(ident);
   	matrix<std::complex<double> > Htemp_(dim_,dim_), Htemp2_(dim_,dim_);	
 
   	matrix<std::complex<double> >* curUnitary = &ident;
	std::complex<double>* freqs=new std::complex<double>[dim_];
	//  freqs_ = new std::complex<double>[dim_];       
   	curUnitary = &ident;
   	for(size_t j = 0; j < num_time_; ++j){
   		Htemp_= *H_drift_;
			
//any rotating frame, commented out to avoid seg faults			
			if(freqs_!=NULL)
 		{	for(size_t d = 0; d < dim_; ++d)
 				freqs[d] = exp(std::complex<double>(0.0,-freqs_[d]*(j*h_))); 
 			MOs::diagmult(Htemp_, freqs, Htemp_);
 			for(size_t d = 0; d < dim_; ++d)
 				freqs[d] = exp(std::complex<double>(0.0,freqs_[d]*(j*h_)));
 			MOs::multdiag(Htemp_, Htemp_, freqs);
 		}
		
// 		//total hamiltonian
   		for(size_t k = 0; k < num_controls_; ++k){	
   			theControls_[k]->getMatrixControl(j,  &Htemp2_);
      
   			Htemp2_.SetOutputStyle(Matrix);
   			Htemp_ = Htemp_ + Htemp2_;
			
			
   		}

   		Htemp_.SetOutputStyle(Matrix);
   	       
   		//timestep propagator Unitaries_[j] and total propagator U_[j]
   		Unitaries_[j]= ExpM::EigenMethod(Htemp_,-i*h_, &(Z[j]), W[j]);	
   		U_[j] = Unitaries_[j]*(*curUnitary);
	    
   		curUnitary = &(U_[j]);
   	}	
}	
//forward propogation for density matrices
 void Evolution::forwardpropagate_Density()
   {
   	matrix<std::complex<double> > ident(dim_,dim_);  MOs::Identity(ident);
   	matrix<std::complex<double> > Htemp_(dim_,dim_), Htemp2_(dim_,dim_);	
	matrix<std::complex<double> >* lastRho = &rho_initial_;
   	matrix<std::complex<double> >* curUnitary = &ident;
 	//std::complex<double>* freqs  = new std::complex<double>[dim_];
   
	lastRho = &rho_initial_;
   	curUnitary = &ident;
   	for(size_t j = 0; j < num_time_; ++j){
   		Htemp_= *H_drift_;
			
	//any rotating frame, commented out to avoid seg faults			
	/*	if(freqs_!=NULL)
 		{	for(size_t d = 0; d < dim_-1; ++d)
 				freqs[d] = exp(std::complex<double>(0.0,-freqs_[d]*(j*h_))); 
 			MOs::diagmult(Htemp_, freqs, Htemp_);
 			for(size_t d = 0; d < dim_-1; ++d)
 				freqs[d] = exp(std::complex<double>(0.0,freqs_[d]*(j*h_)));
 			MOs::multdiag(Htemp_, Htemp_, freqs);
 		}
	*/


		
// 		//total hamiltonian
   		for(size_t k = 0; k < num_controls_; ++k){	
   			theControls_[k]->getMatrixControl(j,  &Htemp2_);
   	
   			Htemp2_.SetOutputStyle(Matrix);
   			Htemp_ = Htemp_ + Htemp2_;
			
			
   		}

   		Htemp_.SetOutputStyle(Matrix);
   		
   		//timestep propagator Unitaries_[j] and total propagator U_[j]
   		Unitaries_[j]= ExpM::EigenMethod(Htemp_,-i*h_, &(Z[j]), W[j]);	
   		Rho1_[j] = Unitaries_[j]*(*lastRho)*(MOs::Dagger(Unitaries_[j]));
		U_[j]=Unitaries_[j]*(*curUnitary);
	       
   		curUnitary = &(U_[j]);
		lastRho=&(Rho1_[j]);
   	}	
}	
//forward propagate for unitaries with taylor methods
void Evolution::forwardpropagate_Taylor()
{	matrix<std::complex<double> > Htemp_(dim_,dim_), Htemp2_(dim_,dim_);
  	static matrix<std::complex<double> > ident(dim_,dim_); 
	MOs::Identity(ident);
	matrix<std::complex<double> >* curUnitary = &ident;
 	
	curUnitary = &ident;
	
	
	for(size_t j = 0; j < num_time_; ++j){			
		U_[j] = (*curUnitary);
		
		Htemp_ = *H_drift_;
		for(size_t k = 0; k < num_controls_; ++k){	
			theControls_[k]->getMatrixControl(j,  &Htemp2_);
			for(size_t q=0; q< dim_; ++q){
				for(size_t p=0; p<dim_; ++p){
					Htemp_(q,p) += Htemp2_(q,p);
				}
			}
		}
		
		int scaling = 3;
		Htemp_=(1.0/pow(2.0, scaling)) *Htemp_; 
		MOs::Identity(Unitaries_[j]);	
		Unitaries_[j] = Unitaries_[j] + (-i*h_)*Htemp_ + 0.5*(-i*h_)*(-i*h_)*Htemp_*Htemp_ + 0.5*0.33333333333*(-i*h_)*(-i*h_)*(-i*h_)*Htemp_*Htemp_*Htemp_;
	
		for(size_t s=0; s<scaling; s++) //squaring
			Unitaries_[j] = Unitaries_[j]*Unitaries_[j];
		
		U_[j] = Unitaries_[j]*(*curUnitary);
	       
		curUnitary=&(U_[j]);
		}
  	
}		
//forward propagate commute for unitaries
void Evolution::forwardpropagate_Commute()
{

  matrix<std::complex<double> > tempmat(dim_,dim_), ident(dim_,dim_), *propag;
	std::complex<double> scalarexponents[dim_];
	std::complex<double> freqs[dim_];
	MOs::Identity(ident);
	propag=&ident;
	int start, inc;
	for(size_t j = 0; j < num_time_; j++){			
	
		if(j%8==0 || j%8==3 || j%8==5 || j%8==6) {		
			start=num_controls_-1; inc=-1;
		} else {	
			start=0; inc=1;
		}
			
		for(int k = start; k < num_controls_ && k>=0;  k+=inc){				
			for (size_t p=0; p < dim_; p++) 	
				scalarexponents[p] = std::exp(-i*h_*theControls_[k]->u_[j]*theControls_[k]->W[p]);
			MOs::multdiag(tempmat,  theControls_[k]->Z, scalarexponents);
			
			UnitariesC_[k][j] = tempmat*theControls_[k]->ZeH;
			
				//rotating frame
				for(size_t d = 0; d < dim_; ++d)
					freqs[d] = exp(std::complex<double>(0.0,-freqs_[d]*(j*h_)));
				MOs::diagmult(UnitariesC_[k][j], freqs, UnitariesC_[k][j]);
				for(size_t d = 0; d < dim_; ++d)
					freqs[d] = exp(std::complex<double>(0.0,freqs_[d]*(j*h_)));
				MOs::multdiag(UnitariesC_[k][j], UnitariesC_[k][j], freqs);
							
			UC_[k][j] = UnitariesC_[k][j] * *propag;
			propag = &(UC_[k][j]);
		}
		U_[j] = *propag;
		
					}

	  }
	



#endif /*Evolution_h */

