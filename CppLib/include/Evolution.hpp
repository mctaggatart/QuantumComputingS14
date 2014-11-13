/*
Author: Jay M. Gambetta and Felix Motzoi, Anastasia McTaggart
This does the actual comutation of the gradient			
Dependences: Controls, but relies on
 "MatrixExponential.hpp" 
 "MatrixOperations.hpp"
 "MatrixOperations.hpp"
 "QuantumOperations.hpp"	
 "QuantumOperations.hpp
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
  Evolution(size_t dim, size_t num_controls, size_t num_time, const double h, matrix<std::complex<double> >** dis,size_t numDis, size_t typeDis);// number of controls, the number of time points
		virtual ~Evolution();

				//getter functions
 matrix<std::complex<double> >*	GetHdrift();
 matrix<std::complex<double> >  GetUDesired() { return U_desired_; }	//desired Unitiary
 Control*		        GetuControl(const size_t k);
 void		                SetHdrift( matrix<std::complex<double> >& H_drift);
  
 matrix<std::complex<double> > GetRhoInitial();
  //setter functions
  void SetNumTimes(const size_t newnumtimes);
  
  void SetOmega(double omega);//used for dissipation, T1, T2 depend on this

  void SetRhoInitial(const matrix<std::complex<double> >& rho_initial); //initial state
  //sets the desired unitary.		
  void SetUDesired(const matrix<std::complex<double> >& U_desired);
  void SetTrueRhoDesired(const matrix<std::complex<double> >& Udes);//sets desired Rho	final
  void Setucontrol(Control* control, const int index=-1); //sets the controls

  void SetOppsDesired(matrix<std::complex<double> >** Udes1); //sets our desired opperator		
 
  void writepopulations(char* outfile); //writes populations to a graph
						
  ptrPropagate pPropagate; //function pointer to propagate function of choice

  //different forms of forward propagate                
  void		forwardpropagate(); //standard for unitaries
  void		forwardpropagate_Density();//for density matrices
                void		forwardpropagate_Taylor();
		void	     	forwardpropagate_Commute();
  matrix<std::complex<double> > dissipator1(double dt, std::complex<double> gamma, matrix<std::complex<double> > A, matrix<std::complex<double> >* rho);//first form of the dissipator 
  matrix<std::complex<double> > dissipator2(double dt, std::complex<double> gamma, matrix<std::complex<double> > A, matrix<std::complex<double> >* rho); //2nd form, used in the code
  matrix<std::complex<double> > GenDis(int j,double dt, std::complex<double> gamma,matrix<std::complex<double> >* rho); //actually does our dissipator
  matrix<std::complex<double> > GenKrausDis(int j,double dt, matrix<std::complex<double> >* rho); //for kraus operators
  //matrix<std::complex<double> > SetDis(matrix<std::complex<double> > dis[],int typeDis_, int numDis_);
                matrix<std::complex<double> >*	U_; //an array of the forward evolutions (rho in paper
              
  matrix<std::complex<double> >*	Rho1_; //our rho at each time step
                matrix<std::complex<double> >** UC_; //an array of the forward evolutions (rho in paper) splitting by control
  double*					freqs_; //used in prior versions
						
	public:
  matrix<std::complex<double> >**	A_; //dissipator array
		size_t	num_time_;	//  the number of time points	
		double	h_;	//the time step
		double	tgate_;
		size_t	dim_;	//dimensions of the Hilbert space
		double	one_on_dim_;	//used alot so I have defined it
                std::complex<double>	T1A,T1B,T1C, T2A, T2B, T2C;	//T1 and T2 decoherance operators

               matrix<std::complex<double> > *H_drift_;	//drift hamiltonian
	       matrix<std::complex<double> > H_drift_exp;	//exponential of drift hamiltonia

	       matrix<std::complex<double> > U_desired_;	//the desired unitary, previously called Rho_Desired.
               matrix<std::complex<double> > true_rho_desired_;	//the desired rho
		matrix<std::complex<double> > rho_initial_;	//the initial rho or unitary
  size_t numDis_;
  size_t typeDis_;
		size_t				 num_controls_;	//number of controls,
		Control**			 theControls_;	//list of pointers to controls		
		size_t				*controlsetflag_;	//a array which when all 1 all controls and drift are set
  double					 omega_;
              
		matrix<std::complex<double> >*	 Unitaries_;	//an array of the unitaries to be applied at each time step
		matrix<std::complex<double> >**	 UnitariesC_;	//an array of the unitaries to be applied at each time step for each control
		std::complex<double>		 i;	// the complex number i
		size_t				 lastcontrol;	//index of the last control to be set
		
		long times0, times1;			//timing the fidelity and gradient steps of the algorithm

		
		matrix<std::complex<double> > *Z;	//transormation matrices for each time step
		double	**W;	//eigenvalues for each time step
};

#include "Control.hpp"
//Evolution takes Dimensions, number controls, number of timesteps, dt, the dissipator in matrix form, the number of dissipators in matrix form, and their number of types.
inline Evolution::Evolution(size_t dim, size_t num_controls, size_t num_time, const double h, matrix<std::complex<double> >** dis, size_t numDis, size_t typeDis) 
  : dim_(dim), num_controls_(num_controls), num_time_(num_time), h_(h), lastcontrol(0), numDis_(numDis), typeDis_(typeDis)
{
  //this method is primarily initializing everything.
  //following lines used in old forms of the code
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
	//A_ is our list of dissipators, and initalizes it
	A_=new matrix<std::complex<double> >*[num_time_];
	//initializes our list of unitaries
	Unitaries_=new  matrix<std::complex<double> >[num_time_];

       	for(int i=0; i<num_time_; ++i){
	  //sets to an appropriate sized array for each time step
	A_[i]= new matrix<std::complex<double> >[numDis_*typeDis_];
	//initialize our list of dissipators
	for(size_t k=0; k < numDis_*typeDis_; ++k){
	  	  A_[i][k].initialize(dim_,dim_);

	}
	}
	for(size_t j=0; j < num_time_; ++j){
	  // initialize other matrices to numtime with proper dimensions
		U_[j].initialize(dim_,dim_);
		Rho1_[j].initialize(dim_, dim_);
		Unitaries_[j].initialize(dim_,dim_);
	}
	
	i=std::complex<double>(0.0,1.0);
	one_on_dim_=1/double(dim_);	
	H_drift_exp.initialize(dim_,dim_);
		rho_initial_.initialize(dim_,dim_);
		//arbitrarily chosen initial value. Should be reset by code already, but check.
		rho_initial_(1,1)=1;
	       

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
	delete [] A_;
	delete [] Unitaries_;
	delete []Z;
	delete []W;
}
inline void Evolution::SetOmega(double omega){
  omega_=omega;
  //omega is an arbitrarily chosen value, should eventually be stochastic/take into account background noise
	//relaxation constant, should be omega/10
  T1A=omega_/10.0;
	  T1B=omega_/10.0;
	T1C=omega_/10.0;
	//phase decoherance
	T2A=omega_/10.0;
	//set to less than T1 b/c physically less phase decoherance than relaxation should be omega/20
	T2B=omega_/10.0;
	T2C=omega_/10.0;
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
	delete [] A_;
	delete [] Unitaries_;
	delete []Z;
	delete []W;
	
	U_=new matrix<std::complex<double> >[num_time_];
	Rho1_=new matrix<std::complex<double> >[num_time_];
	A_= new matrix<std::complex<double> >*[num_time_];
	for(size_t j=0; j<num_time_; ++j){
A_[j]= new matrix<std::complex<double> >[numDis_*typeDis_];
	for(size_t k=0; k < numDis_*typeDis_; ++k){
	   A_[j][k].initialize(dim_,dim_);
	}
	}
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
  matrix<std::complex<double> > rho_des1;
  //set Rho desired to U*RHO*U^Daggar
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


//sets our dissipator value
inline void Evolution::SetOppsDesired(matrix<std::complex<double> >** Udes1){

  for(int k=0; k<num_time_; ++k){
  for(size_t i=0; i<( numDis_)*(typeDis_); ++i){
    //set our actual value of our dissipator at each step
  A_[k][i]=Udes1[k][i];
  }
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
							
	for(size_t j =0; j < num_time_; ++j){
		popout << h_* j;
		//modified to use Rho_[j] instead of U_[j], as rho is the populations
		//pop = U_[j]*rho_initial_*MOs::Dagger(U_[j]); old form
		pop= Rho1_[j];
		
		for(int d=0; d<dim_; d++)
			popout << '\t' << pop(d,d);
		popout << endl;
	}
	popout.close();
}

inline matrix<std::complex<double > > Evolution::dissipator2(double dt,std::complex<double> gamma, matrix<std::complex<double> > A, matrix<std::complex<double> >* rho){
 
  std::complex<double> dt_;
  dt_=dt;
  //alternate form of dissipator. dt*D*Rho*D^Dagger -.5(D^Dagger*D*Rho-.5*Rho*D^Dagger*D
 return (dt)*(A*(*rho)*MOs::Dagger(A)-0.5*(MOs::Dagger(A)*A*(*rho))-0.5*((*rho)*MOs::Dagger(A)*A));


}
/*Method one for making a dissipator*/
inline matrix<std::complex<double > > Evolution::dissipator1(double dt,std::complex<double> gamma, matrix<std::complex<double> > A, matrix<std::complex<double> >* rho){
  std::complex<double> dt_; //turn dt into a complex double for type casting reasons
  dt_=dt;

	matrix<std::complex<double> > ident(dim_,dim_);  MOs::Identity(ident);
	//return I+Sqrt(dT*gamma)*D-.5*D^Dagger*D*Dt*gamma
	return (ident+std::sqrt(dt_*gamma)*A-0.5*MOs::Dagger(A)*A*dt_*gamma);
	

}
//Krauss Dissipators
matrix<std::complex<double> > Evolution::GenKrausDis(int j,double dt, matrix<std::complex<double> >* rho){

  matrix<std::complex<double> > opperator=A_[j][0];
  matrix<std::complex<double> > tempdis= A_[j][0]*(*rho)*MOs::Dagger(A_[j][0]);
  
  for(size_t k=1; k<(typeDis_)*(numDis_); k++){
    //the form of the dissipator  
    tempdis = A_[j][k]*(*rho)*MOs::Dagger(A_[j][k]); //NEED TO MAKE A A GLOBAL VARIABLE SET!!
    }
  
 
   return tempdis;
}
//creates our dissipator
matrix<std::complex<double> > Evolution::GenDis(int j,double dt,std::complex<double> gamma, matrix<std::complex<double> >* rho){

  matrix<std::complex<double> > tempdis=dissipator2(dt, gamma, A_[j][0], rho);
 
  for(size_t k=1; k<(typeDis_)*(numDis_); k++){
    // cout<<k<<"k"<<endl;
    tempdis += dissipator2(dt, gamma, A_[j][k], rho); //NEED TO MAKE A A GLOBAL VARIABLE SET!!
    }
  
  return tempdis;
}

//forward propogate for unitaries
   void Evolution::forwardpropagate()
   {
   	matrix<std::complex<double> > ident(dim_,dim_);  MOs::Identity(ident);
   	matrix<std::complex<double> > Htemp_(dim_,dim_), Htemp2_(dim_,dim_);	
 
   	matrix<std::complex<double> >* curUnitary = &ident;
	std::complex<double>* freqs=new std::complex<double>[dim_];
	//	  freqs_ = new std::complex<double>[dim_];       
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
     matrix<std::complex<double> > PauliZ; //PauliZ done by hand
     matrix<std::complex<double> > Anhil(dim_,dim_); MOs::Destroy(Anhil);//produce dim x dim Sigma Minus matrix (eg, spin lowering)
     matrix<std::complex<double> > ident(dim_,dim_);  MOs::Identity(ident); //produces a dim by dim identity matrix
 PauliZ.initialize(dim_,dim_);
 PauliZ(1,1)=1;
 
	matrix<std::complex<double> > Htemp_(dim_,dim_), Htemp2_(dim_,dim_);	
	matrix<std::complex<double> >* lastRho = &rho_initial_;
   	matrix<std::complex<double> >* curUnitary = &ident;
 	std::complex<double>* freqs  = new std::complex<double>[dim_];

	lastRho = &rho_initial_;

   	curUnitary = &ident;
   	for(size_t j = 0; j < num_time_; j++){
   		Htemp_= *H_drift_;
			
	//any rotating frame, commented out to avoid seg faults			
		if(freqs_!=NULL)
 		{	for(size_t d = 0; d < dim_-1; ++d)
 				freqs[d] = exp(std::complex<double>(0.0,-freqs_[d]*(j*h_))); 
 			MOs::diagmult(Htemp_, freqs, Htemp_);
 			for(size_t d = 0; d < dim_-1; ++d)
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
		//Updates Unitaries	
 U_[j]=Unitaries_[j]*(*curUnitary);
 //updates rho to compare to/use
 *lastRho=Unitaries_[j]*(*lastRho)*MOs::Dagger(Unitaries_[j]);
 //sets our rho at each time step
Rho1_[j]=GenKrausDis(j, h_,lastRho);

   		curUnitary = &(U_[j]);//sets current unitary location
		lastRho=&(Rho1_[j]);//sets rho location
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
//forward propagate commute for unitaries, commented out due to errors
void Evolution::forwardpropagate_Commute()
{
  /*
  matrix<std::complex<double> > tempmat(dim_,dim_), ident(dim_,dim_), *propag;
  //	std::complex<double> scalarexponents[dim_];
  //	std::complex<double> freqs[dim_];
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


  */
}
	



#endif /*Evolution_h */

