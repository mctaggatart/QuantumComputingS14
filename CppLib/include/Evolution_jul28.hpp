
/*
Name: Grape.h
Author: Jay M. Gambetta and Felix Motzoi
			
Dependences: MatrixExponential.hpp, QuantumOperations.hpp
*/

#ifndef Evolution_h
#define Evolution_h
	
class Control;
class Evolution;
	
typedef void (Evolution::* ptrPropagate)();
	
#include "QuantumOperations.hpp"
class Evolution {
	public:
		Evolution(size_t dim, size_t num_controls, size_t num_time, const double h);// number of controls, the number of time points
		virtual ~Evolution();
	
		void SetNumTimes(const size_t newnumtimes);
		
		matrix<std::complex<double> >* GetHdrift();
		void SetHdrift( matrix<std::complex<double> >& H_drift);
		
		matrix<std::complex<double> > GetRhoInitial();
		void SetRhoInitial(const matrix<std::complex<double> >& rho_initial);
		
		Control* GetuControl(const size_t k);
		void Setucontrol(Control* control, const int index=-1);
		
		matrix<std::complex<double> > GetRhoDesired() { return rho_desired_; }
		void SetRhoDesired(const matrix<std::complex<double> >& rho_desired);

		void Evolution::writepopulations(char* outfile);
						
		ptrPropagate pPropagate;

		void forwardpropagate();
		void forwardpropagate_Taylor();
		void forwardpropagate_CommuteAlt();
		void forwardpropagate_Commute();
		
		matrix<std::complex<double> >* U_; //an array of the forward evolutions (rho in paper)
		matrix<std::complex<double> >** UC_; //an array of the forward evolutions (rho in paper) splitting by control
		double* freqs_;
						
	protected:
		size_t num_time_;//  the number of time points	
		double h_; //the time step
		double tgate_;
		size_t dim_;//dimensions of the Hilbert space
		double one_on_dim_;//used alot so I have defined it
		
		matrix<std::complex<double> > *H_drift_; //drift hamiltonian
		matrix<std::complex<double> > H_drift_exp; //exponential of drift hamiltonian
		matrix<std::complex<double> > rho_desired_; //the desired rho or unitary			//rename as U_desired_
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
	Z = new matrix<std::complex<double> >[num_time_];	
	W = new double *[num_time_];
	for(size_t j = 0; j < num_time_; ++j)
	{ Z[j].initialize(dim_,dim_);
	  W[j] = new double[dim_];	
	}
	pPropagate = &Evolution::forwardpropagate;
	theControls_ = new Control*[num_controls_];
	controlsetflag_=new size_t[num_controls_+3]; // one more for Hdrift and one more for rho_desired and rho_initial
	UC_ = new  matrix<std::complex<double> >*[num_controls_];
	UnitariesC_ = new  matrix<std::complex<double> >*[num_controls_];
	for(size_t k=0; k < num_controls_; ++k){
		//controls_[k] = new double[num_time_];
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
	Unitaries_=new  matrix<std::complex<double> >[num_time_];
	for(size_t j=0; j < num_time_; ++j){
		U_[j].initialize(dim_,dim_);
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
		//U_[j].clear();
		//Unitaries_[j].clear();
	//	cout << W[j] << " ";
		delete [] W[j];
	}
	delete [] UC_;
	delete [] UnitariesC_;
	delete [] theControls_;
	delete [] controlsetflag_;
	delete [] U_;
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
	delete [] Unitaries_;
	delete []Z;
	delete []W;
	
	U_=new matrix<std::complex<double> >[num_time_];
	Unitaries_=new  matrix<std::complex<double> >[num_time_];
	for(size_t j=0; j < num_time_; ++j){
		U_[j].initialize(dim_,dim_);
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

inline void Evolution::Setucontrol(Control* control, const int k){
	if(k>=0) lastcontrol=k;
	
	theControls_[lastcontrol]=control;
	
	if(!controlsetflag_[num_controls_])
		theControls_[lastcontrol]->ZeH= MOs::Dagger(theControls_[lastcontrol]->Z)*H_drift_exp;
	
	controlsetflag_[lastcontrol++]=0;
}

		
//Set the target state for transfer
inline void Evolution::SetRhoDesired(const matrix<std::complex<double> >& rho_desired){
	rho_desired_=rho_desired;
	controlsetflag_[num_controls_+1]=0;
	if(verbose==yes)
	{
		rho_desired_.SetOutputStyle(Matrix);
	 	std::cout << "--------------------rho desired is set---------------------------" << std::endl;
		std::cout << rho_desired_ << std::endl;
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


	
inline void Evolution::writepopulations(char* outfile)
{	
	ofstream popout;
	UFs::OpenFile(outfile,popout, 16);
	
	matrix<std::complex<double> > pop = rho_initial_;			
							
	for(size_t j =0; j < num_time_; j++){
		popout << h_* j;
		pop = U_[j]*rho_initial_*MOs::Dagger(U_[j]);
		for(int d=0; d<dim_; d++)
			popout << '\t' << pop(d,d);
		popout << endl;
	}
	popout.close();
}


void Evolution::forwardpropagate()
{
	matrix<std::complex<double> > ident(dim_,dim_); 
	matrix<std::complex<double> > Htemp_(dim_,dim_), Htemp2_(dim_,dim_);
	MOs::Identity(ident);
	matrix<std::complex<double> >* lastrho = &ident; 
	
	std::complex<double> freqs[dim_];
	// the propagators and the foraward evolution
	lastrho = &ident;
	
	for(size_t j = 0; j < num_time_; ++j){		
		U_[j].SetOutputStyle(Matrix);
		Htemp_.SetOutputStyle(Matrix);
		Unitaries_[j].SetOutputStyle(Matrix);
	//	cout << j << endl << U_[j] << endl;

		for(size_t d = 0; d < dim_; ++d)
			freqs[d] = exp(std::complex<double>(0.0,-freqs_[d]*(j*h_)));
		MOs::diagmult(Htemp_, freqs, *H_drift_);
		for(size_t d = 0; d < dim_; ++d)
			freqs[d] = exp(std::complex<double>(0.0,freqs_[d]*(j*h_)));
		MOs::multdiag(Htemp_, Htemp_, freqs);
	
		for(size_t k = 0; k < num_controls_; ++k){	
			theControls_[k]->getMatrixControl(j,  &Htemp2_);
			Htemp_ = Htemp_ + Htemp2_;
		}
	//	cout << "Htemp " << Htemp_ << endl;
		times0-=clock();
		Unitaries_[j]= ExpM::EigenMethod(Htemp_,-i*h_, &(Z[j]), W[j]);	
	//	cout << "Unt " << Unitaries_[j] << endl;
		U_[j] = Unitaries_[j]*(*lastrho);
	times0+=clock();
			//	cout <<j << endl <<  U_[j] << endl;
	//	if(j==516) exit(1);
		lastrho=&(U_[j]);	
	}
}	

void Evolution::forwardpropagate_Taylor()
{
	static matrix<std::complex<double> > ident(dim_,dim_); 
	MOs::Identity(ident);
	matrix<std::complex<double> >* lastrho = &ident; 
	matrix<std::complex<double> > Htemp_, Htemp2_;
	lastrho = &ident;
	
	
	for(size_t j = 0; j < num_time_; ++j){			
		U_[j] = (*lastrho);
		Htemp_ = *H_drift_;
		for(size_t k = 0; k < num_controls_; ++k){	
			theControls_[k]->getMatrixControl(j,  &Htemp2_);
			for(size_t q=0; q< dim_; ++q){
				for(size_t p=0; p<dim_; ++p){
					Htemp_(q,p) += Htemp2_(q,p);
				}
			}
		}
		Htemp_=0.125*Htemp_; 
		MOs::Identity(Unitaries_[j]);	
		Unitaries_[j] = Unitaries_[j] + (-i*h_)*Htemp_ + 0.5*(-i*h_)*(-i*h_)*Htemp_*Htemp_ + 0.5*0.33333333333*(-i*h_)*(-i*h_)*(-i*h_)*Htemp_*Htemp_*Htemp_;
	
		Unitaries_[j] = Unitaries_[j]*Unitaries_[j];
		Unitaries_[j] = Unitaries_[j]*Unitaries_[j];
		Unitaries_[j] = Unitaries_[j]*Unitaries_[j];
		
		U_[j] = Unitaries_[j]*(*lastrho);
		lastrho=&(U_[j]);
	}	
}		

void Evolution::forwardpropagate_CommuteAlt()
{
	static matrix<std::complex<double> > ident(dim_,dim_); 
	MOs::Identity(ident);
	matrix<std::complex<double> > Htemp_;
	matrix<std::complex<double> >* lastrho = &ident; 
	// the propagators and the foraward evolution
	lastrho = &ident;	
	
	for(size_t j = 0; j < num_time_; ++j){			
		U_[j] = (*lastrho);
		
		if(j%2){
			Unitaries_[j]= ExpM::EigenMethod(*H_drift_,-i*h_, &(Z[j]), W[j]);
			U_[j] = Unitaries_[j]*U_[j];
			for(size_t k = 0; k < num_controls_; ++k){	
				theControls_[k]->getMatrixControl(j,  &Htemp_);
				Unitaries_[j]= ExpM::EigenMethod(Htemp_,-i*h_, &(Z[j]), W[j]);
				U_[j] = Unitaries_[j]*U_[j];				
			}
		} else {
			for(int k = num_controls_-1; k >=0 ; --k){	
				theControls_[k]->getMatrixControl(j,  &Htemp_);
				Unitaries_[j]= ExpM::EigenMethod(Htemp_,-i*h_, &(Z[j]), W[j]);
				U_[j] = Unitaries_[j]*U_[j];			
			}
					
			Unitaries_[j]= ExpM::EigenMethod(*H_drift_,-i*h_, &(Z[j]), W[j]);
			U_[j] = Unitaries_[j]*U_[j];
		}
		lastrho=&(U_[j]);
	}	
}		

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
	//		tempmat2 = H_drift_exp * *propag;
	//		propag=&tempmat2;
	//	}
			
		for(int k = start; k < num_controls_ && k>=0;  k+=inc){				
			for (size_t p=0; p < dim_; p++) 	
				scalarexponents[p] = std::exp(-i*h_*theControls_[k]->u_[j]*theControls_[k]->W[p]);
			MOs::multdiag(tempmat,  theControls_[k]->Z, scalarexponents);
			
			times0-=clock();
			UnitariesC_[k][j] = tempmat*theControls_[k]->ZeH;
			times0+=clock();
			
			for(size_t d = 0; d < dim_; ++d)
					freqs[d] = exp(std::complex<double>(0.0,-freqs_[d]*(j*h_)));
				MOs::diagmult(UnitariesC_[k][j], freqs, UnitariesC_[k][j]);
				for(size_t d = 0; d < dim_; ++d)
					freqs[d] = exp(std::complex<double>(0.0,freqs_[d]*(j*h_)));
				MOs::multdiag(UnitariesC_[k][j], UnitariesC_[k][j], freqs);
			
				
			
			times0-=clock();
			UC_[k][j] = UnitariesC_[k][j] * *propag;
			times0+=clock();
			
			propag = &(UC_[k][j]);
		}
		U_[j] = *propag;
	//	if((j)%2==0)// || (j)%4==3) 
	//	{
	//		U_[j] = H_drift_exp * *propag;
	//		propag = &U_[j];
	//	}
		//		UnitariesC_[k][j] =UnitariesC_[k][j]*H_drift_exp;
	
	}
}	



#endif /*Evolution_h */

