/*
Name: Grape.h
Author: Jay M. Gambetta and Felix Motzoi

Dependences: MatrixExponential.hpp, QuantumOperations.hpp
Brief Discription: This program implements grape using my matrix class
Limitations: I dont see the need for both statetransfer and unitarytransfer, these should be combined into one.

Version History
	v0: June  9, 2008.
	v1: Feburary 5th, 2009 uses my new matrix class
	v2: Feburary 25th, 2009 edits to make adding different performance functions easier. It turns out it is a bit slower with the pointer to the function but not much
	v3: March ~1, 2009 added envelope for adding nonlinearity to the controls/time e.g. working in the lab frame
	v4: May 20, 2009 added some new functions

*/
/*Copy over:	GetFidelity
				GetCount
				StateTransfer
				GetPopulations
*/				
#ifndef OptimizeEvolution_h
#define OptimizeEvolution_h
#include "MatrixExponential.hpp"
#include "QuantumOperations.hpp"
#include "Evolution.hpp"
class OptimizeEvolution;

typedef double (OptimizeEvolution::* Fid)() const;
typedef double (OptimizeEvolution::* GradFid)(const size_t j, const size_t k);
typedef void (OptimizeEvolution::* ptrGradients)();

class OptimizeEvolution : public Evolution{
	public:
		OptimizeEvolution(size_t dim, size_t num_controls, size_t num_time, const double dt, char* filename);
		virtual ~OptimizeEvolution();

		void SetNumericalParameters(const double fidelity, const double base_a, const double epsilon, const double tolerance, size_t max_iter);
		void SetNumTimes(const size_t newnumtimes);
		
		//helper functions
		void computegradient();
		void computegradient_Commute();
		void updatecontrols();
		void unwindcontrols();
		
		//ensemble/robust versions
		void robustforwardpropagate();
		void robustgradient();
		double robustfidelity() const;
		
		//optimization function
		void UnitaryTransfer();	
		
		//Fidelity functions adapted from Khaneja GRAPE paper
		double Phi3() const;
		double GradPhi3(const size_t j, const size_t k) const;
		double Phi4() const;
		double GradPhi4(const size_t j, const size_t k) ;
		double Phi4TrCav() const;
		double GradPhi4TrCav(const size_t j, const size_t k) ;		
		double Phi4Sub2() const;
		double GradPhi4Sub2(const size_t j, const size_t k) ;
		
		double LeakagePenalty() const;
		double GradLeakagePenalty(const size_t j, const size_t k) ;
		double RelaxationPenalty() const;
		double GradRelaxationPenalty(const size_t j, const size_t k) ;
		
		//public variables, probably should be made protected
		size_t nconfigs_;					
		double config_;		
		
		GradFid gradPhi;
		Fid	Phi;
		GradFid gradPenalty;
		Fid	Penalty;
		
		OptimizeEvolution** evolutions;		//different evolutions for e.g. robust control or comparison purposes
		size_t num_evol;						//number such evolutions

		ptrGradients pGetGradient;

		size_t gatesteptime;				//when sweeping gate times increment gate time by this amount at each iteration
		size_t ngatetimes;							//number of gate times to sweep through
		double startfid;						//starting fidelity before sweep
		
		//Accessor methods
		double GetTopFidelity(){ return top_fidelity_; }
		size_t GetCount(){ return count_; }			
		void setFilename(char* newfilename) { strcpy(filename_, newfilename);  }
		void writecontrols(char* suffix);
		void writepopulations(char* suffix);
		void writefidelities(char* suffix);
		void writeconfigdata(char* suffix);
		
	public:
		size_t level0index_;				//the index of the 0 state for a qubit in a larger Hilbert space
		matrix<std::complex<double> > gradU;	//rho in the Khaneja paper
		matrix<std::complex<double> > lambda;	//lambda in the Khaneja paper
		double top_fidelity_;				// the highest achieved fidelity
		double fidelity_;					// the desired in fidelity 
		double tolerance_;					// the minimal change in fidelity 
		double base_a_;						// scale parameter for derivative, epsilon =base_a^power
		int power_;							//scales parameter for derivative, base_a^power
		double epsilon_;					//scales parameter for derivative, base_a^power
		double alpha_;		
							
		double **gradient_;					//an array containing the gradiant for each i and j
							
		size_t pos_count_;					// counts the number of consecutive forward steps in the algorithm (positive delta_fid)	
		size_t count_;						// counter used to count calls to grape
		size_t failcount_;					// counter used to count failed calls to grape	
		size_t max_iter_;					// max iteration befor it turns off
		char filename_[80];	
		
};							

#include "AnalyticControl.hpp"


inline OptimizeEvolution::OptimizeEvolution(size_t dim, size_t num_controls, size_t num_time, const double h, char* filename) 
	: Evolution(dim, num_controls, num_time, h), Penalty(NULL), gradPenalty(NULL){
	config_=0;
	if(filename) strcpy(this->filename_, filename);
	level0index_=0;
	num_evol=0;
	gradU.initialize(dim,dim);	
	lambda.initialize(dim,dim);	
	pGetGradient = 	static_cast<ptrPropagate>(&OptimizeEvolution::computegradient);		
	Phi = &OptimizeEvolution::Phi4;
	gradPhi = &OptimizeEvolution::GradPhi4;
	gradient_ = new double *[num_controls_];
	for(size_t k=0; k < num_controls_; ++k)
		gradient_[k] = new double[num_time_];
	char filecommand[80];
	system(strcat(strcpy(filecommand, "mkdir "), filename));	
}

OptimizeEvolution::~OptimizeEvolution(){
	for( size_t k = 0; k < num_controls_; ++k)
		delete [] gradient_[k];
	delete [] gradient_;
}

inline void OptimizeEvolution::SetNumericalParameters(const double fidelity, const double base_a, const double epsilon, const double tolerance, size_t max_iter){

	nconfigs_=1;
	tolerance_=tolerance;
	fidelity_=fidelity;
	tgate_=h_*num_time_;
	max_iter_ = max_iter;
	base_a_ = base_a;
	epsilon_=epsilon*h_;
	alpha_ = 0;//alpha*2.0*h;
	lastcontrol=0;
	
	if(verbose==yes)
	{
	 	std::cout << "--------------------Numerical Parameters--------------------" << std::endl;
		std::cout << "smallest change in phi allowed: " << tolerance_ << std::endl;
		std::cout << "desired fidelity: " << fidelity_ << std::endl;
		std::cout << "max iterations allow before termination: " << max_iter_ << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}

inline void OptimizeEvolution::SetNumTimes(const size_t newnumtimes)
{
	Evolution::SetNumTimes(newnumtimes);
	for( size_t k = 0; k < num_controls_; ++k)
	{	delete [] gradient_[k];	
		gradient_[k] = new double[num_time_];
	}
}



void OptimizeEvolution::robustgradient() 
{
	for(size_t k = 0; k < num_controls_; ++k)
			for(size_t j = 0; j < num_time_; ++j)
				gradient_[k][j]=0;	
	for (size_t s=0; s<num_evol; s++ )
	{
		evolutions[s]->computegradient();
		for(size_t k = 0; k < num_controls_; ++k)
			for(size_t j = 0; j < num_time_; ++j)
				gradient_[k][j]+= evolutions[s]->gradient_[k][j];
	}
}

void OptimizeEvolution::robustforwardpropagate()
{
	for (size_t s=0; s<num_evol; s++ )
	{
		for(size_t k = 0; k < num_controls_; ++k)
		{	evolutions[s]->theControls_[k]->Replicate(theControls_[k]);
			for(size_t j = 0; j < num_time_; ++j)
			{	
	//			cout << evolutions[s]->theControls_[k]->u_[j] << " " << theControls_[k]->u_[j] << endl;
			}
		}
		(evolutions[s]->*(evolutions[s]->pPropagate))();							
	}
}

double OptimizeEvolution::robustfidelity() const
{
	double fidel = 0;
	for (size_t s=0; s<num_evol; s++ )
		fidel+=(evolutions[s]->*evolutions[s]->Phi)();
	return fidel/num_evol;
}


void OptimizeEvolution::computegradient()
{
	matrix<std::complex<double> > Htemp_(dim_, dim_);
	lambda = MOs::Dagger(rho_desired_); //Lagrange multiplier
	double deriv;
	
	for(size_t k = 0; k < num_controls_; ++k)
			for(size_t j = 0; j < num_time_; ++j)
				gradient_[k][j]=0;	
	
	for(int j = num_time_-1; j >=0 ; --j){
			for(size_t k = 0; k < num_controls_; ++k){			
				theControls_[k]->getMatrixGradient(j,  &Htemp_); // k-th control,j-th pixel matrix control
				
				//integral of k,j-th control in diagonal frame of total hamiltonian
				Htemp_ = MOs::Dagger(Z[j])*Htemp_* Z[j];
				for (size_t q = 0; q < dim_; ++q)
					for (size_t p = 0; p < dim_; ++p)
						if(p!=q && (W[j][q]-W[j][p]))Htemp_(q,p) *= (i*cos(-(W[j][q]-W[j][p])*h_)-sin(-(W[j][q]-W[j][p])*h_)-i)/(W[j][q]-W[j][p])/h_;
				Htemp_ = Z[j]*Htemp_*MOs::Dagger(Z[j]);		
				
				gradU= Htemp_*U_[j]; // k,j-th matrix gradient of total matrix propagation
				deriv = (this->*gradPhi)(j,k);
				if(gradPenalty)
					deriv = deriv + (this->*gradPenalty)(j,k);
				theControls_[k]->setSubGradients(j, deriv , gradient_[k]);  //trace and chain rule to other (scalar) pixel controls
		 }
		 if(j) lambda = lambda*Unitaries_[j]; 		 
	}
		
}

void OptimizeEvolution::computegradient_Commute()
{
	int start, inc;
	matrix<std::complex<double> > Htemp_;
	lambda = MOs::Dagger(rho_desired_);		//Lagrange multiplier
	
	for(size_t k = 0; k < num_controls_; ++k)
			for(size_t j = 0; j < num_time_; ++j)
				gradient_[k][j]=0;	
	
	for(int j = num_time_-1; j >=0 ; --j){
		if(j%8==0 || j%8==3 || j%8==5 || j%8==6)	{	start=0; inc=+1; } 
		else		{	start=num_controls_-1; inc=-1;	}
		for(int k = start; k < num_controls_ && k>=0;  k+=inc){
			theControls_[k]->getMatrixGradient(j,  &Htemp_);	 // k-th control,j-th pixel matrix control
			gradU=Htemp_*UC_[k][j];								// k,j-th matrix gradient of total matrix propagation
			theControls_[k]->setSubGradients(j,  (this->*gradPhi)(j,k), gradient_[k]);	//trace and chain rule to other (scalar) pixel controls
			lambda = lambda*UnitariesC_[k][j];
		}				 
	}
}

inline void OptimizeEvolution::updatecontrols()
{
	 ++pos_count_;
	 if(pos_count_ ==3){		//if we get the correct direction 3 times we set power to be bigger
		power_++;
		pos_count_ =0;
	 }	
	 
	for(size_t k =0; k < num_controls_; ++k)
	  for(size_t j = 0; j < num_time_; ++j) 
		  theControls_[k]->u_[j]+=pow(base_a_,power_)*gradient_[k][j];
	 
}

inline void OptimizeEvolution::unwindcontrols()
{
	power_--;
	pos_count_ =0;
	
	for(size_t j = 0; j < num_time_; ++j){
		for(size_t k =0; k < num_controls_; ++k){
			theControls_[k]->u_[j]+=pow(base_a_,power_)*(1-base_a_)*gradient_[k][j];
		}
	}		
	failcount_++;
}


void OptimizeEvolution::UnitaryTransfer(){
	
	for(size_t k = 0; k < num_controls_+2; ++k)
		if(controlsetflag_[k]) UFs::MyError("Grape::UnitaryTransfer(): you have not set the drift and all AnalyticControl hamiltonians:\n");
	ofstream fidelout;
	char datafilename[80], scriptfilename[80];
	if(filename_!=NULL)
	{	strcat(strcpy(datafilename,filename_), "/fidels");
		strcpy(scriptfilename, datafilename);
		UFs::OpenFile(strcat(datafilename, ".dat"), fidelout, 16);
	}//Set the counters to zero	   
	size_t n_conseq_unimprov = count_ = failcount_ = pos_count_ = 0;
	
	double current_fidelity=0.0,  avgimprov = 0.0, delta_fidelity=1.0;
	writecontrols("controls_start");
	
	for(top_fidelity_=0.0, power_=0; abs(delta_fidelity) > tolerance_ && n_conseq_unimprov<20 && count_< max_iter_; count_++)
	{	
		(this->*pPropagate)();	
		current_fidelity=(this->*Phi)();
		if(Penalty)		
			current_fidelity=(current_fidelity+(this->*Penalty)());
		//for(int j=0; j<num_time_; j++)
		//		cout << U_[j].SetOutputStyle(Matrix) << endl;
				
		n_conseq_unimprov += (current_fidelity-top_fidelity_==delta_fidelity); //numerical precision error
		delta_fidelity=current_fidelity-top_fidelity_;
		
		if(delta_fidelity>0)  //the update in the controls			
		{
			top_fidelity_+=delta_fidelity;
			if(top_fidelity_ >= fidelity_) break;
			avgimprov+=delta_fidelity;
			(this->*pGetGradient)();
			updatecontrols();
		}
		else unwindcontrols();
		
		for(int k=0; k<num_controls_; k++)
				theControls_[k]->Interpolate();
		
		size_t echorate = 50;
		if (count_%echorate==0)
		{	std::cout << "step " << count_+1 << ", " << count_-failcount_+1 << " successes, " << failcount_ << " failures, with error " << 1-top_fidelity_ << " and average improvement " << avgimprov/echorate <<  std::endl; 			  
			fidelout << count_+1 << '\t' << 1-current_fidelity << '\t' << pow(base_a_,power_) << endl;	
			avgimprov=0;
			writecontrols("controls_int");
			writepopulations("popul_int");	
		}
   }
   
    if(count_>= max_iter_) cout << "max iterations: ";
	else if(abs(delta_fidelity) <= tolerance_) cout << "tolerance: ";
	else if(top_fidelity_ >= fidelity_) cout << "success: ";
	else cout << "numerical precision: ";
   
	std::cout << "finished on step " << count_ << ", " << count_-failcount_ << "successes, " << failcount_ << " failures, with error " << 1-top_fidelity_ << " and delta fidelity " << delta_fidelity << std::endl;  
	
	fidelout << count_+1 << '\t' << 1-top_fidelity_ << '\t' << pow(base_a_,power_) << endl;
	fidelout.close();
	USs::MakeGnuplotScript(2, scriptfilename, "", "", datafilename, true);
	
	char command[80];
	strcat(strcat(strcpy(command, "gnuplot "), scriptfilename),".p");
	system(command);
					
	writecontrols("controls_final");
	writepopulations("popul_final");
}


inline double OptimizeEvolution::Phi3() const{
	//the measure implemented, phi (PHI_3) is phi = trace[U_desired* UN-1....U_0]/D
	std::complex<double> temp1=0.0;
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += std::conj(rho_desired_(p,q))*U_[num_time_-1](p,q);
		}
	}
	return std::real(temp1)/dim_;
}

inline double OptimizeEvolution::GradPhi3(const size_t j, const size_t k) const{
	//phi_3
	std::complex<double> temp1=0;
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += lambda(q,p)*gradU(p,q);	
		}
	}	
	return epsilon_*std::imag(temp1);
}
inline double OptimizeEvolution::Phi4() const{
	//the measure implemented, phi (PHI_4) is phi = |trace[U_desired* UN-1....U_0 ]|^2/D^2
	std::complex<double> temp1=0.0;
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p)
			temp1 += std::conj(rho_desired_(p,q))*U_[num_time_-1](p,q);		
	} 
	return real(temp1*std::conj(temp1))/dim_/dim_;
}
inline double OptimizeEvolution::GradPhi4(const size_t j, const size_t k) {
	//phi_4
	//2.0*h_*imag(QOs::Expectation( MOs::Dagger(L),H*R)*QOs::Expectation(MOs::Dagger(R),L));
	std::complex<double> temp1=0, temp2=0;
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += lambda(q,p)*gradU(p,q);
			temp2 += std::conj(U_[num_time_-1](p,q))*rho_desired_(p,q);	
		}
	}
	return 2.0*epsilon_*std::imag(temp1*temp2);//*one_on_dim_*one_on_dim_;
}
inline double OptimizeEvolution::Phi4TrCav() const{
	//the measure implemented, phi (PHI_4) is phi = |trace[U_desired* UN-1....U_0 ]|^2/D^2
	std::complex<double> temp1=0.0;
	size_t dimcav = 3;
	matrix<complex<double> > UA = (1.0/dimcav)*MOs::TraceOutA(U_[num_time_-1], dimcav);
	
	U_[num_time_-1].SetOutputStyle(Matrix);
	UA.SetOutputStyle(Matrix);
	cout << "fid trA" << endl;
	cout << U_[num_time_-1] << endl << UA << endl;
	
	
	size_t dimq = dim_/dimcav;
	for(size_t q=0; q< dimq; ++q){
		for(size_t p=0; p<dimq; ++p)
			temp1 += std::conj(rho_desired_(p,q))*UA(p,q);		
	}
	//return real(temp1*std::conj(temp1))/dim_/dim_;
	return real(temp1*std::conj(temp1))/4/4;
}
inline double OptimizeEvolution::GradPhi4TrCav(const size_t j, const size_t k) {
	//phi_4
	//2.0*h_*imag(QOs::Expectation( MOs::Dagger(L),H*R)*QOs::Expectation(MOs::Dagger(R),L));
	std::complex<double> temp1=0, temp2=0;
	size_t dimcav = 3;
	matrix<complex<double> > UA = (1.0/dimcav)*MOs::TraceOutA(U_[num_time_-1], dimcav);
	matrix<complex<double> > gradUA = (1.0/dimcav)*MOs::TraceOutA(gradU, dimcav);
	matrix<complex<double> > lamA = (1.0/dimcav)*MOs::TraceOutA(lambda, dimcav);
	size_t dimq = dim_/dimcav;
	
	for(size_t q=0; q< dimq; ++q){
		for(size_t p=0; p<dimq; ++p){
			temp1 += lamA(q,p)*gradU(p,q);
			temp2 += std::conj(U_[num_time_-1](p,q))*rho_desired_(p,q);	
		}
	}
	return 2.0*epsilon_*std::imag(temp1*temp2);//*one_on_dim_*one_on_dim_;
}


inline double OptimizeEvolution::Phi4Sub2() const{
	//the measure implemented, phi (PHI_4_sub) is phi = |sum_{i=0,1} <i|U_desired* UN-1....U_0|i>|^2/2^2
	std::complex<double> temp1=0.0;
	for(size_t q=level0index_; q< level0index_+2; ++q){
		for(size_t p=0; p<dim_; ++p){
		temp1 += std::conj(rho_desired_(p,q))*U_[num_time_-1](p,q);		
		}
	}
	return std::real(temp1*std::conj(temp1))*0.25;
}

inline double OptimizeEvolution::GradPhi4Sub2(const size_t j, const size_t k){
	//Phi_4_subsystem
	//2.0*h_*imag(QOs::Expectation( MOs::Dagger(L),H*R)*QOs::Expectation(MOs::Dagger(R),L));
	std::complex<double> temp1=0, temp2=0;
	for(size_t q=level0index_; q< level0index_+2; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += lambda(q,p)*gradU(p,q);
			temp2 += std::conj(U_[num_time_-1](p,q))*rho_desired_(p,q);	
		}
	}
	return 2.0*epsilon_*std::imag(temp1*temp2);//*one_on_dim_*one_on_dim_;
}
//*
inline double OptimizeEvolution::LeakagePenalty() const{
	//the measure implemented, phi (PHI_4_sub) is phi = |sum_{i=0,1} <i|U_desired* UN-1....U_0|i>|^2/2^2
	std::complex<double> ret=0;
	for(size_t j=0; j<num_time_; j++) 
	{	std::complex<double> temp1=0.0;
		for(size_t q=level0index_; q< level0index_+2; ++q){
			for(size_t p=level0index_; p<level0index_+2; ++p){
			temp1 += std::conj(U_[j](p,q))*U_[j](p,q);		
			}
		}
		ret+=temp1*0.5;
	}	
	return 0.001*std::real(ret)/num_time_;
}

inline double OptimizeEvolution::GradLeakagePenalty(const size_t j, const size_t k){
	std::complex<double> ret=0;
	//matrix<complex<double> > Htemp_(dim_, dim_);
	//theControls_[k]->getMatrixGradient(j,  &Htemp_);
	matrix<complex<double> > UA=gradU;
	//UA = Htemp_*Unitaries_[j];
	for(size_t j=0; j<num_time_; j++) 
	{	std::complex<double> temp1=0.0, temp2=0;		
		for(size_t q=level0index_; q< level0index_+2; ++q){
			for(size_t p=level0index_; p<level0index_+2; ++p){
				temp1 += std::conj(U_[j](p,q))*UA(p,q);	
			}
		}
		if(j<num_time_-1) UA = Unitaries_[j+1]*UA;
		ret+=temp1;
	}	
	return 0.0004*std::imag(ret)/num_time_;
}/*/

inline double OptimizeEvolution::LeakagePenalty() const{
	//the measure implemented, phi (PHI_4_sub) is phi = |sum_{i=0,1} <i|U_desired* UN-1....U_0|i>|^2/2^2
	std::complex<double> ret=0;
	for(size_t j=0; j<num_time_; j++) 
	{	std::complex<double> temp1=0.0;
		for(size_t q=level0index_; q< level0index_+2; ++q){
			for(size_t p=level0index_; p<level0index_+2; ++p){
			temp1 += std::conj(Unitaries_[j](p,q))*Unitaries_[j](p,q);		
			}
		}
		ret+=temp1*0.5;
	}	
	return 0.01*std::real(ret)/num_time_;
}

inline double OptimizeEvolution::GradLeakagePenalty(const size_t j, const size_t k){
	std::complex<double> ret=0;
	matrix<complex<double> > Htemp_(dim_, dim_);
	theControls_[k]->getMatrixGradient(j,  &Htemp_);
	matrix<complex<double> > UA=gradU;
	UA = Htemp_*Unitaries_[j];
	//for(size_t j=0; j<num_time_; j++) 
	{	std::complex<double> temp1=0.0, temp2=0;		
		for(size_t q=level0index_; q< level0index_+2; ++q){
			for(size_t p=level0index_; p<level0index_+2; ++p){
				temp1 += std::conj(Unitaries_[j](p,q))*UA(p,q);	
			}
		}
	//	if(j<num_time_-1) UA = UA * Unitaries_[j+1];
		ret+=temp1;
	}	
	return 0.004*std::imag(ret)/num_time_;
}*/



inline double OptimizeEvolution::RelaxationPenalty() const{
	//the measure implemented, phi (PHI_4_sub) is phi = |sum_{i=0,1} <i|U_desired* UN-1....U_0|i>|^2/2^2
	size_t dimq = 2;
	matrix<complex<double> > initrho(dim_, dim_), a(dim_/dimq,dim_/dimq), IdentQ(dimq,dimq);
	 MOs::Identity(IdentQ);
	//initrho=(1/dim_)*MOs::Identity(initrho); //maximally mixed state
	initrho(0,0)=initrho(0,1)=initrho(1,0)=initrho(1,1)=0.5;
	std::complex<double> ret=0;
	MOs::Destroy(a);
	a = (h_/100) * MOs::TensorProduct(a, IdentQ);
	matrix<complex<double> > ad = MOs::Dagger(a);
	
	for(size_t j=0; j<num_time_; j++) 
	{
		initrho = Unitaries_[j] * initrho * MOs::Dagger(Unitaries_[j]);
		initrho += a * initrho * ad - 0.5 * ad * a * initrho - 0.5 * initrho * ad * a; 
	}	
//	cout << initrho.SetOutputStyle(Matrix) << endl ;
	initrho=initrho*initrho;
//	cout << initrho.SetOutputStyle(Matrix) << endl ;
//	cout << initrho(0,0)+initrho(1,1) << endl;
	return 1-pow(abs(initrho(0,0)+initrho(1,1)),2);
}

inline double OptimizeEvolution::GradRelaxationPenalty(const size_t j, const size_t k){
	std::complex<double> ret=0;
	matrix<complex<double> > Htemp_(dim_, dim_);
	theControls_[k]->getMatrixGradient(j,  &Htemp_);
	matrix<complex<double> > UA=gradU;
	//UA = Htemp_*Unitaries_[j];
	for(size_t j=0; j<num_time_; j++) 
	{	std::complex<double> temp1=0.0, temp2=0;		
		
		
		
		if(j<num_time_-1) UA = Unitaries_[j+1]*UA;
		ret+=temp1;
	}	
	return std::imag(ret)/num_time_;
}




void OptimizeEvolution::writecontrols(char* suffix="")
{
	char scriptfilename[80], datafilename[80];
	ofstream dataout;	
	if(filename_!=NULL && suffix!=NULL)
		strcat(strcat(strcpy(datafilename,filename_),"/"), suffix);		
	else return;
	strcpy(scriptfilename, datafilename);
	UFs::OpenFile(strcat(datafilename, ".dat"),dataout, 16);
	
	USs::MakeGnuplotScript(num_controls_, scriptfilename, "", "", datafilename, false);
	
	for(size_t j =0; j < num_time_; j++){
		dataout <<  h_*(j+0.5);
		for(int k=0; k<num_controls_; k++)
				dataout << "\t" << theControls_[k]->u_filt[j];
		dataout << "\n";		
	}
	
	dataout.close();
	
	char command[80];
	strcat(strcat(strcpy(command, "gnuplot "), scriptfilename),".p");
	system(command);
}
/*	char scriptfilename[80], datafilename[80], command[80];
	ofstream dataout;	
	if(filename_==NULL && suffix==NULL)
		return;
	
	system(strcat(strcpy(command, "cd "), filename_));
	system("pwd");
	
	strcpy(datafilename,suffix);
	strcat(datafilename, ".dat");
	strcpy(scriptfilename,suffix);
	USs::MakeGnuplotScript(num_controls_, scriptfilename, "", "", datafilename);
	UFs::OpenFile(datafilename,dataout, 16);
	
	for(size_t j =0; j < num_time_; j++){
		dataout <<  h_*(j+0.5);
		for(int k=0; k<num_controls_; k++)
				dataout << "\t" << theControls_[k]->u_filt[j];
		dataout << "\n";		
	}
	
	dataout.close();
	
	strcat(strcat(strcpy(command, "gnuplot "), scriptfilename),".p");
	system(command);
	system("cd ..");
	system("pwd");
*/

void OptimizeEvolution::writepopulations(char* suffix="")
{
	char scriptfilename[80], datafilename[80],startpop[2];
	ofstream dataout;	
	strcpy(startpop, "__");
	for(size_t l=0; l<2; l++)
	{
		if(filename_!=NULL && suffix!=NULL)
			strcat(strcat(strcpy(datafilename,filename_),"/"), suffix);		
		else return;
		startpop[1]='0'+l;
		strcat(datafilename, startpop);
		strcpy(scriptfilename, datafilename);
		
		UFs::OpenFile(strcat(datafilename, ".dat"),dataout, 16);		
		USs::MakeGnuplotScript(dim_, scriptfilename, "", "", datafilename, false);
				
		
		for(size_t j =0; j < num_time_; j++){  //loop thru dim
			dataout <<  h_*(j+0.5);
			for(size_t k=0; k<dim_; k++)
				dataout<< "\t" << pow(abs(U_[j](l,k)),2);
			dataout << "\n";		
		}
		
		dataout.close();
		
		char command[80];
		strcat(strcat(strcpy(command, "gnuplot "), scriptfilename),".p");
		system(command);
	}
}

void OptimizeEvolution::writefidelities(char* suffix="")
{

}

void OptimizeEvolution::writeconfigdata(char* suffix="")
{
	//dataout << Hdrift
	//control hams
	//tgate
	//target fid
	//target gate
	//nipix, nsubpix
	//initial cond? redund?
}


#endif /* OptimizeEvolution_h */

