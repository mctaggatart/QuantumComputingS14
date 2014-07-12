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
#ifndef Grape_h
#define Grape_h
#include "MatrixExponential.hpp"
#include "QuantumOperations.hpp"
#include "ControlOptimizer.hpp"

class Grape : public ControlOptimizer{
	public:
		Grape(size_t dim, size_t num_controls, size_t num_time);// number of controls, the number of time points
		virtual ~Grape();
	
		void SetPhysicalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, double qubit_freq_, double anharmonicity_);
		
		void UnitaryTransfer(double (Grape::* phi)() const, double (Grape::* gradphi)(const size_t j, const size_t k) , 
				void (Grape::* ptrSetHnonlinearity)(size_t j), void (Grape::* ptrSetHgradnonlin)(const size_t j));
		void UnitaryTransfer(double (Grape::* phi)() const, double (Grape::* gradphi)(const size_t j, const size_t k) )
		{   UnitaryTransfer(phi, gradphi, &Grape::SetLinear, &Grape::SetGradLinear);
		}
		
		//getshortestcontrols, getcontrols, getrobustcontrols
		void sweeptimes(const char* outfile); //should pass ptr to one-iteration function, then also add "Ascent()" fn
		void sweeptimesandcompare();
		void sweepenergies(const double min, const double max, const matrix<std::complex<double> >& Shift, const char* outfile); //should pass ptr to one-iteration function, then also add "Ascent()" fn
		void averageconfigurations();
		
		double Phi3() const;
		double GradPhi3(const size_t j, const size_t k) const;
		double Phi4() const;
		double GradPhi4(const size_t j, const size_t k) const;
		double Phi4TraceCav() const;
		double GradPhi4TraceCav(const size_t j, const size_t k) const;
		double Phi4Sub2() const;
		double GradPhi4Sub2(const size_t j, const size_t k);
		
		double GradPhi4Sub2Filter(const size_t j, const size_t k) const;
		double GradPhi4Sub2RWFilter(const size_t j, const size_t k) const;
		void SetLinear(size_t j);
		void SetGradLinear(size_t j);

		void SetLab2Quadratures(size_t j);
		void SetGradLab2Quadratures(size_t j);
		void SetCounterRotating(size_t j);
		void SetGradCounterRotating(size_t j);
		void SetCounterRotatingFilt(size_t j);
		void SetGradCounterRotatingFilt(size_t j);
		void SetLab2QuadFiltered(size_t j);
		void SetGradLab2QuadFiltered(size_t j);		
	
		size_t nconfigs_;
		
	protected:
		
		//system specific parameters:
		double qubit_freq_;
		double anharmonicity_;
		double drive_freq_;
		double filter_norm_;
		double filter_rate_;
		double config_;		
		
};

inline Grape::Grape(size_t dim, size_t num_controls, size_t num_time) : ControlOptimizer(dim, num_controls, num_time){
	
	config_=0;
	drive_freq_=0;
	
}

Grape::~Grape(){

}

inline void Grape::SetPhysicalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, 
											double qubit_freq, double anharmonicity)	
{	qubit_freq_=qubit_freq;
	anharmonicity_=anharmonicity;	
	nconfigs_=1;
	SetNumericalParameters(fidelity, h, base_a, epsilon, tolerance, max_iter);
}

void Grape::UnitaryTransfer(double (Grape::* phi)() const, double (Grape::* gradphi)(const size_t j, const size_t k), 
						void (Grape::* ptrSetHnonlinearity)(size_t j), void (Grape::* ptrSetHgradnonlin)(const size_t j)){
	static matrix<std::complex<double> > ident(dim_,dim_); 
	MOs::Identity(ident);
	matrix<std::complex<double> >* lastrho = &ident; 
		
	//Some flags to make sure all parameters are set
	size_t test=0;
	for(size_t k = 0; k < num_controls_+2; ++k)
	{
		test+=controlsetflag_[k];
		if(verbose==yes)
		{
		 	std::cout << k << "th control flag is " << controlsetflag_[k] << std::endl;
		}
	}
	if(test != 0)
		UFs::MyError("Grape::UnitaryTransfer(): you have not set the drift and all control hamiltonians:\n");
			
	//Set the counters to zero
	count_ =0;
	failcount_ =0;
	consec_failcount_=0;
	pos_count_ =0;
	//the power scale for epsilot (epsilon = base_a^power_)
	power_ =0;
	// fidelity ranges from 0 to 1 with 0 be orthogonal and 1 being the same time
	double current_fidelity=0.0;
	double delta_fidelity=1;
	top_fidelity_=0.0;	//Starting the algorithim	
	while(delta_fidelity*delta_fidelity > tolerance_ ) //TODO: min power setting, add success condition
	{			
		if(count_> max_iter_){	
			std::cout << "Grape::UnitaryTransfer(), algorithim did not converge in the allow max iterations" << std::endl;
			break;
		}
		count_++;
		if(num_controls_>2) drive_freq_ = qubit_freq_;// - controls_[2][1];	
		current_fidelity=0;
		
		for(size_t k =0; k < num_controls_; ++k){
		  for(size_t j = 0; j < num_time_; ++j)
		  {	tempgradient_[k][j]=0;
		  }
		}	
		

	//	for(phase_=0; phase_<1.99*M_PI; phase_+=2*M_PI/nphases) 
		for(config_=0; config_<nconfigs_; config_++) 
		{	
			// the propagators and the foraward evolution
			lastrho = &ident;


			for(size_t j = 0; j < num_time_; ++j){			
				(this->*ptrSetHnonlinearity)(j);// updates H_controls_tdep_ for time dependent controls where time dependence is define by ptrSetHnonlinearity (hard coded)
				Htemp_ = H_drift_;
				for(size_t k = 0; k < num_controls_; ++k){				
					for(size_t q=0; q< dim_; ++q){
						for(size_t p=0; p<dim_; ++p){
							Htemp_(q,p) += H_controls_tdep_[k](q,p);
						}
					}
				}			
				//std::cout << Htemp_ << std::endl;	
			//	if(iscrap && count_>=38846) cout <<  "eig\n";
				Unitaries_[j]= ExpM::EigenMethod(Htemp_,-i*h_);
			//	if(iscrap && count_>=38846) cout <<  "eiged\n";
				rho_[j] = Unitaries_[j]*(*lastrho);
				lastrho=&(rho_[j]);
			}	
	//		if(iscrap) cout <<  "united\n";
			//Unitaries_[0].SetOutputStyle(Matrix);
			//std::cout << Unitaries_[0] << std::endl;
			//the average configuration fidelity
			current_fidelity+=(this->*phi)();
			//	if(iscrap) cout <<  "phied\n";
			// the backward evolution
			lambda_[num_time_-1] = rho_desired_;
			for(size_t j = num_time_-1; j > 0; --j){
		//		if(iscrap && count_>=38846) cout <<  "lamb\n";
				lambda_[j-1] = MOs::Dagger(Unitaries_[j])*lambda_[j];
			}
		//	if(iscrap) cout <<  "dagged\n";
			double tempgrad;
			for(size_t k =0; k < num_controls_; ++k){
				for(size_t j = 0; j < num_time_; j+=nsubpixels_[k]){
					tempgrad=0;
					for(size_t l = 0; l < nsubpixels_[k]; ++l){
						(this->*ptrSetHgradnonlin)(j+l);
						tempgrad+=(this->*gradphi)(j+l,k);
					}
					//		std::cout << "wha\n";
					for(size_t l = 0; l < nsubpixels_[k]; ++l)
					{	tempgradient_[k][j+l]+=tempgrad;
					}
				}
			} 								
		}	 
		current_fidelity/=nconfigs_;
		delta_fidelity=current_fidelity-top_fidelity_;
			
		if(delta_fidelity>0)  //the update in the controls			
		{
			// std::cout <<count_ << std::endl;	
				 
			 top_fidelity_+=delta_fidelity;
			 ++pos_count_;//if we get to correct directions we set power to be bigger
			 consec_failcount_=0;	
			
			 if(top_fidelity_ >= fidelity_)
			{	
				std::cout << top_fidelity_ <<"success, ";
				break;
			}
		
			 if(pos_count_ ==2){
				power_++;
				pos_count_ =0;
			 }
			 for(size_t k =0; k < num_controls_; ++k){
			  for(size_t j = 0; j < num_time_; ++j)
			  {	 gradient_[k][j]=tempgradient_[k][j];
				 controls_[k][j]+=pow(base_a_,power_)*gradient_[k][j];  
			  }
			 }
			 if (count_%200==1)
					std::cout << "step " << count_ << ", " << failcount_ << " failures, with fidelity " << top_fidelity_ << " and delta fidelity " << delta_fidelity << std::endl; 			  
		}
		else {
				power_--;
				pos_count_ =0;
	//			std::cout << controls_[0][2]-pow(base_a_,power_)*(base_a_)*gradient_[0][2] << "+" <<pow(base_a_,power_)*(1) << '*'  << gradient_[0][2]<< " ,fid="<<current_fidelity << " \n " ;
				
				for(size_t j = 0; j < num_time_; ++j){
					for(size_t k =0; k < num_controls_; ++k){
						controls_[k][j]+=pow(base_a_,power_)*(1-base_a_)*gradient_[k][j];
					}
				}		
				failcount_++;
				consec_failcount_++;		
		}
	
  }
		
	std::cout << "finished on step " << count_ << ", " << failcount_ << " failures, with fidelity " << top_fidelity_ << " and delta fidelity " << delta_fidelity << std::endl;  
}

void Grape::sweeptimes(const char* outfile){  //datafile
	
	ofstream dataout;
	UFs::OpenFile(outfile,dataout, 16);
	
	size_t rwa_num_time=theControls_[0]->npixels_;
	double lambda = sqrt(3.0/2);
	cout <<lambda*lambda << endl;
	double rab = M_PI;
	double lamb2 =  1/sqrt(2);
	double Delta = anharmonicity_;
	double Delta2 = anharmonicity_;
	H_drift_(0,0)=Delta2-Delta;	
	cout << theControls_[0]->Hcontrol_ << endl;
	for(size_t itime=0; itime<4; itime++, rwa_num_time+=theControls_[0]->nsubpixels_)
	{	
		double bestfidel=0;
		double rwa_dt = tgate_/double(rwa_num_time);
		SetNumTimes(rwa_num_time);
		
		for(size_t tries=0; tries<1; tries++)
		{	
			
		//	if(tries==62) iscrap=true;
			std::cout << tgate_ << " " << tries<< "\n";
			
			//allocate controls (..dynamically)
			theControls_[0]->ShiftedGaussian( 0, tgate_, rab, tgate_/2);
			//theControls_[0]->SquarePulse( 0, tgate_, M_PI/tgate_);
			//theControls_[0]->writefile("one0.dat");
			double norm = theControls_[0]->Normalize(rab);
			//if(tries>0) theControls_[0]->SquarePulse( 0, tgate_, M_PI/tgate_);
//				theControls_[1]->GaussianDerivative(0, tgate_, -rab/Delta*norm, tgate_/2);
			//
			//if(tries==2)
			for(size_t j = 0; j < num_time_; ++j)
			{	double x = theControls_[0]->u_[j];
//				//1ch
//				theControls_[2]->u_[j]=(lambda*lambda-4)*(x*x)/Delta/4;
//								-((27 -16 *lambda*lambda + 2*lambda*lambda*lambda*lambda) *x*x*x*x )/(Delta*Delta*Delta)/16;
//				theControls_[3]->u_[j]=lambda*theControls_[1]->u_[j]*x/2/Delta;
//				theControls_[0]->u_[j]+=(lambda*lambda-4)*(x*x*x)/Delta/Delta/8;
//						-(x*x*x*x*x* (104 - 84 *lambda*lambda + 13*lambda*lambda*lambda*lambda))/(8*Delta*Delta*Delta*Delta)/16;
					
//				theControls_[1]->u_[j]*=(1+3*x*x*(65-28*lambda*lambda)/(12*Delta*Delta*Delta)/4);
				//2ch  Adi
//				theControls_[2]->u_[j]=(lambda*lambda/Delta-lamb2*lamb2/Delta2)*(x*x)/4;

//				theControls_[2]->u_[j]-=(lambda*lambda-0.5)*(lambda*lambda+1.5)*(x*x*x*x)/anharmonicity_/anharmonicity_/anharmonicity_/8;
//				theControls_[2]->u_[j]-=(lambda*lambda-0.5)*theControls_[1]->u_[j]*theControls_[1]->u_[j]/anharmonicity_/anharmonicity_/anharmonicity_;
//				theControls_[0]->u_[j]+=(lambda*lambda/Delta/Delta+lamb2*lamb2/Delta2/Delta2)*(x*x*x)/8;

//				theControls_[0]->u_[j]-=x*(0.0625*x*x*x*x*(0.625/4 - 0.5*lambda*lambda + 0.625*lambda*lambda*lambda*lambda + 0.5* (-0.5 - 1.25*lambda*lambda)) - 0.5* (0.5 + lambda*lambda)* (theControls_[1]->u_[j]*theControls_[1]->u_[j])/4)/anharmonicity_/anharmonicity_/anharmonicity_/anharmonicity_;
//				theControls_[1]->u_[j]=0;
				
				//2Ch Drg
//				theControls_[2]->u_[j]=(x*x)*( -lamb2*lamb2/Delta2 +(lambda*lambda - 4)/Delta)/4;
//				theControls_[0]->u_[j]-=(x*x*x)*( (4-lambda*lambda)/Delta/Delta - lamb2*lamb2/Delta2/Delta2 )/8;				
			}
			
			
			//else 
			if(tries>1)
			{
				theControls_[0]->RandomControl(0, M_PI/tgate_*2.0);
				
				cout << theControls_[0]->Normalize(rab)  << " norm\n";
				
				if(tries%2)
						theControls_[1]->RandomControl(-M_PI/anharmonicity_/tgate_/3, M_PI/anharmonicity_/tgate_/3);
			}
			
			UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2); 
			cout << tgate_ <<"\t" << 1-top_fidelity_ << "\n";
			if(top_fidelity_>bestfidel) bestfidel = top_fidelity_;
			//if(top_fidelity_>fidelity_) break;
		}
		
		//cout << 1-abs(rho_[num_time_-1](1,2)) << endl;;
		dataout <<  tgate_ <<"\t" << 1-bestfidel << "\n";
		cout << "best for " << tgate_ <<" was " << 1-bestfidel << "\n";
		//if (top_fidelity_==fidelity_)
		theControls_[0]->writefile("one0.dat");
		theControls_[1]->writefile("one1.dat");
		theControls_[2]->writefile("one2.dat");
	}
	dataout.close();
}

void Grape::sweepenergies(const double min, const double max, const matrix<std::complex<double> >& Shift, const char* outfile)
{  //datafile
	
	ofstream dataout;
	UFs::OpenFile(outfile,dataout, 16);
	
	size_t rwa_num_time=theControls_[0]->npixels_;
	matrix<std::complex<double> > Drift_bckup (H_drift_);
	
	double lambda = sqrt(3.0/2);
	cout <<lambda*lambda << endl;
	double rab = M_PI;
	double lamb2 =  1/sqrt(2);
	double Delta = anharmonicity_;
	double Delta2 = anharmonicity_;
	//H_drift_(0,0)=Delta2-Delta;	
	
	cout << theControls_[0]->Hcontrol_ << endl;
	for(double shift=min; shift<max; shift +=  (max-min)/40 )
	{	
		cout << "Hello\n";
		H_drift_ = Drift_bckup + shift*Shift;
		if(abs(anharmonicity_)>abs(shift)) {
			Delta = shift;
			Delta2 = anharmonicity_;
		} else {
			Delta = anharmonicity_;
			Delta2 = shift;
		}
		double bestfidel=0;
		//double rwa_dt = tgate_/double(rwa_num_time);
		//SetNumTimes(rwa_num_time);
		
		for(size_t tries=0; tries<1; tries++)
		{	
		//	if(tries==62) iscrap=true;
			std::cout << tgate_ << " " << tries<< "\n";
			
			//allocate controls (..dynamically)
			theControls_[0]->ShiftedGaussian( 0, tgate_, rab, tgate_/2);
			//theControls_[0]->SquarePulse( 0, tgate_, M_PI/tgate_);
			//theControls_[0]->writefile("one0.dat");
			double norm = theControls_[0]->Normalize(rab);
			//if(tries>0) theControls_[0]->SquarePulse( 0, tgate_, M_PI/tgate_);
//			theControls_[1]->GaussianDerivative(0, tgate_, -rab/Delta*norm, tgate_/2);
			//
			//if(tries==2)
			for(size_t j = 0; j < num_time_; ++j)
			{	double x = theControls_[0]->u_[j];
//				//1ch
	//			theControls_[2]->u_[j]=(lambda*lambda-4)*(x*x)/Delta/4;
		//						-((27 -16 *lambda*lambda + 2*lambda*lambda*lambda*lambda) *x*x*x*x )/(Delta*Delta*Delta)/16;
		//		theControls_[3]->u_[j]=lambda*theControls_[1]->u_[j]*x/2/Delta;
	//			theControls_[0]->u_[j]+=(lambda*lambda-4)*(x*x*x)/Delta/Delta/8;
		//				-(x*x*x*x*x* (104 - 84 *lambda*lambda + 13*lambda*lambda*lambda*lambda))/(8*Delta*Delta*Delta*Delta)/16;
					
		//		theControls_[1]->u_[j]*=(1+3*x*x*(65-28*lambda*lambda)/(12*Delta*Delta*Delta)/4);
				//2ch  Adi
				theControls_[2]->u_[j]=(lambda*lambda/Delta-lamb2*lamb2/Delta2)*(x*x)/4;

//				theControls_[2]->u_[j]-=(lambda*lambda-0.5)*(lambda*lambda+1.5)*(x*x*x*x)/anharmonicity_/anharmonicity_/anharmonicity_/8;
//				theControls_[2]->u_[j]-=(lambda*lambda-0.5)*theControls_[1]->u_[j]*theControls_[1]->u_[j]/anharmonicity_/anharmonicity_/anharmonicity_;
				theControls_[0]->u_[j]+=(lambda*lambda/Delta/Delta+lamb2*lamb2/Delta2/Delta2)*(x*x*x)/8;

//				theControls_[0]->u_[j]-=x*(0.0625*x*x*x*x*(0.625/4 - 0.5*lambda*lambda + 0.625*lambda*lambda*lambda*lambda + 0.5* (-0.5 - 1.25*lambda*lambda)) - 0.5* (0.5 + lambda*lambda)* (theControls_[1]->u_[j]*theControls_[1]->u_[j])/4)/anharmonicity_/anharmonicity_/anharmonicity_/anharmonicity_;
//				theControls_[1]->u_[j]=0;
				
				//2Ch Drg
//				theControls_[2]->u_[j]=(x*x)*( -lamb2*lamb2/Delta2 +(lambda*lambda - 4)/Delta)/4;
//				theControls_[0]->u_[j]-=(x*x*x)*( (4-lambda*lambda)/Delta/Delta - lamb2*lamb2/Delta2/Delta2 )/8;				
			}
			
			
			//else 
			if(tries>1)
			{
				theControls_[0]->RandomControl(0, M_PI/tgate_*2.0);
				
				cout << theControls_[0]->Normalize(rab)  << " norm\n";
				
				if(tries%2)
						theControls_[1]->RandomControl(-M_PI/anharmonicity_/tgate_/3, M_PI/anharmonicity_/tgate_/3);
			}
			
			UnitaryTransfer(&Grape::Phi4Sub2, &Grape::GradPhi4Sub2); 
			cout << tgate_ <<"\t" << 1-top_fidelity_ << "\n";
			if(top_fidelity_>bestfidel) bestfidel = top_fidelity_;
			//if(top_fidelity_>fidelity_) break;
		}
		
		//cout << 1-abs(rho_[num_time_-1](1,2)) << endl;;
		dataout <<  shift <<"\t" << 1-bestfidel << "\n";
		cout << "best for " << shift <<" was " << 1-bestfidel << "\n";
		//if (top_fidelity_==fidelity_)
		theControls_[0]->writefile("one0.dat");
		theControls_[1]->writefile("one1.dat");
		theControls_[2]->writefile("one2.dat");
	}
	dataout.close();
}
		

inline double Grape::Phi3() const{
	//the measure implemented, phi (PHI_3) is phi = trace[U_desired* UN-1....U_0]/D
	std::complex<double> temp1=0.0;
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += std::conj(rho_desired_(p,q))*rho_[num_time_-1](p,q);
		}
	}
	return std::real(temp1)*one_on_dim_;
}
inline double Grape::GradPhi3(const size_t j, const size_t k) const{
	//phi_3
	std::complex<double> temp1=0;
	matrix<std::complex<double> > Rtemp;
	Rtemp=H_controls_tdep_[k]*rho_[j];
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 +=std::conj(lambda_[j](p,q))*Rtemp(p,q);	
		}
	}	
	return epsilon_*std::imag(temp1);
}
inline double Grape::Phi4() const{
	//the measure implemented, phi (PHI_4) is phi = |trace[U_desired* UN-1....U_0 ]|^2/D^2
	std::complex<double> temp1=0.0;
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += std::conj(rho_desired_(p,q))*rho_[num_time_-1](p,q);
		}
	}
	return real(temp1*std::conj(temp1))*one_on_dim_*one_on_dim_;
}
inline double Grape::GradPhi4(const size_t j, const size_t k) const{
	//phi_4
	//2.0*h_*imag(QOs::Expectation( MOs::Dagger(L),H*R)*QOs::Expectation(MOs::Dagger(R),L));
	std::complex<double> temp1=0, temp2=0;
	matrix<std::complex<double> > Rtemp;
	Rtemp=H_controls_tdep_[k]*rho_[j];
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 +=std::conj(lambda_[j](p,q))*Rtemp(p,q);
			temp2 += std::conj(rho_[j](p,q))*lambda_[j](p,q);	
		}
	}
	return 2.0*epsilon_*std::imag(temp1*temp2);//*one_on_dim_*one_on_dim_;
}
inline double Grape::Phi4Sub2() const{
	//the measure implemented, phi (PHI_4_sub) is phi = |sum_{i=0,1} <i|U_desired* UN-1....U_0|i>|^2/2^2
	std::complex<double> temp1=0.0;
	for(size_t q=1; q< 3; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += std::conj(rho_desired_(p,q))*rho_[num_time_-1](p,q);		}
	}
	return std::real(temp1*std::conj(temp1))*0.25;
}
inline double Grape::GradPhi4Sub2(const size_t j, const size_t k){
	//Phi_4_subsystem
	//2.0*h_*imag(QOs::Expectation( MOs::Dagger(L),H*R)*QOs::Expectation(MOs::Dagger(R),L));
	std::complex<double> temp1=0, temp2=0;
	//matrix<std::complex<double> > Rtemp(dim_, dim_);
	DFtemp_=H_controls_tdep_[k]*rho_[j];
	for(size_t q=1; q< 3; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 +=std::conj(lambda_[j](p,q))*DFtemp_(p,q);
			temp2 += std::conj(rho_[j](p,q))*lambda_[j](p,q);	
		}
	}
	return 2.0*epsilon_*std::imag(temp1*temp2);//*one_on_dim_*one_on_dim_;
}
inline double Grape::Phi4TraceCav() const{
	//the measure implemented, phi (PHI_4_sub) is phi = |sum_{i=0,1} <i|U_desired* UN-1....U_0|i>|^2/2^2
	std::complex<double> temp1=0.0;
	matrix<std::complex<double> > Rtemp=MOs::TraceOutB(rho_[num_time_-1],5);///5=dimCav
	for(size_t q=0; q< 2; ++q){
		for(size_t p=0; p<dimQ_; ++p){
			temp1 += std::conj(rho_desired_(p,q))*Rtemp(p,q);
		}
	}
	return std::real(temp1*std::conj(temp1))*0.25;
}
inline double Grape::GradPhi4TraceCav(const size_t j, const size_t k) const{
	//Phi_4_subsystem
	//2.0*h_*imag(QOs::Expectation( MOs::Dagger(L),H*R)*QOs::Expectation(MOs::Dagger(R),L));
	std::complex<double> temp1=0, temp2=0;
	matrix<std::complex<double> > Rtemp;
	Rtemp=MOs::TraceOutB(H_controls_tdep_[k]*rho_[j],5);///5=dimCav
	for(size_t q=0; q< 2; ++q){
		for(size_t p=0; p<dimQ_; ++p){
			temp1 +=std::conj(lambda_[j](p,q))*Rtemp(p,q);
			temp2 += std::conj(rho_[j](p,q))*lambda_[j](p,q);	
		}
	}
	return 2.0*epsilon_*std::imag(temp1*temp2);//*one_on_dim_*one_on_dim_;
}
inline double Grape::GradPhi4Sub2Filter(const size_t j, const size_t k) const{
	double gradient=0.0;
	if(j<nsubpixels_[k] || j>num_time_-nsubpixels_[k]) return gradient;
	std::complex<double> temp1=0, temp2=0;
	matrix<std::complex<double> > Rtemp;
	int endtime = j + 100;
	if(j>num_time_-100) endtime = num_time_;
	for(size_t l = (j>100)*(j-100); l <endtime; ++l){
	  Rtemp=cos(((double)l)*drive_freq_*h_ - (k%2)*M_PI/2)*H_controls_tdep_[k]*rho_[l]*exp(-abs((int)j-(int)l)*h_/filter_rate_)*(h_/filter_rate_/2.0 /filter_norm_);
	  for(size_t q=0; q< 2; ++q){
	  	 for(size_t p=0; p<dim_; ++p){
			temp1 +=std::conj(lambda_[l](p,q))*Rtemp(p,q);
			temp2 += std::conj(rho_[l](p,q))*lambda_[l](p,q);	
		 }
	  }
	  gradient +=2.0*epsilon_*std::imag(temp1*temp2);
	}
	return gradient;
}

inline double Grape::GradPhi4Sub2RWFilter(const size_t j, const size_t k) const{
	double gradient=0.0;
	double prefac;
	if(j<nsubpixels_[k]*12 || j>num_time_-12*nsubpixels_[k]) return gradient;
	int span=100;
	int endtime = j + span;
	if(j>num_time_-span) endtime = num_time_;
	
	double phase = config_*2*M_PI/nconfigs_; 
	
	for(size_t l = (j>span)*(j-span); l <endtime; ++l) 	
	//int l=j;
	{
		prefac = cos(phase + ((double)j)*(drive_freq_)*h_ - (k%2)*M_PI/2)*exp(-abs((int)j-(int)l)*h_/filter_rate_)*(h_/filter_rate_/2.0 /filter_norm_); 
		Udiag_[k][0][0]=Udiag_[k][1][0]=1;
		Udiag_[k][0][1]=exp(-std::complex<double>(0.0,phase +(double)j*drive_freq_*h_));
		Udiag_[k][0][2]=exp(-std::complex<double>(0.0,2.0*(phase + (double)j*drive_freq_*h_)));
		Udiag_[k][0][3]=exp(-std::complex<double>(0.0,3.0*(phase + (double)j*drive_freq_*h_)));
		Udiag_[k][1][1]=exp(std::complex<double>(0.0,(phase + (double)j*drive_freq_*h_)));
		Udiag_[k][1][2]=exp(std::complex<double>(0.0,2.0*(phase + (double)j*drive_freq_*h_)));	
		Udiag_[k][1][3]=exp(std::complex<double>(0.0,3.0*(phase + (double)j*drive_freq_*h_)));	
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				H_controls_tdep_[k](p,q)= Udiag_[k][0][q]*H_controls_[k](p,q)*Udiag_[k][1][p]*prefac;
			
			}
		}		std::complex<double> temp1=0, temp2=0;
	matrix<std::complex<double> > Rtemp;
	Rtemp=H_controls_tdep_[k]*rho_[j];
	for(size_t q=0; q< 2; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 +=std::conj(lambda_[j](p,q))*Rtemp(p,q);
			temp2 += std::conj(rho_[j](p,q))*lambda_[j](p,q);	
		}
	}
	
	gradient+= 2.0*epsilon_*std::imag(temp1*temp2);//*one_on_dim_*one_on_dim_;
	}
	return gradient;
}

void Grape::SetLinear(size_t j)
{
/*	std::complex<double> lambda = H_controls_[0](1,2);
	std::complex<double> lambdai = H_controls_[1](1,2);
	if(config_==1)
	{
		H_controls_[0](1,2)=H_controls_[0](2,1)=lambda*0.93;
		H_controls_[1](1,2)=lambdai*0.93;
	} else if(config_==2){
		H_controls_[0](1,2)=H_controls_[0](2,1)=lambda*1.07;
		H_controls_[1](1,2)=lambdai*1.07;
	}
	H_controls_[1](2,1)=-H_controls_[1](1,2);
*/	
	for(size_t k = 0; k < num_controls_; ++k){	
		H_controls_tdep_[k]=controls_[k][j]*H_controls_[k];
		
	}
//	H_controls_[0](1,2)=H_controls_[0](2,1)=lambda;
//	H_controls_[1](1,2)=lambdai;
//	H_controls_[1](2,1)=-H_controls_[1](1,2);

}
void Grape::SetGradLinear(size_t j)
{
/*	std::complex<double> lambda = H_controls_[0](1,2);
	std::complex<double> lambdai = H_controls_[1](1,2);
	if(config_==1)
	{
		H_controls_[0](1,2)=H_controls_[0](2,1)=lambda*0.93;
		H_controls_[1](1,2)=lambdai*0.93;
	} else if(config_==2){
		H_controls_[0](1,2)=H_controls_[0](2,1)=lambda*1.07;
		H_controls_[1](1,2)=lambdai*1.07;
	}
	H_controls_[1](2,1)=-H_controls_[1](1,2);
*/	
	for(size_t k = 0; k < num_controls_; ++k){	
		H_controls_tdep_[k]=H_controls_[k];
	}
//	H_controls_[0](1,2)=H_controls_[0](2,1)=lambda;
//	H_controls_[1](1,2)=lambdai;
//	H_controls_[1](2,1)=-H_controls_[1](1,2);
}

void Grape::SetLab2Quadratures(size_t j)
{
	for(size_t k = 0; k < num_controls_; ++k){	
		if((k+1)%3)
			H_controls_tdep_[k]= controls_[k][j]*H_controls_[k]*cos(((double)j)*(drive_freq_)*h_ - (k%2)*M_PI/2);
		else  H_controls_tdep_[k]=controls_[k][j]*H_controls_[k];	
	}
}
void Grape::SetGradLab2Quadratures(size_t j)
{
	for(size_t k = 0; k < num_controls_; ++k){	
		if((k+1)%3)
			H_controls_tdep_[k] = H_controls_[k]*cos(((double)j)*(drive_freq_/*+controls_[(k/3)*3+2][j]*/)*h_ - (k%2)*M_PI/2);
		else  H_controls_tdep_[k] = H_controls_[k];
	}	
}

void Grape::SetLab2QuadFiltered(size_t j)
{
	double cntrl;
	//std::cout<< j << " " << controls_[0][j] <<"\n";
	for(size_t k = 0; k < num_controls_; ++k) {
		cntrl=0;
		if((k+1)%3)
		{
		for(size_t l = 0; l < num_time_; ++l)					
		{														
			cntrl += controls_[k][l]*exp(-abs((int)j-(int)l)*h_/filter_rate_)*h_/filter_rate_/2 /filter_norm_;	
		}			
		H_controls_tdep_[k]= cntrl*cos(((double)j)*(drive_freq_/*+controls_[(k/3)*3+2][j]*/)*h_ - (k%2)*M_PI/2)*H_controls_[k];
		}
		else  H_controls_tdep_[k]=0.0*H_controls_[k];
		//std::cout<< controls_[k][j] << " " << cntrl << "\n";
		//if(!k) std::cout<< j << " " << cntrl << "\n";
	}	
}

void Grape::SetGradLab2QuadFiltered(size_t j)
{
		for(size_t k = 0; k < num_controls_; ++k) {
			H_controls_tdep_[k]= H_controls_[k];
		}
}

void Grape::SetCounterRotating(size_t j)
{
	double phase = config_*2*M_PI/nconfigs_; 
		
	for(size_t k = 0; k < num_controls_; ++k) {	
		for(size_t l=0; l<dim_; l++) {
			Udiag_[k][0][l]=exp(-std::complex<double>(0.0,l*(phase +(double)j*drive_freq_*h_)));
			Udiag_[k][1][l]=exp(std::complex<double>(0.0,l*(phase + (double)j*drive_freq_*h_)));
		}		
		if((k+1)%3)
		{	
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				H_controls_tdep_[k](p,q)= Udiag_[k][0][q]*controls_[k][j]*H_controls_[k](p,q)*cos(phase + ((double)j)*(drive_freq_)*h_ - (k%2)*M_PI/2)*Udiag_[k][1][p];
			}
		}	
		}
		else  
		H_controls_tdep_[k]=(controls_[k][j])*H_controls_[k];
	}
	
}
void Grape::SetGradCounterRotating(size_t j)
{
	double phase = config_*2*M_PI/nconfigs_; 
	
	int l=j;
	
	//for(double phase=0; phase<2*M_PI; phase+=2*M_PI/7) 
	{
		
	for(size_t k = 0; k < num_controls_; ++k) {	
		for(size_t l=0; l<dim_; l++) {
			Udiag_[k][0][l]=exp(-std::complex<double>(0.0,l*(phase +(double)j*drive_freq_*h_)));
			Udiag_[k][1][l]=exp(std::complex<double>(0.0,l*(phase + (double)j*drive_freq_*h_)));
		}
		if((k+1)%3)
		{	
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				H_controls_tdep_[k](p,q)= Udiag_[k][0][q]*H_controls_[k](p,q)*cos(phase + ((double)j)*(drive_freq_)*h_ - (k%2)*M_PI/2)*Udiag_[k][1][p];
				//H_controls_tdep_[k](p,q)= Udiag_[k][0][q]*H_controls_[k](p,q)*cos(phase + ((double)j)*(drive_freq_)*h_ - (k%2)*M_PI/2)*Udiag_[k][1][p]*exp(-abs((int)j-(int)l)*h_/filter_rate_)*(h_/filter_rate_/2.0 /filter_norm_);
			
			}
		}	
		}
		else  H_controls_tdep_[k]=H_controls_[k];
	/*	{
			for(size_t q=0; q< dim_; ++q){
			  for(size_t p=0; p<dim_; ++p){
					H_controls_tdep_[k](p,q)= Udiag_[k][0][q]*( -controls_[(k/3)*3][j]*H_controls_[(k/3)*3](p,q)*sin( ((double)j)*(qubit_freq_+controls_[k][j])*h_ ) 
						+ controls_[(k/3)*3+1][j]*H_controls_[(k/3)*3+1](p,q)*cos( ((double)j)*(qubit_freq_+controls_[k][j])*h_)   )*Udiag_[k][1][p]*((double)j)*h_;
			  }
			}	
		
		}*/
	}
	}
	
}
void Grape::SetCounterRotatingFilt(size_t j)
{

	double phase = config_*2*M_PI/nconfigs_; 
	
	double cntrl;
	for(size_t k = 0; k < num_controls_; ++k) {	
		Udiag_[k][0][0]=Udiag_[k][1][0]=1;
		Udiag_[k][0][1]=exp(-std::complex<double>(0.0,phase +(double)j*drive_freq_*h_));
		Udiag_[k][0][2]=exp(-std::complex<double>(0.0,2.0*(phase + (double)j*drive_freq_*h_)));
		Udiag_[k][0][3]=exp(-std::complex<double>(0.0,3.0*(phase + (double)j*drive_freq_*h_)));
		Udiag_[k][1][1]=exp(std::complex<double>(0.0,(phase + (double)j*drive_freq_*h_)));
		Udiag_[k][1][2]=exp(std::complex<double>(0.0,2.0*(phase + (double)j*drive_freq_*h_)));	
		Udiag_[k][1][3]=exp(std::complex<double>(0.0,3.0*(phase + (double)j*drive_freq_*h_)));		
		cntrl=0;
		if((k+1)%3)
		{	
		int endtime = j + 300;
		if(j>num_time_-300) endtime = num_time_;
		for(size_t l = (j>300)*(j-300); l <endtime; ++l)
		//size_t l=j;
		//for(size_t l =0; l <num_time_; ++l)
		{														
			cntrl += controls_[k][l]*exp(-abs((int)j-(int)l)*h_/filter_rate_)*h_/filter_rate_/2 /filter_norm_;	
		}
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				H_controls_tdep_[k](p,q)= Udiag_[k][0][q]*cntrl*H_controls_[k](p,q)*cos(phase + ((double)j)*(drive_freq_)*h_ - (k%2)*M_PI/2)*Udiag_[k][1][p];
			}
		}	
		}
		else  H_controls_tdep_[k]=controls_[k][j]*H_controls_[k];
	}
}
void Grape::SetGradCounterRotatingFilt(size_t j)
{
		for(size_t k = 0; k < num_controls_; ++k) {
			H_controls_tdep_[k]= H_controls_[k];
		}
}



#endif /* Grape_h */

