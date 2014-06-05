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
#ifndef GrapeTransfer_h
#define GrapeTransfer_h
#include "MatrixExponential.hpp"
#include "QuantumOperations.hpp"
#include "ControlOptimizer.hpp"

class GrapeTransfer : public ControlOptimizer{
	public:
		Grape(size_t dim, size_t num_controls, size_t num_time);// number of controls, the number of time points
		virtual ~Grape();
	
		void SetPhysicalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, double qubit_freq_, double anharmonicity_);
		
		void StateTransfer(double (Grape::* phi)() const, double (Grape::* gradphi)(const size_t j, const size_t k) const, 
				void (Grape::* ptrSetHnonlinearity)(size_t j), void (Grape::* ptrSetHgradnonlin)(const size_t j));
		void StateTransfer(double (Grape::* phi)() const, double (Grape::* gradphi)(const size_t j, const size_t k) const)
		{   StateTransfer(phi, gradphi, &Grape::SetLinear, &Grape::SetGradLinear);
		}
		
		double Phi0() const;
		double GradPhi0(const size_t j,const size_t k) const;
		
		void SetLinear(size_t j);
		void SetGradLinear(size_t j);

		
	protected:
		
		//system specific parameters:
		double qubit_freq_;
		double anharmonicity_;
		double drive_freq_;
		double filter_norm_;
		double filter_rate_;
		
};
inline GrapeTransfer::GrapeTransfer(size_t dim, size_t num_controls, size_t num_time) : ControlOptimizer(dim, num_controls, num_time){
	
	config_=0;
	drive_freq_=0;
	
}
Grape::~Grape(){

}

inline void GrapeTransfer::SetPhysicalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, 
											double qubit_freq, double anharmonicity)	
{	qubit_freq_=qubit_freq;
	anharmonicity_=anharmonicity;	
	nconfigs_=1;
	SetNumericalParameters(fidelity, h, base_a, epsilon, tolerance, max_iter);
}

//Performs the state transfer with given Hamiltonian(s)
void GrapeTransfer::StateTransfer(double (Grape::* phi)() const, double (Grape::* gradphi)(const size_t j, const size_t k) const, 
   void (Grape::* ptrSetHnonlinearity)(size_t j) , void (Grape::* ptrSetHgradnonlin)(const size_t j)){
	

	//Some flags to make sure all parameters are set
	size_t test=0;
	for(size_t k = 0; k < num_controls_+3; ++k)
	{
		test+=controlsetflag_[k];
		if(verbose==yes)
		{
		 	std::cout << k << "th control flag is " << controlsetflag_[k] << std::endl;
		}
	}
	if(test != 0)
		UFs::MyError("Grape::StateTransfer(): you have not set the drift and all control hamiltonians:\n");
		
	//Set the counters to zero
	count_ =0;
	pos_count_ =0;
	//the power scale for epsilot (epsilon = base_a^power_)
	power_ =0;
	// fidelity ranges from 0 to 1 with 0 be orthogonal and 1 being the same time
	double current_fidelity=0.0, current_delta_fidelity=1.0, last_fidelity=-1.0;
	
	//Starting the algorithim	
	while(current_delta_fidelity*current_delta_fidelity > tolerance_){

		if(count_> max_iter_){	
			std::cout << "Grape::StateTransfer(), algorithim did not converge in the allow max iterations" << std::endl;
			break;
		}
		count_ = count_ + 1;
		
		// the propagators and the foraward evolution
		for(size_t j = 0; j < num_time_; ++j){
			Htemp_ = H_drift_;
			for(size_t k = 0; k < num_controls_; ++k){
				Htemp_ += controls_[k][j]*H_controls_[k];
			}
			// std::cout << Htemp_ << std::endl;
			Unitaries_[j]=ExpM::EigenMethod(Htemp_,-i*h_);
			if(j==0){
				rho_[0] = Unitaries_[0]*rho_initial_*MOs::Dagger(Unitaries_[0]);
				//rho_[0] = rhoupdate(Unitaries_[0],rho_initial);
			}
			else{
				rho_[j] = Unitaries_[j]*rho_[j-1]*MOs::Dagger(Unitaries_[j]);
			}
		}
		
		//the fidelities and break if we reach the required fidelity
		current_fidelity=(this->*phi)();
		if(current_fidelity >= fidelity_)
			break;
		current_delta_fidelity=current_fidelity-last_fidelity;
		last_fidelity=current_fidelity;
		// std::cout << current_fidelity << '\t' <<  current_delta_fidelity << std::endl;

		//the update in the controls
		if(current_delta_fidelity>=0){
			// the backward evolution
			lambda_[num_time_-1] = rho_desired_;
			for(size_t j = num_time_-1; j > 0; --j){
				//Udg.rho.U
				lambda_[j-1] = MOs::Dagger(Unitaries_[j])*lambda_[j]*Unitaries_[j];
			}
			for(size_t j = 0; j < num_time_; ++j){
				for(size_t k =0; k < num_controls_; ++k){
					// std::cout << power_ << std::endl;
					gradient_[k][j]=(this->*gradphi)(j,k);
					controls_[k][j]+=pow(base_a_,power_)*(gradient_[k][j]-alpha_*controls_[k][j]);
				}
			}
			++pos_count_;//if we get to correct directions we set power to be bigger
			if(pos_count_ ==2){
				power_++;
				pos_count_ =0;
			}
		}
		else{
			power_--;
			pos_count_ =0;
			for(size_t j = 0; j < num_time_; ++j){
				for(size_t k =0; k < num_controls_; ++k){
					controls_[k][j]+=pow(base_a_,power_)*(1-base_a_)*gradient_[k][j];
					// std::cout << power_ << std::endl;
				}
			}
					
		}
		
		
		if (count_%100==0){
			std::cout << "step " << count_ << " which has an fidelity of " << current_fidelity << " and a delta fidelity of " << current_delta_fidelity << std::endl; 
		}
	}
	std::cout << "step " << count_ << " which has an fidelity of " << current_fidelity << " and a delta fidelity of " << current_delta_fidelity << std::endl;  
}




inline double GrapeTransfer::Phi0() const{
	//the measure implemented, phi (PHI_0) is phi = trace[rho_desired* UN-1....U_0 rho_initial U_0\dg ... U_N-1\dg]
	std::complex<double> temp1=0.0;
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += conj(rho_desired_(p,q))*rho_[num_time_-1](p,q);
		}
	}
	return real(temp1);
}
inline double GrapeTransfer::GradPhi0(const size_t j, const size_t k) const{
	//phi_0
	std::complex<double> temp1=0;
	matrix<std::complex<double> > Rtemp;
	Rtemp=QOs::SuperHami(H_controls_tdep_[k],rho_[j]);
	//Fill in the upper triangle
	for (size_t p=0; p< dim_; ++p){
		for(size_t q=p+1; q < dim_; ++q){
			Rtemp(p,q)=std::conj(Rtemp(q,p));
		}
	}	
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 +=std::conj(lambda_[j](p,q))*Rtemp(p,q);	
		}
	}
	// std::cout << h_*QOs::Expectation(L,-i*H*R+i*R*H)- h_*temp1 << std::endl;
	// std::cout << temp1 << std::endl;
	// 	Rtemp.SetOutputStyle(Matrix);
	// std::cout << Rtemp << std::endl;		
	return epsilon_*real(temp1);
}

void GrapeTransfer::SetLinear(size_t j)
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
void GrapeTransfer::SetGradLinear(size_t j)
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

#endif /* GrapeTransfer_h */

