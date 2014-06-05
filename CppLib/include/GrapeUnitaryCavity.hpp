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
#ifndef GrapeUnitaryCavity_h
#define GrapeUnitaryCavity_h
#include "MatrixExponential.hpp"
#include "QuantumOperations.hpp"
#include "GrapeUnitary.hpp"

class GrapeUnitaryCavity : public GrapeUnitary{
	public:
		GrapeUnitaryCavity(size_t dim, size_t num_controls, size_t num_time);// number of controls, the number of time points
		virtual ~GrapeUnitaryCavity();
	
		void SetPhysicalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, double qubit_freq_, double anharmonicity_);
		
		double Phi4TraceCav() const;
		double GradPhi4TraceCav(const size_t j, const size_t k) ;
		double Phi43levCav() const;
		double GradPhi43levCav(const size_t j, const size_t k) ;		
					
};

inline GrapeUnitaryCavity::GrapeUnitaryCavity(size_t dim, size_t num_controls, size_t num_time) : GrapeUnitary(dim, num_controls, num_time){
	
	config_=0;
	drive_freq_=0;
	
}

GrapeUnitaryCavity::~GrapeUnitaryCavity(){

}

inline void GrapeUnitaryCavity::SetPhysicalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, 
											double qubit_freq, double anharmonicity)	
{	qubit_freq_=qubit_freq;
	anharmonicity_=anharmonicity;	
	nconfigs_=1;
	SetNumericalParameters(fidelity, h, base_a, epsilon, tolerance, max_iter);
}

inline double GrapeUnitaryCavity::Phi4TraceCav() const{
	//the measure implemented, phi (PHI_4_sub) is phi = |sum_{i=0,1} <i|U_desired* UN-1....U_0|i>|^2/2^2
	std::complex<double> temp1=0.0;
//	matrix<std::complex<double> > Rtemp=MOs::TraceOutB(rho_[num_time_-1],5);///5=dimCav
	double topdim=4;
	double dimsq=topdim*topdim;
	for(size_t p=0; p<dim_; ++p){
		for(size_t q=0; q< topdim; ++q)
			temp1 += std::conj(rho_desired_(p,q))*rho_[num_time_-1](p,q);
//		temp1 += std::conj(rho_desired_(p,3))*rho_[num_time_-1](p,3);
	}
	return std::real(temp1*std::conj(temp1))/dimsq;
}
inline double GrapeUnitaryCavity::GradPhi4TraceCav(const size_t j, const size_t k) {
	//Phi_4_subsystem
	//2.0*h_*imag(QOs::Expectation( MOs::Dagger(L),H*R)*QOs::Expectation(MOs::Dagger(R),L));
	std::complex<double> temp1=0, temp2=0;
	matrix<std::complex<double> > Rtemp;
	Rtemp=H_controls_tdep_[k]*rho_[j];
	
	double topdim=4;
	for(size_t p=0; p<dim_; ++p){
		for(size_t q=0; q< topdim; ++q){
			temp1 +=std::conj(lambda_[j](p,q))*Rtemp(p,q);
			temp2 += std::conj(rho_[j](p,q))*lambda_[j](p,q);	
		}
	//	temp1 +=std::conj(lambda_[j](p,3))*Rtemp(p,3);
	//	temp2 += std::conj(rho_[j](p,3))*lambda_[j](p,3);
	}
	return 2.0*epsilon_*std::imag(temp1*temp2);//*one_on_dim_*one_on_dim_;
}
inline double GrapeUnitaryCavity::Phi43levCav() const{
	//the measure implemented, phi (PHI_4_sub) is phi = |sum_{i=0,1} <i|U_desired* UN-1....U_0|i>|^2/2^2
	std::complex<double> temp1=0.0;
//	matrix<std::complex<double> > Rtemp=MOs::TraceOutB(rho_[num_time_-1],5);///5=dimCav
	double topdim=5;
	int dimQ = 2;
	double dimsq=topdim*topdim;
	for(size_t p=0; p<dim_; ++p){
		//for(size_t q=0; q< 2; ++q)
		temp1 += std::conj(rho_desired_(p,0))*rho_[num_time_-1](p,0);
	//	temp1 += std::conj(rho_desired_(p,1))*rho_[num_time_-1](p,1);
		temp1 += std::conj(rho_desired_(p,dimQ+1))*rho_[num_time_-1](p,dimQ+1);
	}
	return std::real(temp1*std::conj(temp1))/4;
}
inline double GrapeUnitaryCavity::GradPhi43levCav(const size_t j, const size_t k) {
	//Phi_4_subsystem
	//2.0*h_*imag(QOs::Expectation( MOs::Dagger(L),H*R)*QOs::Expectation(MOs::Dagger(R),L));
	std::complex<double> temp1=0, temp2=0;
	matrix<std::complex<double> > Rtemp;
	Rtemp=H_controls_tdep_[k]*rho_[j];
	
	double topdim=5;
	int dimQ = 2;
	for(size_t p=0; p<dim_; ++p){
		//for(size_t q=0; q< 2; ++q){
			temp1 +=std::conj(lambda_[j](p,0))*Rtemp(p,0);
			temp2 += std::conj(rho_[j](p,0))*lambda_[j](p,0);	
		//}
		
	//		temp1 +=std::conj(lambda_[j](p,1))*Rtemp(p,1);
	//		temp2 += std::conj(rho_[j](p,1))*lambda_[j](p,1);
		temp1 +=std::conj(lambda_[j](p,dimQ+1))*Rtemp(p,dimQ+1);
		temp2 += std::conj(rho_[j](p,dimQ+1))*lambda_[j](p,dimQ+1);
	}
	return 2.0*epsilon_*std::imag(temp1*temp2);//*one_on_dim_*one_on_dim_;
}

/*
inline double GrapeUnitaryCavity::Phi4TraceCav() const{
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
inline double GrapeUnitaryCavity::GradPhi4TraceCav(const size_t j, const size_t k) {
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
}*/


#endif /* GrapeUnitaryCavity_h */

