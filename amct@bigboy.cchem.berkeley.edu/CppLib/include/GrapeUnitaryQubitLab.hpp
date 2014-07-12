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
#ifndef GrapeUnitaryQubitLab_h
#define GrapeUnitaryQubitLab_h
#include "MatrixExponential.hpp"
#include "QuantumOperations.hpp"
#include "GrapeUnitaryQubit.hpp"

class GrapeUnitaryQubitLab : public GrapeUnitaryQubit{
	public:
		GrapeUnitaryQubitLab(size_t dim, size_t num_controls, size_t num_time);// number of controls, the number of time points
		virtual ~GrapeUnitaryQubitLab();
	
		void SetPhysicalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, double qubit_freq_, double anharmonicity_);
				
		void SetLab2Quadratures(size_t j);
		void SetGradLab2Quadratures(size_t j);
		void SetCounterRotating(size_t j);
		void SetGradCounterRotating(size_t j);
		
		
	protected:
		
		
};

inline GrapeUnitaryQubitLab::GrapeUnitaryQubitLab(size_t dim, size_t num_controls, size_t num_time) : GrapeUnitaryQubit(dim, num_controls, num_time){
	
	config_=0;
	drive_freq_=0;
	
}

GrapeUnitaryQubitLab::~GrapeUnitaryQubitLab(){

}

inline void GrapeUnitaryQubitLab::SetPhysicalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, 
											double qubit_freq, double anharmonicity)	
{	qubit_freq_=qubit_freq;
	anharmonicity_=anharmonicity;	
	nconfigs_=1;
	SetNumericalParameters(fidelity, h, base_a, epsilon, tolerance, max_iter);
}

void GrapeUnitaryQubitLab::SetLab2Quadratures(size_t j)
{
	for(size_t k = 0; k < num_controls_; ++k){	
		if((k+1)%3)
			H_controls_tdep_[k]= controls_[k][j]*(*H_controls_[k])*cos(((double)j)*(drive_freq_)*h_ - (k%2)*M_PI/2);
		else  H_controls_tdep_[k]=controls_[k][j]*(*H_controls_[k]);	
	}
}
void GrapeUnitaryQubitLab::SetGradLab2Quadratures(size_t j)
{
	for(size_t k = 0; k < num_controls_; ++k){	
		if((k+1)%3)
			H_controls_tdep_[k] = (*H_controls_[k])*cos(((double)j)*(drive_freq_/*+controls_[(k/3)*3+2][j]*/)*h_ - (k%2)*M_PI/2);
		else  H_controls_tdep_[k] = *H_controls_[k];
	}	
}

void GrapeUnitaryQubitLab::SetCounterRotating(size_t j)
{
	double phase = config_*2*M_PI/nconfigs_; 
				
	for(size_t k = 0; k < num_controls_; ++k) {
		for(size_t l=0; l<dim_; l++) {
			Udiag_[k][0][l]=exp(-std::complex<double>(0.0,l*(phase +(double)j*drive_freq_*h_)));
			Udiag_[k][1][l]=exp(std::complex<double>(0.0,l*(phase + (double)j*drive_freq_*h_)));
		}		
		if((k+1)%3)
		{	
//		cout << k << " " << j << " " << cos(phase + ((double)j)*(drive_freq_)*h_ - (k%2)*M_PI/2) << 	" "<< Udiag_[k][0][1] << " " << Udiag_[k][1][1] << endl;
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				H_controls_tdep_[k](p,q)= Udiag_[k][0][q]*controls_[k][j]*(*H_controls_[k])(p,q)*cos(phase + ((double)j)*(drive_freq_)*h_ - (k%2)*M_PI/2)*Udiag_[k][1][p];
//				cout << k << " " << j << " " << p+2*q << " " << (*H_controls_[k])(p,q) << endl;
			}
		}	
		}
		else  
		H_controls_tdep_[k]=(controls_[k][j])*(*H_controls_[k]);
	}
	
}
void GrapeUnitaryQubitLab::SetGradCounterRotating(size_t j)
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
				H_controls_tdep_[k](p,q)= Udiag_[k][0][q]*(*H_controls_[k])(p,q)*cos(phase + ((double)j)*(drive_freq_)*h_ - (k%2)*M_PI/2)*Udiag_[k][1][p];
				//H_controls_tdep_[k](p,q)= Udiag_[k][0][q]*(*H_controls_[k])(p,q)*cos(phase + ((double)j)*(drive_freq_)*h_ - (k%2)*M_PI/2)*Udiag_[k][1][p]*exp(-abs((int)j-(int)l)*h_/filter_rate_)*(h_/filter_rate_/2.0 /filter_norm_);
			
			}
		}	
		}
		else  H_controls_tdep_[k]=*H_controls_[k];
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

#endif /* GrapeUnitaryQubitLab_h */

