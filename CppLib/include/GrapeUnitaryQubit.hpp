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
#ifndef GrapeUnitaryQubit_h
#define GrapeUnitaryQubit_h
#include "MatrixExponential.hpp"
#include "QuantumOperations.hpp"
#include "OptimizeEvolution.hpp"

class GrapeUnitaryQubit : public OptimizeEvolution{
	
	
	
	public:
		GrapeUnitaryQubit(size_t dim, size_t num_controls, size_t num_time);// number of controls, the number of time points
		virtual ~GrapeUnitaryQubit();
	
		void SetPhysicalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, double qubit_freq_, double anharmonicity_);
		
		inline void setlevel0index(size_t newindex) { assert(newindex<dim_); level0index_=newindex; }
		double Phi4Sub2() const;
		double GradPhi4Sub2(const size_t j, const size_t k);
		
	protected:
		
};

inline GrapeUnitaryQubit::GrapeUnitaryQubit(size_t dim, size_t num_controls, size_t num_time) : OptimizeEvolution(dim, num_controls, num_time){
	
	config_=0;
	drive_freq_=0;
	level0index_=0;
	
}

GrapeUnitaryQubit::~GrapeUnitaryQubit() {

}

inline void GrapeUnitaryQubit::SetPhysicalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, 
											double qubit_freq, double anharmonicity)	
{	qubit_freq_=qubit_freq;
	anharmonicity_=anharmonicity;	
	nconfigs_=1;
	SetNumericalParameters(fidelity, h, base_a, epsilon, tolerance, max_iter);
}





#endif /* GrapeUnitaryQubit_h */

