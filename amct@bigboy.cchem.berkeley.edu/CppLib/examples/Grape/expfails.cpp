/*
Name: GrapeUnitaryEnvelope.cpp
Author: felix motzoi

Dependences: GrapeUnitaryEnvelope.hpp  
Brief Discription: This demonstrates grape for lab frame vs rotating frame with large pixels

Version History
	v0: Using AnalyticControl class, Oct. 6 2009.

*/
#include "MatrixExponential.hpp"

using namespace std;

int main (int argc, char const *argv[]){

	verbose=no;
	cout << "Running program " << argv[0] << endl;
	srand ( time(NULL) );
	
	size_t dimQ=2;
	size_t refdim=dimQ*dimQ*dimQ*dimQ*dimQ;
	
	//Typical operators	
		complex<double> ii(0.0,1.0);
		
		cout << "declaring Hilbert space\n";
		
		matrix<complex<double> > HDriftRef(refdim, refdim);		
		matrix<complex<double> > Hcouple(dimQ*dimQ, dimQ*dimQ);
		
		matrix<complex<double> > IdentQ(dimQ, dimQ), IdentQ2(dimQ*dimQ, dimQ*dimQ), IdentQ3(dimQ*dimQ*dimQ, dimQ*dimQ*dimQ);
		MOs::Identity(IdentQ); IdentQ2=MOs::TensorProduct(IdentQ,IdentQ); IdentQ3=MOs::TensorProduct( IdentQ,IdentQ2);
		matrix<complex<double> > o(dimQ,dimQ), od(dimQ,dimQ);
	
		MOs::Destroy(o); od = MOs::Dagger(o); 
		
		Hcouple =  ( MOs::TensorProduct(o, od) + MOs::TensorProduct(od, o) );
		HDriftRef =  MOs::TensorProduct(IdentQ3, Hcouple)
				+  1.0*(M_PI) * MOs::TensorProduct(IdentQ2, MOs::TensorProduct( Hcouple , IdentQ) )
				+  1.1*(M_PI) *MOs::TensorProduct(IdentQ2, MOs::TensorProduct( o,MOs::TensorProduct(IdentQ,od)) + MOs::TensorProduct( od,MOs::TensorProduct(IdentQ,o)))
				+  0.8*(M_PI) *MOs::TensorProduct(IdentQ, ( MOs::TensorProduct(o, MOs::TensorProduct( IdentQ2,od)) + MOs::TensorProduct(od, MOs::TensorProduct( IdentQ2,o))))
				+  0.9*(M_PI) *MOs::TensorProduct(IdentQ, ( MOs::TensorProduct(Hcouple,IdentQ2)))
				+  0.7*(M_PI) *MOs::TensorProduct(IdentQ, ( MOs::TensorProduct(o, MOs::TensorProduct( IdentQ,MOs::TensorProduct(od, IdentQ))) + MOs::TensorProduct(od, MOs::TensorProduct( IdentQ,MOs::TensorProduct(o, IdentQ)))))
				
				+  0.5*(M_PI) * MOs::TensorProduct(Hcouple, IdentQ3) 
				;
		
		HDriftRef.SetOutputStyle(Matrix);	
		cout << HDriftRef << endl; 

/*testing exp of large herimitian matrix*/
		matrix<double>  p(refdim,refdim);		
		for(size_t j=0; j<refdim; j++)
			for(size_t k=0; k<refdim; k++)
				p(k,j) = abs(HDriftRef(k,j));
				
		cout << "exponent of HDriftRef is " << endl << endl; 
		HDriftRef = ExpM::EigenMethod(HDriftRef,-(1.0/5));
		cout << HDriftRef << endl; 
		return 0;		
/**/
	
}
