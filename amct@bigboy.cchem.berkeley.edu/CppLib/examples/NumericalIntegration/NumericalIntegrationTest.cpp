/*
Name: FILENAME.cpp
Author: Jay Gambetta

Dependences: NumericalIntegration.hpp
Brief Discription: This program demonstrates my numerical integration class. 

Limitations: None

Version History
    v0: May 11th, 2009.

*/
#include "NumericalIntegration.hpp"
using namespace std; 

template <class T>
class Function_child: public Function<double> {
	//This is the parent class of the equations. It has essentially empty elements which we will override with specific examples
	public:
		inline Function_child(){};
		T function(const T x); //the function to be evaluated
		virtual ~Function_child(){};
};

template <class T>
inline T Function_child<T>::function(const T x){
	return exp(-2.0*x)+ cos(10.0*x);
}

int main (int argc, char const *argv[]){
	Function_child<double> sys;
	//Function<double> * pointersys = &sys;
	double a = 0.0;
	double b = 10.0;
	size_t n = 1000.0;
	NI_solver<double > solver(a, b, &sys);
	solver.SetNumericalParameters(Constant, n);
	cout << " Trapezoidal: " << solver.Trapezoidal() << " Simpson: " << solver.Simpson() << " Simpson3/8: " << solver.Simpson3on8() << " Booles: " << solver.Booles() << '\t' << (1.0-exp(-2.0*10.0))/2.0 + 1.0/10.0*(sin(10.0*10.0)-0.0)<< endl;
	return 0;
}
