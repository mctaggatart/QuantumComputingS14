/*
Name: NumericalIntegration.hpp
Author: Jay Gambetta

Dependences: 
Brief Discription: This is my numerical integration class. I have set it up so that it takes in a function f(x) and performs numerical integration on this function. It does this by passing the function an address. 

It takes in a address to a function of the form f(x) and then does numerical integration on it. 

Limitations:

Version History
	v0: May 11th, 2009
*/


#ifndef NumericalIntegration_h
#define NumericalIntegration_h
#include "UsefullFunctions.hpp"

enum Steppertype  {Constant, Adaptive};

template <class T>
class Function {
	//This is the parent class of the equations. It has essentially empty elements which we will override with specific examples
	public:
		inline Function(){};
		virtual T function(const T x) = 0; //the function to be evaluated
		virtual ~Function(){};
};


template <class T>
class NI_solver
{
	public:
		//NI_solver(const T a, const T b, T (* function_ptr)(const T x));
		NI_solver(const T a, const T b, Function<T>* function_ptr);
		~NI_solver ();
		void SetNumericalParameters(const enum Steppertype steppertype, size_t& n);
		
		//Solvers. 
		T Trapezoidal();
		T Simpson();
		T Simpson3on8();
		T Booles();

	private: 												
		T a_, b_, h_; //time, start time, end time, saved stepes 
		size_t n_; //number of points in the integral
		enum Steppertype steppertype_; // the stepper used
		//T (* function_ptr_)(const double t);
		Function<T>* function_ptr_; // pointer to the function to be integrated

};
template <class T>
inline NI_solver<T>::NI_solver(const T a, const T  b, Function<T>* function_ptr){
	// declares the single integral solver
	a_ = a;
	b_ = b;
	function_ptr_=function_ptr;
	if(verbose==yes)
	{
		std::cout << "--------------------Simulation Parameters-------------------" << std::endl;
		std::cout << "a: " << a_ << std::endl;
		std::cout << "b: " << b_ << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}
template <class T>
inline NI_solver<T>::~NI_solver(){
	if(verbose==yes)
		std::cout << "Closing the solver: ~NI_solver" << std::endl; 
}
template <class T>
inline void NI_solver<T>::SetNumericalParameters(const enum Steppertype steppertype, size_t& n){
	n_=n;
	h_ = (b_-a_)/T(n_);
	steppertype_ = steppertype; // Sets the stepper type to the input type
	if(steppertype_ != Constant)
		UFs::MyError("NI_solver::SetNumericalParameters: you have not set numerical parameters correctly:\nUse \nxxx.SetNumericalParameters(Constant,h);\n");
	
	if(verbose==yes)
	{
	 	std::cout << "--------------------Numerical Parameters--------------------" << std::endl;
		std::cout << "step size: " << h_ << std::endl;
		std::cout << "A constant steeper must be used. " <<  std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}

template <class T>
inline T NI_solver<T>::Trapezoidal(){
// 	This routine integrate a function of one variable using the trapezoidal method.

    	T sum = function_ptr_->function(a_)*0.5;
	for (size_t i =1; i < n_; i++){
		sum += function_ptr_->function(a_+T(i)*h_);
	}
	sum += function_ptr_->function(a_+T(n_)*h_)*0.5;
   	return sum*h_;
} 
template <class T>
inline T NI_solver<T>::Simpson(){
//	This routine integrate a function of one variable using the Simpson method. 	
	T sum =  function_ptr_->function(a_)*0.5;
    	T sum_mid = 0.0;
    	for (size_t i =1; i <= n_; i++){
		sum += function_ptr_->function(a_+T(i)*h_);
		sum_mid += function_ptr_->function(a_+(T(i)-0.5)*h_);
	}
    	sum -=function_ptr_->function(a_+T(n_)*h_)*0.5;
    	return (sum+2.0*sum_mid)*h_/3.0;
}   
template <class T>
inline T NI_solver<T>::Simpson3on8(){
//	This routine integrate a function of one variable using the Simpson method. 	
	T sum =  function_ptr_->function(a_)*0.5;
    	T sum_mid = 0.0;
      	for (size_t i =1; i <= n_; i++){
		sum += function_ptr_->function(a_+T(i)*h_);
		sum_mid += function_ptr_->function(a_+(T(i)-1.0/3.0)*h_);
		sum_mid += function_ptr_->function(a_+(T(i)-2.0/3.0)*h_);
	}
    	sum -=function_ptr_->function(a_+T(n_)*h_)*0.5;
    	return (2.0*sum+3.0*sum_mid)*h_*0.125;
}   
template <class T>
inline T NI_solver<T>::Booles(){
//	This routine integrate a function of one variable using the Simpson method. 	
	T sum =  function_ptr_->function(a_)*0.5;
    	T sum_mid = 0.0;
      	for (size_t i =1; i <= n_; i++){
		sum += function_ptr_->function(a_+T(i)*h_);
		sum_mid += 16.0*function_ptr_->function(a_+(T(i)-0.25)*h_);
		sum_mid += 6.0*function_ptr_->function(a_+(T(i)-0.5)*h_);
		sum_mid += 16.0*function_ptr_->function(a_+(T(i)-0.75)*h_);
	}
    	sum -=function_ptr_->function(a_+T(n_)*h_)*0.5;
    	return (1.0/45.0)*(7.0*sum+sum_mid)*h_;
}   
#endif /* NumericalIntegration_h */
