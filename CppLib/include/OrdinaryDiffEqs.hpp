/*
Name: OrdinaryDiffEqs.hpp
Author: Jay Gambetta

Dependences: UsefullFunctionss.hpp
Brief Discription: This Header contians all the ordinary differential equation solvers I have written, the classes  it takes in must have access by y[i] and have a function y.size(); which returns the length of the arrary. In works under the following procedure It has a class ODE_solver which cotains all the solvers that will be used by the class. 

To allow inputing of different functions i use another class called ODE_equations it contains a contructor, distructor and some virtual functions derives and savesubset. In the actually code 
	ODE_equations_child sys;
	ODE_equations<T >* pointersys = &sys; 
	which produces a pointer which is of type ODE_epuations and points to the address of the child euquations which contain the actuall functions

	I am thinking about removing the need for the ODE_equations class and replacing this with just an input to the solver, i.e. the address of a function derivs and savesubset. 


Limitations: the constant steppers dont stop on the correct time and the adaptive do some 
	unnessary steps (put some of the lines in the if loop so that they are only done if we are goint to end before a set time)

Version History
	v0: June 15th, 2007. Started the class file
	v1: Janurary 3rd, 2009. Change to be consistents with the new library
*/


#ifndef OrdinaryDiffEqs_h
#define OrdinaryDiffEqs_h
#include "UsefullFunctions.hpp"

enum Steppertype  {Constant, Adaptive};
enum Savetype {Complete, Testing, Subset};

template <class T>
class ODE_equations {
	//This is the parant class of the equations. It has essentially empty elements which we will override with specific examples
	public:
		inline ODE_equations(){};
		virtual void derivs(const double x, const T& y, T& dy) = 0; //updates dy with A(y) ie dot y = A(y)
		virtual void SaveSubset(std::vector<double>& data, const T& y) = 0; // saves y to data under user specified processing
		virtual ~ODE_equations(){};
};

template <class T>
class ODE_solver {
	//This is the solver class - it constructs a solve for a specific ODE given by the class ODE_equations
	//It contains sever types of solves
	public:
		// Constructs the class with a given starttime, endtime, savesteps (number of output points), A set of equations and a filename to save to
		ODE_solver(const double starttime, const double endtime, const T& yin, ODE_equations<T>* pointerequation, const size_t savesteps, const std::string outfile);
		virtual ~ODE_solver(); 	
		void SetNumericalParameters(const enum Steppertype steppertype, double& h);
		void SetNumericalParameters(const enum Steppertype steppertype, double& h, double const eps_abs, double const eps_rel, double const scale_y, double const scale_dydt);
		void SetDataSaveType(const enum Savetype savetype);
		
		//Solvers. 
		void EulerConstStep();
		void EulerConstStepRomberg();
		void RK4ConstStep();
		void RK4ConstStepRomberg();
		void EulerAdapStep();
		void RK4AdapStep();
		void RK45AdapStep();
		void RKCKAdapStep();
		
		//Functions used by some of the sovlers
		void RK4TwoStep(double const, T& yout_a, T& yout_b, double const);
		void RK45OneStep(double const x, T& yout, T& yerr, double const h);
		void RKCKOneStep(double const x, T& yout, T& yerr, double const h);
		
		//Need to add a reset for y and t

	private:
		size_t savesteps_; // the number of steps saved, 												
		double t, starttime_, endtime_, h_, dtsave_; //time, start time, end time, saved stepes 
		double eps_abs_, eps_rel_, scale_y_, scale_dydt_; // adaptive parameters
		T y_; // the state
		size_t dim_; // the dimensions of the system to be solved
		ODE_equations<T>* pointerequation_; // pointer to the ordinary differential equation(s) to be solved
		std::ofstream data_out_;  // the output data stream
		enum Steppertype steppertype_;
		enum Savetype savetype_;
		void SaveData();
};
template <class T>
inline ODE_solver<T>::ODE_solver(const double starttime, const double endtime, const T& yin, ODE_equations<T>* pointerequation, const size_t savesteps, const std::string outfile){
	starttime_ = starttime;
	endtime_ = endtime;
	savesteps_ = savesteps;
	dtsave_ = (endtime_-starttime_)/double(savesteps_); 
	dim_=yin.size();
	y_ = yin;
	t =  starttime;
	pointerequation_ = pointerequation; //pointer to the ODE equations
	savetype_=Complete;  //Sets the save type to default only save time
	UFs::OpenFile(outfile, data_out_, 16); 	// opens a file for saving the data
	if(verbose==yes)
	{
		std::cout << "--------------------Simulation Parameters-------------------" << std::endl;
		std::cout << "start time: " << starttime_ << std::endl;
		std::cout << "end time: " << endtime_ << std::endl;
		std::cout << "data will be save evey dt " << dtsave_ << " giving " << savesteps_ << " data points" << std::endl;
		std::cout << "data save file: " << outfile << std::endl;
		std::cout << "the inner dimensions of the vector is " << dim_ << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}
template <class T>
inline ODE_solver<T>::~ODE_solver(){
	if(verbose==yes)
		std::cout << "Closing the solver: ~ODE_solver" << std::endl; 
	data_out_.close();
}
template <class T>
inline void ODE_solver<T>::SetNumericalParameters(const enum Steppertype steppertype, double& h){
	h_ = h;
	steppertype_ = steppertype; // Sets the stepper type to the input type
	if(steppertype_ != Constant)
		UFs::MyError("ODE_solver::SetNumericalParameters: you have not set numerical parameters correctly:\nUse \nxxx.SetNumericalParameters(Constant,h);\n");
	
	if(verbose==yes)
	{
	 	std::cout << "--------------------Numerical Parameters--------------------" << std::endl;
		std::cout << "step size: " << h_ << std::endl;
		std::cout << "A constant steeper must be used. " <<  std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}
template <class T>
inline void ODE_solver<T>::SetNumericalParameters(const enum Steppertype steppertype, double& h, double const eps_abs, double const eps_rel, double const scale_y, double const scale_dydt){
	h_ = h;
	eps_abs_ = eps_abs;
	eps_rel_ = eps_rel;
	scale_y_ = scale_y;
	scale_dydt_ = scale_dydt;
	steppertype_ = steppertype; // Sets the stepper type to the input type
	if(steppertype_ != Adaptive)
		UFs::MyError("ODE_solver::SetNumericalParameters: you have not set numerical parameters correctly:\nUse \nxxx.SetNumericalParameters(Addative,h.eps_abs,eps_rel,scale_y,scale_dydt);\n");
	
	if(verbose==yes)
	{
		std::cout << "--------------------Numerical Parameters--------------------" << std::endl;
		std::cout << "initial step size: " << h_ << std::endl;
		std::cout << "the absolute error: " << eps_abs_ << std::endl;
		std::cout << "the realtive error: " << eps_rel_ << std::endl;
		std::cout << "the scaling factor for the y: " << scale_y_ << std::endl;
		std::cout << "the scaling factor of the dydx: " << scale_dydt_ << std::endl;
		std::cout << "An adaptive steeper must be used. " <<  std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}
template <class T>
inline void ODE_solver<T>::SetDataSaveType(enum Savetype savetype){
	savetype_ = savetype;
	if(verbose==yes)
	{
		std::cout << "---------------------Data Saved Type-------------------------" << std::endl;
		std::cout << "The data has been set to save type: "<< savetype_ << std::endl;
		std::cout << "-------------------------------------------------------------" << std::endl;
	}
}
template <class T>
inline void ODE_solver<T>::SaveData(){
	data_out_ << t << "\t";
	std::vector<double> data;
	
	switch (savetype_){
	case Complete: 
		for (size_t j=0; j < dim_; j++){
			data_out_ << y_[j] << "\t";
		}
		break;
	case Subset:
		//Saves the operators defined in ODE_equations_child::SaveSubset
		pointerequation_->SaveSubset(data, y_);
		for (size_t j =0; j < data.size(); j++){
			data_out_ << data[j] << "\t";
		}
		break;
	case Testing: 
	data_out_ << h_ << "\t";
		for (size_t j=0; j < dim_; j++){
			data_out_ << y_[j] << "\t";
		}
		break;
	default:
		UFs::MyError("Routine ODE_solver::SaveData: No valide save data type set:\n");
	}
	data_out_ << std::endl;
}	
// ==================== Constant Solvers =========================
template <class T>
inline void ODE_solver<T>::EulerConstStep(){
	if(steppertype_ != Constant)
		UFs::MyError("ODE_solver::EulerConstStep(): you have not set the numerical parameters:\nUse \nxxx.SetNumericalParameters(Constant,h);\n");
	T a(dim_);
	double tstop = starttime_;
	//ouput intital conditons
	ODE_solver<T>::SaveData();
	for(size_t k=0; k < savesteps_; k++){
		tstop += dtsave_;
		while (t < tstop){
			pointerequation_->derivs(t,y_,a);
  			for (size_t j=0; j <dim_; j++){
				y_[j] = y_[j]+h_*a[j];
			}
			t += h_;
		}
		ODE_solver<T>::SaveData();
	}
}
template <class T>
inline void ODE_solver<T>::EulerConstStepRomberg(){
	if(steppertype_ != Constant)
		UFs::MyError("ODE_solver::EulerConstStepRomberg(): you have not set the numerical parameters:\nUse \nxxx.SetNumericalParameters(Constant,h);\n");
	T a(dim_), ytemp(dim_);
	double h2 = 2*h_;
	const double Beta = 1.0; // 1/ (2^gamma-1)
	
	double tstop = starttime_;
	//ouput intital conditons
	ODE_solver<T>::SaveData();
	for(size_t k=0; k < savesteps_; k++){
		tstop += dtsave_;
		while (t < tstop){
			pointerequation_->derivs(t,y_,a);
			for (size_t j=0; j < dim_; j++){
				ytemp[j] = y_[j]+h2*a[j];  //ytemp acting as y1
				y_[j] = y_[j]+h_*a[j];  //y acting as y2
			}
			pointerequation_->derivs(t+h_,y_,a);
			for (size_t j=0; j < dim_; j++){
				y_[j] = y_[j] + h_*a[j];  //y acting as y2
				y_[j] = y_[j] +  (y_[j] - ytemp[j])*Beta;
			}
			t += h2;
		}
		ODE_solver<T>::SaveData();
	}
}
template <class T>
inline void ODE_solver<T>::RK4ConstStep(){
	if(steppertype_ != Constant)
		UFs::MyError("ODE_solver::RK4ConstStep(): you have not set the numerical parameters:\nUse \nxxx.SetNumericalParameters(Constant,h);\n");
	size_t j;
	T k1(dim_), k_temp1(dim_), k_temp2(dim_), ytemp(dim_);
	double h_half= 0.5*h_, h6 = h_/6.0;
	
	double tstop = starttime_;
	//ouput intital conditons
	ODE_solver<T>::SaveData();
	for(size_t k=0; k < savesteps_; k++){
		tstop += dtsave_;
		while (t < tstop){
			//step 1 (half step) 
			pointerequation_->derivs(t,y_,k1);
			for (j=0; j < dim_; j++){
			 	ytemp[j] = y_[j] + h_half*k1[j];
			}
			//step 2
			pointerequation_->derivs(t+h_half,ytemp,k_temp1); // returns k2 = f(t+h/2, y+h*k1/2) k_temp1 act as k2 here
			for (j=0; j < dim_; j++){
				ytemp[j] = y_[j] + h_half*k_temp1[j];
			}
			//step 3
			pointerequation_->derivs(t+h_half,ytemp,k_temp2); //returns k3 = f(t+h/2. y+h*k2/2) k_temp2 acts as k3 here 
			for (j=0; j < dim_; j++){
				ytemp[j] = y_[j] + h_*k_temp2[j];
				k_temp2[j] += k_temp1[j];     // k_temp2 is acting as k2 + k3 
			}
			//step 4
			pointerequation_->derivs(t+h_,ytemp,k_temp1); // returns k4 = f(t+h, y+h*k3) k_temp1 is acting as k4 now
			for (j=0; j < dim_; j++){
				y_[j] = y_[j] + h6*(k1[j] + k_temp1[j] + 2.0*k_temp2[j]);
			} // yout = y + h*(k1+2.0*k2+2.0*k3+k4)/6
			t += h_;
		}
		ODE_solver<T>::SaveData();
	}
}
template <class T>
inline void ODE_solver<T>::RK4TwoStep(double const x, T& yout_a, T& yout_b, double const h){
	
	size_t i;
	double h_half= 0.5*h, h_on_6 = h/6.0, h2 = 2.0*h, h_on_3 = h/3.0, h_1_and_half = 1.5*h;
	T k1(dim_), k_temp1_a(dim_), k_temp2_a(dim_), ytemp_a(dim_), k_temp1_b(dim_), k_temp2_b(dim_), ytemp_b(dim_);
     
	// b is the y2 case and the a is the y1 case
	// k1 = f(x,y)
	pointerequation_->derivs(x,y_,k1);
	for (i=0; i < dim_; i++){
		
		ytemp_b[i] = y_[i] + h_half*k1[i];
		ytemp_a[i] = y_[i] + h*k1[i]; 	
	}
	pointerequation_->derivs(x+h_half,ytemp_b,k_temp1_b); // returns k2 = f(x+h/2, y+h*k1/2) k_temp1_b act as k2 here 
	pointerequation_->derivs(x+h,ytemp_a,k_temp1_a); // returns k2 = f(x+h/2, y+h*k1/2) k_temp1_a act as k2 here 
	for (i=0; i < dim_; i++){

		ytemp_b[i] = y_[i] + h_half*k_temp1_b[i];
		ytemp_a[i] = y_[i] + h*k_temp1_a[i];  
	}
	pointerequation_->derivs(x+h_half,ytemp_b,k_temp2_b); //returns k3 = f(x+h/2. y+h*k2/2) k_temp2_b acts as k3 here
	pointerequation_->derivs(x+h,ytemp_a,k_temp2_a); //returns k3 = f(x+h/2. y+h*k2/2) k_temp2_a acts as k3 here 
	for (i=0; i < dim_; i++){
		
		ytemp_b[i] = y_[i] + h*k_temp2_b[i];
		k_temp2_b[i] += k_temp1_b[i];     // k_temp2_b is acting as k2 + k3 

		ytemp_a[i] = y_[i] + h2*k_temp2_a[i];
		k_temp2_a[i] += k_temp1_a[i];     // k_temp2_a is acting as k2 + k3 
	
	}
	pointerequation_->derivs(x+h,ytemp_b,k_temp1_b); // returns k4 = f(x+h, y+h*k3) k_temp1_b is acting as k4 now
	pointerequation_->derivs(x+h2,ytemp_a,k_temp1_a); // returns k4 = f(x+h, y+h*k3) k_temp1_a is acting as k4 now
	for (i=0; i < dim_; i++){
		
		yout_b[i] = y_[i] + h_on_6*(k1[i] + k_temp1_b[i] + 2.0*k_temp2_b[i]);
		yout_a[i] = y_[i] + h_on_3*(k1[i] + k_temp1_a[i] + 2.0*k_temp2_a[i]);

	} // yout = y + h*(k1+2.0*k2+2.0*k3+k4)/6
	
	//second step
	pointerequation_->derivs(x+h,yout_b,k1); // returns k1 = f(x, y) k1 act as k1 here 
	for (i=0; i< dim_; i++){
		
	 	ytemp_b[i] = yout_b[i] + h_half*k1[i];
	}
	pointerequation_->derivs(x+h_1_and_half,ytemp_b,k_temp1_b); // returns k2 = f(x+h/2, y+h*k1/2) k_temp1_b act as k2 here 
	
	for (i=0; i< dim_; i++){
		
		ytemp_b[i] = yout_b[i]+h_half*k_temp1_b[i];
	}
	pointerequation_->derivs(x+h_1_and_half,ytemp_b,k_temp2_b); //returns k3 = f(x+h/2. y+h*k2/2) k_temp2_b acts as k3 here 
	for (i=0; i< dim_; i++){
		
		ytemp_b[i] = yout_b[i] + h*k_temp2_b[i];
		k_temp2_b[i] += k_temp1_b[i];     // k_temp2_b is acting as k2 + k3 
	}
	pointerequation_->derivs(x+h2,ytemp_b,k_temp1_b); // returns k4 = f(x+h, y+h*k3) k_temp1_b is acting as k4 now
	for (i=0; i < dim_; i++){
		
		yout_b[i] = yout_b[i] + h_on_6*(k1[i] + k_temp1_b[i] + 2.0*k_temp2_b[i]);
	} // yout = y + h*(k1+2.0*k2+2.0*k3+k4)/6


}
template <class T>
inline void ODE_solver<T>::RK4ConstStepRomberg(){
	if(steppertype_ != Constant)
		UFs::MyError("ODE_solver::RK4ConstStepRomberg(): you have not set the numerical parameters:\nUse \nxxx.SetNumericalParameters(Constant,h);\n");
	T y_1(dim_), y_2(dim_);
	const double h2 = 2*h_;
	const double Beta = 1/15.0; // 1/ (2^gamma-1)
	
	double tstop = starttime_;
	//ouput intital conditons
	ODE_solver<T>::SaveData();
	for(size_t k=0; k < savesteps_; k++){
		tstop += dtsave_;
		while (t < tstop){
		
			ODE_solver<T>::RK4TwoStep(t,y_1,y_2,h_);
		
			for (size_t j=0; j < dim_; j++){
				y_[j] = y_2[j] +  (y_2[j] - y_1[j])*Beta;
			}
			t += h2;
		}
		ODE_solver<T>::SaveData();
	}
}
// ==================== Apative Solvers ==========================
template <class T>
inline void ODE_solver<T>::EulerAdapStep(){
	
	// eps_abs - the absolute error
	// eps_rel - the relative error
	// scale_y - the scaling factor for the y
	// scale_dydx - the scaling factor for the dydx
	if(steppertype_ != Adaptive)
		UFs::MyError("ODE_solver::EulerAdapStep: you have not set numerical parameters:\nUse \nxxx.SetNumericalParameters(Addative,h.eps_abs,eps_rel,scale_y,scale_dydt);\n");
		
	size_t j;
	T a(dim_), yerr(dim_), ytemp(dim_), derr(dim_);
	const double Beta = 1.0; // 1/(2^gamma-1)
	const double SAFTY = 0.9; 
	const double PSHRNK = -1; // - 1/(gamma)
	const double PGROW = -0.5; // -1/(gamma+1)
	const double TINY = 1.0e-16;
	double errmax, htemp, h2;
	double ttemp, tstop_minus_tiny, hlast, h2last; 
	
	double tstop = starttime_;
	//ouput intital conditons
	ODE_solver<T>::SaveData();
	for(size_t k=0; k < savesteps_; k++){
		tstop += dtsave_;
		tstop_minus_tiny = tstop - TINY;
		//the stepper
		while (t < tstop_minus_tiny){

			for(;;){
				//Take a Euler step 
				pointerequation_->derivs(t,y_,a);
				h2 = 2*h_;	
				for (j=0; j <dim_; j++){
					yerr[j] = y_[j]+h2*a[j]; //yerr acting as y1
					ytemp[j]= y_[j]+h_*a[j];  //ytemp acting as y2
					// the desired error
					// THIS PART IS NOT UNIVERSAL The abs subroutine computes the complex absolute value of x
					// Compute absolute value
					derr[j] = eps_abs_ + eps_rel_*(scale_y_*abs(y_[j])+scale_dydt_*abs(a[j])*h_);	
				}
				pointerequation_->derivs(t+h_,ytemp,a);
				//getting the accuracy and ytemp
				errmax  = 0.0; 
				for (j=0; j < dim_; j++){
					ytemp[j]= ytemp[j]+h_*a[j]; //ytemp acting as y2
					yerr[j] = ytemp[j] - yerr[j]; // the error
					// THIS PART IS NOT UNIVERSAL
					errmax = USs::Max<double>(errmax, abs(yerr[j]/derr[j])); // need to scale the error
					ytemp[j] = ytemp[j] + yerr[j]*Beta; // ytemp acting as the estimate
				}
				//If the observed error yerr  exceeds the desired error derr  by more than 10% for any component 
				//then the method reduces the step-size by an appropriate factor
				if (errmax <= 1.1) break; // step succeeded. Jump to calculation of the next step
			
				// Reduction of the stepper which goes as much as 10
				htemp = SAFTY*h_*pow(errmax,PSHRNK);
				h_=USs::Max<double>(htemp,0.1*h_);
				if (h_ < 0) 	UFs::MyError("Routine ODE_solver::EulerAdapStep: h went below zero");
				if (t+h_ == t) UFs::MyError("Routine ODE_solver::EulerAdapStep: stepsize underflow");
			}
		
			//Apply the step
			ttemp = t+2.0*h_;
			if (ttemp < tstop_minus_tiny){
				for (j=0; j < dim_; j++){
					y_[j] = ytemp[j];
				}
				t = ttemp;
			}
			else {
				hlast = (tstop - t)*0.5;
				h2last = 2*hlast; 
				pointerequation_->derivs(t,y_,a);
				for (j=0; j < dim_; j++){
					ytemp[j] = y_[j]+h2last*a[j];  //ytemp acting as y1
					y_[j] = y_[j]+hlast*a[j];  //y acting as y2
				}
				pointerequation_->derivs(t+hlast,y_,a);
				for (j=0; j < dim_; j++){
					y_[j] = y_[j] + hlast*a[j];  //y acting as y2
					y_[j] = y_[j] +  (y_[j] - ytemp[j])*Beta;
				}
				t = tstop;
			}		

			//Calculation of the next step
			//If the observed error yerr  is less than 50% of the desired error derr for the maximum ratio errmax 
			//then the algorithm takes the opportunity to increase the step-size to bring the error in 
			//line with the desired level,
			if (errmax <= 0.5){
				htemp = SAFTY*h_*pow(errmax,PGROW);
				//Maximum growth is 5
				h_=USs::Min<double>(htemp,5.0*h_);
				if (h_ < 0)  UFs::MyError("Routine ODE_solver::EulerAdapStep:: h went below zero");
			}
	 	}
		ODE_solver<T>::SaveData();
	//	cout << "h is " << h_;
	}
}
template <class T>
inline void ODE_solver<T>::RK4AdapStep(){
	
	
	// eps_abs - the absolute error
	// eps_rel - the relative error
	// scale_y - the scaling factor for the y
	// scale_dydx - the scaling factor for the dydx
	if(steppertype_ != Adaptive)
		UFs::MyError("ODE_solver::RK4AdapStep: you have not set numerical parameters:\nUse \nxxx.SetNumericalParameters(Addative,h.eps_abs,eps_rel,scale_y,scale_dydt);\n");
			
	size_t j;
	T a(dim_), yerr(dim_), ytemp(dim_), derr(dim_);
	const double Beta = 1/15.0; // 1/(2^gamma-1)
	const double SAFTY = 0.9; 
	const double PSHRNK = -0.25; // - 1/(gamma)
	const double PGROW = -0.2; // -1/(gamma+1)
	const double TINY = 1.0e-16;
	double errmax, htemp;
	double ttemp, tstop_minus_tiny, hlast; 
	
	double tstop = starttime_;
	//ouput intital conditons
	ODE_solver<T>::SaveData();
	for(size_t k=0; k < savesteps_; k++){
		tstop += dtsave_;
		tstop_minus_tiny = tstop - TINY;
		//the stepper
		while (t < tstop_minus_tiny){

			for(;;){
				//Take one step step 
				pointerequation_->derivs(t,y_,a);
				// the desired error
				for (j=0; j < dim_; j++){
					derr[j] = eps_abs_ + eps_rel_*(scale_y_*abs(y_[j])+scale_dydt_*abs(a[j])*h_);	
				}
				ODE_solver::RK4TwoStep(t,yerr,ytemp,h_); // /yerr acting as y1 and ytemp acting as y2
				//getting the accuracy and ytemp
				errmax  = 0.0; 
				for (j=0; j <dim_; j++){
					yerr[j] = ytemp[j] - yerr[j]; // the error
					errmax = USs::Max<double>(errmax, abs(yerr[j]/derr[j])); // need to scale the error
					ytemp[j] = ytemp[j] + yerr[j]*Beta; // ytemp acting as the estimate
				}
				//If the observed error yerr  exceeds the desired error derr  by more than 10% for any component 
				//then the method reduces the step-size by an appropriate factor
				if (errmax <= 1.1) break; // step succeeded. Jump to calculation of the next step
			
				// Reduction of the stepper which goes as much as 10
				htemp = SAFTY*h_*pow(errmax,PSHRNK);
				h_=USs::Max<double>(htemp,0.1*h_);
				if (h_ < 0) 	UFs::MyError("Routine ODE_solver::RK4AdapStep: : h went below zero");
				if (t+h_ == t) UFs::MyError("Routine ODE_solver::RK4AdapStep: : stepsize underflow");
			}
		
			//Apply the step
			ttemp = t+2.0*h_;
			if (ttemp < tstop_minus_tiny){
				for (j=0; j < dim_; j++){
					y_[j] = ytemp[j];
				}
				t = ttemp;
			}
			else {
				hlast = (tstop - t)*0.5;
				pointerequation_->derivs(t,y_,a);
				ODE_solver<T>::RK4TwoStep(t,yerr,ytemp,hlast); // /yerr acting as y1 and ytemp acting as y2
				for (j=0; j < dim_; j++){
					y_[j] = ytemp[j] +  (ytemp[j] - yerr[j])*Beta;
				}
				t = tstop;
			}		

			//Calculation of the next step
			//If the observed error yerr  is less than 50% of the desired error derr for the maximum ratio errmax 
			//then the algorithm takes the opportunity to increase the step-size to bring the error in 
			//line with the desired level,
			if (errmax <= 0.5){
				htemp = SAFTY*h_*pow(errmax,PGROW);
				//Maximum growth is 5
				h_=USs::Min<double>(htemp,5.0*h_);
				if (h_ < 0)  UFs::MyError("Routine ODE_solver::RK4AdapStep: : h went below zero");
			}
		}
		ODE_solver<T>::SaveData();
		//	cout << "h is " << h_;
	}
}
// ==================== Apative Solvers embedded =================
template <class T>
inline void ODE_solver<T>::RK45OneStep(double const x, T& yout, T& yerr, double const h){
	
	static const double b2=0.25, b3=3.0/8.0, b4=12.0/13.0, b5=1.0, b6=0.5, b21=0.25, b31=3.0/32.0, b32=9.0/32.0, b41=1932.0/2197.0, b42=-7200.0/2197.0, b43=7296.0/2197.0, b51=439.0/216.0, b52=-8.0, b53=3680.0/513.0, b54=-845.0/4104.0, b61=-8.0/27.0, b62=2.0, b63=-3544.0/2565.0, b64=1859.0/4104.0, b65=-11.0/40.0, c1=16.0/135.0, c3=6656.0/12825.0, c4=28561.0/56430.0, c5=-9.0/50.0, c6=2.0/55.0, dc1=c1-25.0/216.0, dc3=c3-1408.0/2565.0, dc4=c4-2197.0/4104.0, dc5=c5-1.0/5.0, dc6=c6;

	size_t i;
	
	T ak1(dim_), ak2(dim_), ak3(dim_), ak4(dim_), ak5(dim_), ak6(dim_), ytemp(dim_);
    
	//first step
	pointerequation_->derivs(t,y_,ak1);
	for (i=0;i<dim_;i++){
		ytemp[i] = y_[i] + b21*h*ak1[i];
	}
	// second step
	pointerequation_->derivs(x+b2*h,ytemp,ak2);
	for (i=0;i<dim_;i++){
		ytemp[i] = y_[i] + h*(b31*ak1[i] + b32*ak2[i]);
	}
	// third step
	pointerequation_->derivs(x+b3*h,ytemp,ak3);
	for (i=0; i<dim_; i++){
		ytemp[i] = y_[i] + h*(b41*ak1[i] + b42*ak2[i] + b43*ak3[i]); 	
	}
	// fourth step
	pointerequation_->derivs(x+b4*h,ytemp,ak4);  
	for (i=0; i<dim_; i++){
		ytemp[i] = y_[i] + h*(b51*ak1[i] + b52*ak2[i] + b53*ak3[i] + b54*ak4[i]);  
	}
	// fifth step 
	pointerequation_->derivs(x+b5*h,ytemp,ak5); 
	for (i=0; i<dim_; i++){
		ytemp[i] = y_[i] + h*(b61*ak1[i] + b62*ak2[i] + b63*ak3[i] + b64*ak4[i] + b65*ak5[i]); 
	}
	// sixth step
	pointerequation_->derivs(x+b6*h,ytemp,ak6);
	for (i=0; i <dim_; i++){
		yout[i] = y_[i] + h*(c1*ak1[i] + c3*ak3[i] + c4*ak4[i] + c5*ak5[i] + c6*ak6[i]); // The Accumulate increments with proper weights
		yerr[i] = h*(dc1*ak1[i] + dc3*ak3[i] + dc4*ak4[i] + dc5*ak5[i] + dc6*ak6[i]); // Estimate error as difference between fourth and fifth order methods.
	} 
}
template <class T>
inline void ODE_solver<T>::RK45AdapStep(){
	
	// eps_abs - the absolute error
	// eps_rel - the relative error
	// scale_y - the scaling factor for the y
	// scale_dydx - the scaling factor for the dydx
	if(steppertype_ != Adaptive)
		UFs::MyError("ODE_solver::RK45AdapStep: you have not set numerical parameters:\nUse \nxxx.SetNumericalParameters(Addative,h.eps_abs,eps_rel,scale_y,scale_dydt);\n");
		
	size_t j;
	T a(dim_), yerr(dim_), ytemp(dim_), derr(dim_);
	const double SAFTY = 0.9; 
	const double PSHRNK = -0.25; // - 1/(gamma)
	const double PGROW = -0.2; // -1/(gamma+1)
	const double TINY = 1.0e-16;
	double errmax, htemp;
	double ttemp, tstop_minus_tiny, hlast; 
	
	double tstop = starttime_;
	//ouput intital conditons
	ODE_solver<T>::SaveData();
	for(size_t k=0; k < savesteps_; k++){
		tstop += dtsave_;
		tstop_minus_tiny = tstop - TINY;
		//the stepper
		while (t < tstop_minus_tiny){
			
			for(;;){
				//Take one step step 
				pointerequation_->derivs(t,y_,a);
				// the desired error
				for (j=0; j < dim_; j++){
					derr[j] = eps_abs_ + eps_rel_*(scale_y_*abs(y_[j])+scale_dydt_*abs(a[j])*h_);	
				}
				ODE_solver<T>::RK45OneStep(t,ytemp,yerr,h_);
				//getting the accuracy and ytemp
				errmax  = 0.0; 
				for (j=0; j < dim_; j++){
					errmax = USs::Max<double>(errmax, abs(yerr[j]/derr[j])); // need to scale the error
				}
				//If the observed error yerr exceeds the desired error derr by more than 10% for any component 
				//then the method reduces the step-size by an appropriate factor
				if (errmax <= 1.1) break; // step succeeded. Jump to calculation of the next step
				
				// Reduction of the stepper which goes as much as 10
				htemp = SAFTY*h_*pow(errmax,PSHRNK);
				h_=USs::Max<double>(htemp,0.1*h_);
				if (h_ < 0) UFs::MyError("Routine ODE_solver::RK45AdapStep: h went below zero");
				if (t+h_ == t) UFs::MyError("Routine ODE_solver::RK45AdapStep: stepsize underflow");
			}
		
			//Apply the step
			ttemp = t+h_;
			if (ttemp < tstop_minus_tiny){
				for (j=0; j < dim_; j++){
					y_[j] = ytemp[j];
				}
				t = ttemp;
			}
			else {
				hlast = (tstop - t);
				pointerequation_->derivs(t,y_,a);
				ODE_solver<T>::RK45OneStep(t,ytemp,yerr,hlast);
				for (j=0; j < dim_; j++){
					y_[j] = ytemp[j];
				}
				t = tstop;
			}		

			//Calculation of the next step
			//If the observed error yerr  is less than 50% of the desired error derr for the maximum ratio errmax 
			//then the algorithm takes the opportunity to increase the step-size to bring the error in 
			//line with the desired level,
			if (errmax <= 0.5){
				htemp = SAFTY*h_*pow(errmax,PGROW);
				//Maximum growth is 5
				h_=USs::Min<double>(htemp,5.0*h_);
				if (h_ < 0)  UFs::MyError("Routine ODE_solver::RK45AdapStep: h went below zero");
			}
		}
		ODE_solver<T>::SaveData();
		//	cout << "h is " << h_;
	}
}
template <class T>
inline void ODE_solver<T>::RKCKOneStep(double const x, T& yout, T& yerr, double const h){
	
	static const double b2=0.2, b3=0.3, b4=0.6, b5=1.0, b6=0.875, b21=0.2, b31=3.0/40.0, b32=9.0/40.0, b41=0.3, b42=-0.9, b43=1.2, b51=-11.0/54.0, b52=2.5, b53=-70.0/27.0, b54=35.0/27.0, b61=1631.0/55296.0, b62=175.0/512.0, b63=575.0/13824.0, b64=44275.0/110592.0, b65=253.0/4096.0, c1=37.0/378.0, c3=250.0/621.0, c4=125.0/594.0, c6=512.0/1771.0, dc1=c1-2825.0/27648.0, dc3=c3-18575.0/48384.0, dc4=c4-13525.0/55296.0, dc5=-277.0/14336.0, dc6=c6-0.25;

	size_t i;
	
	T ak1(dim_), ak2(dim_), ak3(dim_), ak4(dim_), ak5(dim_), ak6(dim_), ytemp(dim_);
    
	//first step
	pointerequation_->derivs(t,y_,ak1);
	for (i=0;i<dim_;i++){
		ytemp[i] = y_[i] + b21*h*ak1[i];
	}
	// second step
	pointerequation_->derivs(x+b2*h,ytemp,ak2);
	for (i=0;i<dim_;i++){
		ytemp[i] = y_[i] + h*(b31*ak1[i] + b32*ak2[i]);
	}
	// third step
	pointerequation_->derivs(x+b3*h,ytemp,ak3);
	for (i=0; i<dim_; i++){
		ytemp[i] = y_[i] + h*(b41*ak1[i] + b42*ak2[i] + b43*ak3[i]); 	
	}
	// fourth step
	pointerequation_->derivs(x+b4*h,ytemp,ak4);  
	for (i=0; i<dim_; i++){
		ytemp[i] = y_[i] + h*(b51*ak1[i] + b52*ak2[i] + b53*ak3[i] + b54*ak4[i]);  
	}
	// fifth step 
	pointerequation_->derivs(x+b5*h,ytemp,ak5); 
	for (i=0; i<dim_; i++){
		ytemp[i] = y_[i] + h*(b61*ak1[i] + b62*ak2[i] + b63*ak3[i] + b64*ak4[i] + b65*ak5[i]); 
	}
	// sixth step
	pointerequation_->derivs(x+b6*h,ytemp,ak6);
	for (i=0; i <dim_; i++){
		yout[i] = y_[i] + h*(c1*ak1[i] + c3*ak3[i] + c4*ak4[i] + c6*ak6[i]); // The Accumulate increments with proper weights
		yerr[i] = h*(dc1*ak1[i] + dc3*ak3[i] + dc4*ak4[i] + dc5*ak5[i] + dc6*ak6[i]); // Estimate error as difference between fourth and fifth order methods.
	} 
}   
template <class T>
inline void ODE_solver<T>::RKCKAdapStep(){
	
	// eps_abs - the absolute error
	// eps_rel - the relative error
	// scale_y - the scaling factor for the y
	// scale_dydx - the scaling factor for the dydx
	if(steppertype_ != Adaptive)
		UFs::MyError("ODE_solver::RKCKAdapStep: you have not set numerical parameters:\nUse \nxxx.SetNumericalParameters(Addative,h.eps_abs,eps_rel,scale_y,scale_dydt);\n");
		
	size_t j;
	T a(dim_), yerr(dim_), ytemp(dim_), derr(dim_);
	const double SAFTY = 0.9; 
	const double PSHRNK = -0.25; // - 1/(gamma)
	const double PGROW = -0.2; // -1/(gamma+1)
	const double TINY = 1.0e-16;
	double errmax, htemp;
	double ttemp, tstop_minus_tiny, hlast; 
	
	double tstop = starttime_;
	//ouput intital conditons
	ODE_solver<T>::SaveData();
	for(size_t k=0; k < savesteps_; k++){
		tstop += dtsave_;
		tstop_minus_tiny = tstop - TINY;
		//the stepper
		while (t < tstop_minus_tiny){
			
			for(;;){
				//Take one step step 
				pointerequation_->derivs(t,y_,a);
				// the desired error
				for (j=0; j < dim_; j++){
					derr[j] = eps_abs_ + eps_rel_*(scale_y_*abs(y_[j])+scale_dydt_*abs(a[j])*h_);	
				}
				ODE_solver<T>::RKCKOneStep(t,ytemp,yerr,h_);
				//getting the accuracy and ytemp
				errmax  = 0.0; 
				for (j=0; j < dim_; j++){
					errmax = USs::Max<double>(errmax, abs(yerr[j]/derr[j])); // need to scale the error
				}
				//If the observed error yerr exceeds the desired error derr by more than 10% for any component 
				//then the method reduces the step-size by an appropriate factor
				if (errmax <= 1.1) break; // step succeeded. Jump to calculation of the next step
				
				// Reduction of the stepper which goes as much as 10
				htemp = SAFTY*h_*pow(errmax,PSHRNK);
				h_=USs::Max<double>(htemp,0.1*h_);
				if (h_ < 0) UFs::MyError("Routine ODE_solver::RKCKAdapStep: h went below zero");
				if (t+h_ == t) UFs::MyError("Routine ODE_solver::RKCKAdapStep: stepsize underflow");
			}
		
			//Apply the step
			ttemp = t+h_;
			if (ttemp < tstop_minus_tiny){
				for (j=0; j < dim_; j++){
					y_[j] = ytemp[j];
				}
				t = ttemp;
			}
			else {
				hlast = (tstop - t);
				pointerequation_->derivs(t,y_,a);
				ODE_solver<T>::RKCKOneStep(t,ytemp,yerr,hlast);
				for (j=0; j < dim_; j++){
					y_[j] = ytemp[j];
				}
				t = tstop;
			}		

			//Calculation of the next step
			//If the observed error yerr  is less than 50% of the desired error derr for the maximum ratio errmax 
			//then the algorithm takes the opportunity to increase the step-size to bring the error in 
			//line with the desired level,
			if (errmax <= 0.5){
				htemp = SAFTY*h_*pow(errmax,PGROW);
				//Maximum growth is 5
				h_=USs::Min<double>(htemp,5.0*h_);
				if (h_ < 0)  UFs::MyError("Routine ODE_solver::RKCKAdapStep: h went below zero");
			}
		}
		ODE_solver<T>::SaveData();
		//	cout << "h is " << h_;
	}
}


#endif /* OrdinaryDiffEqs_h */

