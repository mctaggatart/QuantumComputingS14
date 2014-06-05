/*
Name: RobustGrape.h
Author: Jay M. Gambetta and Felix Motzoi, Botan Khani

Dependences: MatrixExponential.hpp, QuantumOperations.hpp
Brief Discription: This program implements grape using my matrix class
Limitations: I dont see the need for both statetransfer and unitarytransfer, these should be combined into one.

Version History
	v0: June  9, 2008.
	v1: Feburary 5th, 2009 uses my new matrix class
	v2: Feburary 25th, 2009 edits to make adding different performance functions easier. It turns out it is a bit slower with the pointer to the function but not much
	v3: March ~1, 2009 added envelope for adding nonlinearity to the controls/time e.g. working in the lab frame
	v4: May 20, 2009 added some new functions
	v5: July 31, 2009 added functionality for optical lattices

*/
/*Copy over:	GetFidelity
				GetCount
				StateTransfer
				GetPopulations
*/				
#ifndef RobustGrape_h
#define RobustGrape_h
#include "MatrixExponential.hpp"
#include "QuantumOperations.hpp"
#include "StatisticalFunctions.hpp"
class RobustGrape {
	public:
		RobustGrape(size_t dim, size_t num_controls, size_t num_dof, size_t num_time);// number of controls, the number of time points
		virtual ~RobustGrape();
	
		void SetNumericalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, size_t* nsubpixels);
		void SetHdrift(const matrix<std::complex<double> >& H_drift, const size_t point);
		void SetRhoDesired(const matrix<std::complex<double> >& rho_desired);
		void SetRhoInitial(const matrix<std::complex<double> >& rho_initial);
		void SetHcontrol(const matrix<std::complex<double> >& H_control, const size_t point, const size_t k);
		void Setucontrol(const std::vector<double> u_control, const size_t k);
		void Normalizeucontrol(const size_t k, double value);
		
		matrix<std::complex<double> > GetHdrift(const size_t point);
		matrix<std::complex<double> > GetRhoDesired();
		matrix<std::complex<double> > GetRhoInitial();
		matrix<std::complex<double> > GetHcontrol(const size_t point, const size_t k);
		matrix<std::complex<double> > GetPopulations (const size_t point, size_t initial_level);		
		double* Getucontrol(const size_t k);
		double GetFidelity(const size_t point, double (RobustGrape::* phi)(const size_t) const);
		size_t GetCount();
		
		void StateTransfer(double (RobustGrape::* phi)(const size_t point) const, double (RobustGrape::* gradphi)(const size_t point, const size_t j, const size_t k) const, 
				void (RobustGrape::* ptrSetHnonlinearity)(size_t j), void (RobustGrape::* ptrSetHgradnonlin)(const size_t j));
		//void StateTransfer(double (RobustGrape::* phi)(const size_t point) const, double (RobustGrape::* gradphi)(const size_t j, const size_t k) const)
		//{   StateTransfer(phi, gradphi, &RobustGrape::SetLinear, &RobustGrape::SetLinear);
		//}
/*		void UnitaryTransfer(double (RobustGrape::* phi)() const, double (RobustGrape::* gradphi)(const size_t j, const size_t k) const, 
				void (RobustGrape::* ptrSetHnonlinearity)(size_t j), void (RobustGrape::* ptrSetHgradnonlin)(const size_t j));
		void UnitaryTransfer(double (RobustGrape::* phi)() const, double (RobustGrape::* gradphi)(const size_t j, const size_t k) const)
		{   UnitaryTransfer(phi, gradphi, &RobustGrape::SetLinear, &RobustGrape::SetLinear);
		}*/
		
		double Phi0(const size_t point) const;
		double GradPhi0(const size_t point, const size_t j,const size_t k) const;
		double Phi3(const size_t point) const;
		double GradPhi3(const size_t point, const size_t j, const size_t k) const;
		double Phi4(const size_t point) const;
		double GradPhi4(const size_t point, const size_t j, const size_t k) const;
		double Phi4Sub2(const size_t point) const;
		double GradPhi4Sub2(const size_t point, const size_t j, const size_t k) const;
		void SetLinear(const size_t point, size_t j){}// this is an empty function which tells the the 
		void SetLab2Quadratures(const size_t point, size_t j);
		
	private:
		size_t dim_;//dimensions of the Hilbert space
		double one_on_dim_;//used alot so I have defined it
		matrix<std::complex<double> >* H_drift_; //drift hamiltonian
		matrix<std::complex<double> > rho_desired_; //the disered rho or unitary
		matrix<std::complex<double> > rho_initial_; //the initial rho or unitary
		
		size_t num_controls_; //number of controls,
		size_t num_dof_; //number of system Hamiltonians to control
		matrix<std::complex<double> >** H_controls_; // list of the control hamiltonians
		double **controls_; // list of vectors of u[i][j] where i counts control i and j counts time
		size_t num_time_;//  the number of time points	
		size_t max_iter_; //max iteration befor it turns off
		size_t count_; // counter used to count calls to grape
		size_t failcount_; // counter used to count failed calls to grape		
		double tolerance_; // the minimal change in fidelity 
		double fidelity_; // the desired in fidelity 
		double base_a_; // scale parameter for derivative, epsilon =base_a^power
		double h_; //the time step
		int power_;//scales parameter for derivative, base_a^power
		double epsilon_; //scales parameter for derivative, base_a^power
		double alpha_;
		size_t *nsubpixels_; // how many subpixels to calculate under the (controlable) envelope
		size_t pos_count_; // counts the number of forward steps in the algorithm
		// double (Grape::* gradphiHelper_)(const size_t j, const size_t k) const;
				
		size_t *controlsetflag_;//a array which when all 1 all controls and drift are set
		double ***gradient_;//an array containing the gradiant for each i and j
		matrix<std::complex<double> > Htemp_; //a tempory matrix to store the total hamiltonian
		matrix<std::complex<double> > **H_controls_tdep_; //time-dependent form of the hamiltonian
		std::complex<double> ***Udiag_;
		matrix<std::complex<double> >** rho_; //an array of the forward  states (rho in paper)
		matrix<std::complex<double> >** lambda_;//an array of the backward  states (lambda in paper)
		matrix<std::complex<double> >** Unitaries_;//an array of the unitaries to be applied
		std::complex<double> i; // the complex number i
		
};
inline RobustGrape::RobustGrape(size_t dim, size_t num_controls, size_t num_dof, size_t num_time) : dim_(dim), num_controls_(num_controls), num_dof_(num_dof), num_time_(num_time){
	gradient_ = new double **[num_dof_];
	H_drift_= new  matrix<std::complex<double> > [num_dof_];
	H_controls_=new  matrix<std::complex<double> > *[num_dof_];
	H_controls_tdep_=new  matrix<std::complex<double> > *[num_dof_];
	
	for (size_t point=0; point < num_dof_; ++point){
		gradient_[point] = new double *[num_controls_];
		H_drift_[point].initialize(dim_,dim_);
		H_controls_[point]=new  matrix<std::complex<double> > [num_controls_];
		H_controls_tdep_[point]=new  matrix<std::complex<double> > [num_controls_];
		
		for(size_t k=0; k < num_controls_; ++k){
			H_controls_[point][k].initialize(dim_,dim_);
			H_controls_tdep_[point][k].initialize(dim_,dim_);
			gradient_[point][k] = new double [num_time_];
		}
	}
	
	controls_ = new double *[num_controls_];
	controlsetflag_=new size_t[num_controls_+3]; // one more for Hdrift and one more for rho_desired and rho_initial
	for(size_t k=0; k < num_controls_; ++k){
		controls_[k] = new double[num_time_];
		controlsetflag_[k]=1;//set for the Hcontrols
	}
	controlsetflag_[num_controls_]=1; //set for the Hdrift
	controlsetflag_[num_controls_+1]=1; //set for the rho_desired
	controlsetflag_[num_controls_+2]=1; //set for the rho_initial
	
	rho_      =new  matrix<std::complex<double> > *[num_dof_];
	lambda_   =new  matrix<std::complex<double> > *[num_dof_];
	Unitaries_=new  matrix<std::complex<double> > *[num_dof_];

	for (size_t point=0; point < num_dof_; ++point){
		rho_[point]      =new  matrix<std::complex<double> >[num_time_];
		lambda_[point]   =new  matrix<std::complex<double> >[num_time_];
		Unitaries_[point]=new  matrix<std::complex<double> >[num_time_];

		for(size_t j=0; j < num_time_; ++j){
			rho_[point][j].initialize(dim_,dim_);
			lambda_[point][j].initialize(dim_,dim_);
			Unitaries_[point][j].initialize(dim_,dim_);
		}
	}
	nsubpixels_= new size_t[num_controls];
 	
	i=std::complex<double>(0.0,1.0);
	one_on_dim_=1/double(dim_);
	
	//Is this necessary?
	Udiag_ = new std::complex<double> **[num_controls_];
	for(size_t k=0; k < num_controls_; ++k){
		Udiag_[k] = new std::complex<double> *[2];
		Udiag_[k][0] = new std::complex<double>[dim_];
		Udiag_[k][1] = new std::complex<double>[dim_];
	}
	
}
inline RobustGrape::~RobustGrape(){
	for( size_t k = 0; k < num_controls_; ++k){
		delete [] controls_[k];
		delete [] gradient_[k];
		delete [] Udiag_[k][0];
		delete [] Udiag_[k][1];
		delete [] Udiag_[k];
	}
	delete [] Udiag_;
	delete [] controls_;
	delete [] gradient_;
	delete [] controlsetflag_;
	delete [] H_controls_;
	delete [] rho_;
	delete [] lambda_;
	delete [] nsubpixels_;
	delete [] Unitaries_;
}
inline void RobustGrape::SetNumericalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, size_t* nsubpixels){
	
	tolerance_=tolerance;
	fidelity_=fidelity;
	h_ = h;
	max_iter_ = max_iter;
	base_a_ = base_a;
	epsilon_=epsilon*h;
	for(int k=0; k<num_controls_; k++)
		nsubpixels_[k]=nsubpixels[k];
	alpha_ = 0;//alpha*2.0*h;
	
	
	if(verbose==yes)
	{
	 	std::cout << "--------------------Numerical Parameters--------------------" << std::endl;
		std::cout << "smallest change in phi allowed: " << tolerance_ << std::endl;
		std::cout << "desired fidelity: " << fidelity_ << std::endl;
		std::cout << "step size: " << h_ << std::endl;
		std::cout << "max iterations allow before termination: " << max_iter_ << std::endl;
		std::cout << "number of subpixels in envelope: " << nsubpixels_ << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}
//Set the drift Hamiltonian
inline void RobustGrape::SetHdrift(const matrix<std::complex<double> >& H_drift, size_t point){
	H_drift_[point]=H_drift;
	controlsetflag_[num_controls_]=0;
	if(verbose==yes)
	{
		H_drift_[point].SetOutputStyle(Matrix);
	 	std::cout << "--------------------Hdrift is set---------------------------" << std::endl;
		std::cout << H_drift_[point] << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}
//Set the target state for transfer
inline void RobustGrape::SetRhoDesired(const matrix<std::complex<double> >& rho_desired){
	rho_desired_=rho_desired;
	controlsetflag_[num_controls_+1]=0;
	if(verbose==yes)
	{
		rho_desired_.SetOutputStyle(Matrix);
	 	std::cout << "--------------------rho desired is set---------------------------" << std::endl;
		std::cout << rho_desired_ << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}
//Set the initial state for transfer
inline void RobustGrape::SetRhoInitial(const matrix<std::complex<double> >& rho_initial){
	rho_initial_=rho_initial;
	controlsetflag_[num_controls_+2]=0;
	if(verbose==yes)
	{
		rho_initial_.SetOutputStyle(Matrix);
	 	std::cout << "--------------------rho initial is set---------------------------" << std::endl;
		std::cout << rho_initial_ << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}
//Set the control Hamiltonians
inline void RobustGrape::SetHcontrol(const matrix<std::complex<double> >& H_control, const size_t point, const size_t k){
	H_controls_[point][k]=H_control;
	H_controls_tdep_[point][k]=H_control;
	controlsetflag_[k]=0;
	if(verbose==yes)
	{	
		H_controls_[point][k].SetOutputStyle(Matrix);
	 	std::cout << "------------------The " << k << "th Hcontrol is set-------------------" << std::endl;
		std::cout << H_controls_[point][k] << std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}
//Set the controls
inline void RobustGrape::Setucontrol(std::vector<double> u_control, const size_t k){
	
	for(size_t j = 0; j < num_time_; ++j){
		controls_[k][j]=u_control[j];
	}
	if(verbose==yes)
	{
	 	std::cout << "------------------The " << k << "th ucontrol is set-------------------" << std::endl;
		for(size_t j = 0; j < num_time_; ++j){
			std::cout <<  controls_[k][j] << '\t';
		}
		std::cout<< std::endl;
		std::cout << "------------------------------------------------------------" << std::endl;
	}
}
inline void RobustGrape::Normalizeucontrol(const size_t k, double value)
{
	double sum=0;
	for(size_t j = 0; j < num_time_; ++j){
		sum+=controls_[k][j];
	}
	sum=value/sum/h_;
	std::cout<<sum<<std::endl;
	for(size_t j = 0; j < num_time_; ++j){
		controls_[k][j]*=sum;
	}
}

//Return system values
inline matrix<std::complex<double> > RobustGrape::GetHdrift(const size_t point){
	return H_drift_[point];
}
inline matrix<std::complex<double> > RobustGrape::GetRhoDesired(){
	return rho_desired_;
}
inline matrix<std::complex<double> > RobustGrape::GetRhoInitial(){
	return rho_initial_;
}
inline matrix<std::complex<double> > RobustGrape::GetHcontrol(const size_t point, const size_t k){
	return H_controls_[point][k];
}
inline double* RobustGrape::Getucontrol(const size_t k){
	return controls_[k];
}
//Returns the number of optimization steps in GRAPE
inline size_t RobustGrape::GetCount(){
	return count_;
}
//Returns a full matrix containing the density matrix at the final time
inline matrix<std::complex<double> > RobustGrape::GetPopulations (const size_t point, size_t initial_level){
	//the measure implemented, phi (PHI_0) is phi = trace[rho_desired* UN-1....U_0 rho_initial U_0\dg ... U_N-1\dg]
	char outfile[50];
	//Some flags to make sure all parameters are set
	size_t test=0;
	for(size_t k = 0; k < num_controls_; ++k)
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
		
	// the propagators and the foraward evolution
	for(size_t j = 0; j < num_time_; ++j){
		Htemp_ = H_drift_[point];
		for(size_t k = 0; k < num_controls_; ++k){
			Htemp_ += controls_[k][j]*H_controls_[point][k];
		}
		// std::cout << Htemp_ << std::endl;
		Unitaries_[point][j]=ExpM::EigenMethod(Htemp_,-i*h_);
		if(j==0){
			rho_[point][0] = Unitaries_[point][0]*rho_initial_*MOs::Dagger(Unitaries_[point][0]);
		}
		else{
			rho_[point][j] = Unitaries_[point][j]*rho_[point][j-1]*MOs::Dagger(Unitaries_[point][j]);
		}
	}
		
	//Initialize the matrices used to compute the populations
	matrix<std::complex<double> > populations;
	populations.initialize(num_time_,dim_);

	matrix<std::complex<double> > level, propagator, currentLevel, level2, step, tempMatrix, tempMatrix2;
	level.initialize(dim_,dim_);
	MOs::Null(level);
	
	propagator.initialize(dim_,dim_);
	MOs::Null(propagator);
	
	currentLevel.initialize(dim_, 1);
	MOs::Null(currentLevel);
	
	level2.initialize(dim_, 1);
	MOs::Null(level2);
	
	step.initialize(dim_, 1);
	MOs::Null(step);
	
	tempMatrix.initialize(1, 1);
	MOs::Null(tempMatrix);
	
	tempMatrix2.initialize(1, dim_);
	MOs::Null(tempMatrix2);

	
	for (int i=0; i < dim_; i++)
	{
		level(i, dim_-i-1)=1;
	}
	//Create the initial state selection vector 
	level2(dim_ - initial_level -1, 0)=1;
	
	for (int k=0; k < num_time_; k++)	//Go through each time step
	{
		step=propagator*level2;
		
		for (int j=0; j < dim_; j++)	//Go through each quantum level
		{
			MOs::Null(currentLevel);		
			currentLevel(dim_ -j -1, 0)=1;
			tempMatrix2=MOs::Dagger(currentLevel * step);
			tempMatrix=MOs::Dagger(currentLevel * step)*(currentLevel * step);
			populations(k,j)=tempMatrix(0,0);	//Time=row, Level=col			//This is the occupancy of each level
			//populations[j](0,k)=tempMatrix(0,0);	//Time=Columns, Level=row 
		}
	
		if (k!=num_time_) propagator = rho_[point][k] * propagator;
		if (verbose==yes)
		{
			std::cout << "--------------------Propogator------------------------------" << std::endl;
			USs::OutputMatrix(propagator);
			std::cout << "------------------------------------------------------------" << std::endl;
		}
	} 
	
	return populations;
	//Output the populations to file
//	string line;
//	ofstream moredataout (outfile, ios::app);
//	for(size_t j = 0; j < num_time_; ++j)	//Go through each time step
//	{
//		line=j*h_;
//		for (size_t i=0; i < dim_; i++)			//Go through each dimension of the Hilbert space
//		{
//			line=line "\t" + populations[i](0,j);
//		}
//		//Output line to file
//		moredataout << line << '\n';
//	}
	
//	dataout.close();
}
//Updates the Unitaries and returns the fidelity
inline double RobustGrape::GetFidelity(const size_t point, double (RobustGrape::* phi)(const size_t) const){
	//the measure implemented, phi (PHI_0) is phi = trace[rho_desired* UN-1....U_0 rho_initial U_0\dg ... U_N-1\dg]
	
	//Some flags to make sure all parameters are set
	size_t test=0;
	for(size_t k = 0; k < num_controls_; ++k)
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
		
	// the propagators and the foraward evolution
	for(size_t j = 0; j < num_time_; ++j){
		Htemp_ = H_drift_[point];
		for(size_t k = 0; k < num_controls_; ++k){
			Htemp_ += controls_[k][j]*H_controls_[point][k];
		}
		// std::cout << Htemp_ << std::endl;
		Unitaries_[point][j]=ExpM::EigenMethod(Htemp_,-i*h_);
		if(j==0){
			rho_[point][0] = Unitaries_[point][0]*rho_initial_*MOs::Dagger(Unitaries_[point][0]);
		}
		else{
			rho_[point][j] = Unitaries_[point][j]*rho_[point][j-1]*MOs::Dagger(Unitaries_[point][j]);
		}
	}
		
	//the fidelities and break if we reach the required fidelity
	current_fidelity=(this->*phi)(point); //Phi(rho_desired_,rho_[num_time_-1]);
	return current_fidelity;
}
//Performs the state transfer with given Hamiltonian(s)
void RobustGrape::StateTransfer(double (RobustGrape::* phi)(const size_t) const, double (RobustGrape::* gradphi)(const size_t point, const size_t j, const size_t k) const, 
   void (RobustGrape::* ptrSetHnonlinearity)(size_t j) , void (RobustGrape::* ptrSetHgradnonlin)(const size_t j)){
	

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
	double *current_fidelities;
	double *prev_fidelities;
	current_fidelities=new double[num_dof_];
	prev_fidelities=new double[num_dof_];
	for (size_t j=0; j<num_dof_; ++j)
	{
		prev_fidelities[j]=-1.0;
	}
	double mean_fidelity=0, sum_error=0, STD=0;
	double current_fidelity=0, current_delta_fidelity=1.0, avg_current_delta_fidelity=1.0, last_fidelity=-1.0;
	
	//Starting the algorithim	
	while(avg_current_delta_fidelity*avg_current_delta_fidelity > tolerance_){

		if(count_> max_iter_){	
			std::cout << "Grape::StateTransfer(), algorithim did not converge in the allow max iterations" << std::endl;
			break;
		}
		count_ = count_ + 1;
		mean_fidelity=0.0;
		
		for (size_t point=0; point < num_dof_; ++point)
		{
			// the propagators and the foraward evolution
			for(size_t j = 0; j < num_time_; ++j){
								
				Htemp_ = H_drift_[point];
				for(size_t k = 0; k < num_controls_; ++k){
					Htemp_ += controls_[k][j]*H_controls_[point][k];
				}
				// std::cout << Htemp_ << std::endl;
				Unitaries_[point][j]=ExpM::EigenMethod(Htemp_,-i*h_);
				if(j==0){
					rho_[point][0] = Unitaries_[point][0]*rho_initial_*MOs::Dagger(Unitaries_[point][0]);
					//rho_[0] = rhoupdate(Unitaries_[0],rho_initial);
				}
				else{
					rho_[point][j] = Unitaries_[point][j]*rho_[point][j-1]*MOs::Dagger(Unitaries_[point][j]);
				}
			}
		//	QOs::OutputMatrix(rho_[point][num_time_-1]);
			//the fidelities and break if we reach the required fidelity
			current_fidelities[point]=(this->*phi)(point);
		//	std::cout << "\n Current Fidelity: " << current_fidelities[point];
		}
		
		//Calculate the average of the fidelities
		mean_fidelity=STs::Mean(current_fidelities, num_dof_); 
		STD=STs::STD(current_fidelities, num_dof_); 
		sum_error=num_dof_-STs::Sum(current_fidelities, num_dof_);
		if (verbose==yes){std::cout << "\n Mean Fidelity: " << mean_fidelity << "\t STD Fidelity: " << STD;}

		size_t isControlComplete=1,isControlIncomplete=1, isBelowPrecision=1;
		current_delta_fidelity=0;
		avg_current_delta_fidelity=0;
		for (size_t point=0; point < num_dof_; ++point)
		{
			if(current_fidelities[point] < fidelity_){							//Have all points reached desired fidelity?
				isControlComplete=0;
			}
			if(current_fidelities[point] >= prev_fidelities[point]){			//Have any points increased fidelity?
				isControlIncomplete=0;
			}
			
			//if (verbose==yes){std::cout << "\n Delta Fidelity: " << current_fidelities[point] - prev_fidelities[point];}
			
			current_delta_fidelity=current_fidelities[point]-prev_fidelities[point];
			avg_current_delta_fidelity+=current_delta_fidelity/num_dof_;
			if (fabs(current_delta_fidelity) > tolerance_){						//Have precision requirements been met?
				isBelowPrecision=0;
			}
			prev_fidelities[point]=current_fidelities[point];
		}
		
//		last_fidelity=current_fidelity;
		// std::cout << current_fidelity << '\t' <<  current_delta_fidelity << std::endl;

		//the update in the controls			
		if(isControlIncomplete==0){
			++pos_count_;//if we get to correct directions we set power to be bigger
			if(pos_count_ ==2){
				power_++;
			//	power_=0;
				pos_count_ =0;
			}
		}
		else{
			if (power_ > -52) power_--;
			pos_count_ =0;
		}
		
		//Experiment with this, is it better to modify the control for each point below precision or
		//when all points are below precision?
		for (size_t point=0; point < num_dof_; ++point)
		{
			//the update in the controls			
			if(isControlIncomplete==0){
				// the backward evolution
				lambda_[point][num_time_-1] = rho_desired_;
				for(size_t j = num_time_-1; j > 0; --j){
					//Udg.rho.U
					lambda_[point][j-1] = MOs::Dagger(Unitaries_[point][j])*lambda_[point][j]*Unitaries_[point][j];
				}
				for(size_t j = 0; j < num_time_; ++j){
					for(size_t k =0; k < num_controls_; ++k){							//Modified to take the average of the control
						// std::cout << power_ << std::endl;
						gradient_[point][k][j]=(this->*gradphi)(point,j,k);
						controls_[k][j]+=pow(base_a_,power_)*(gradient_[point][k][j]-alpha_*controls_[k][j])*(1-current_fidelities[point])/sum_error;
					}
				}
				
			}
			else{
				
				for(size_t j = 0; j < num_time_; ++j){
					for(size_t k =0; k < num_controls_; ++k){ 	//Modified to take the average of the control
						controls_[k][j]+=pow(base_a_,power_)*(1-base_a_)*(gradient_[point][k][j]-alpha_*controls_[k][j])*(1-current_fidelities[point])/sum_error;
					}
				}
					
			} //end if
		}//end for
		
		if (isControlComplete==1) std::cout << "\n Control complete";
		if (isControlIncomplete==1) std::cout << "\n Fidelity error";
		if (isBelowPrecision==1) {std::cout << "\n Precision error"; break;}
		
		if (count_%100==0){
			std::cout << "step " << count_ << " which has an mean fidelity of " << mean_fidelity << " and a mean delta fidelity of " << avg_current_delta_fidelity << std::endl; 
		}
	}
	std::cout << "step " << count_ << " which has an mean fidelity of " << mean_fidelity << " and a mean delta fidelity of " << avg_current_delta_fidelity << std::endl;  
}
/* Under Construction...
void RobustGrape::UnitaryTransfer(double (RobustGrape::* phi)() const, double (RobustGrape::* gradphi)(const size_t j, const size_t k) const, 
						void (RobustGrape::* ptrSetHnonlinearity)(size_t j), void (RobustGrape::* ptrSetHgradnonlin)(const size_t j)){
	
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
	pos_count_ =0;
	//the power scale for epsilot (epsilon = base_a^power_)
	power_ =0;
	// fidelity ranges from 0 to 1 with 0 be orthogonal and 1 being the same time
	double current_fidelity=0.0, current_delta_fidelity=1.0;
	double top_fidelity=0.0;
	
	//Starting the algorithim	
	while(current_delta_fidelity*current_delta_fidelity > tolerance_){
		
		if(count_> max_iter_){	
			std::cout << "Grape::UnitaryTransfer(), algorithim did not converge in the allow max iterations" << std::endl;
			break;
		}
		count_ = count_ + 1;
		
		// the propagators and the foraward evolution
		for(size_t j = 0; j < num_time_; ++j){
			(this->*ptrSetHnonlinearity)(j);// updates H_controls_tdep_ for time dependent controls where time dependence is define by ptrSetHnonlinearity (hard coded)
			Htemp_ = H_drift_;
			for(size_t k = 0; k < num_controls_; ++k){				
				for(size_t q=0; q< dim_; ++q){
					for(size_t p=0; p<dim_; ++p){
						Htemp_(q,p) += controls_[k][j]*H_controls_tdep_[k](q,p);
					}
				}
			}
			Unitaries_[j]= ExpM::EigenMethod(Htemp_,-i*h_);

			if(j==0){
				rho_[0] = Unitaries_[0];
			}
			else{
				rho_[j] = Unitaries_[j]*rho_[j-1];
			}
		}	
	
		//the fidelities and break if we reach the required fidelity
		current_fidelity=(this->*phi)();
		// Phi(U_desired,rho_[num_time_-1]);
		
		// //** comment
		// 	rho_[num_time_-1].SetOutputStyle(Matrix);
		// 	std::cout << rho_[num_time_-1] << std::endl;
		// 	matrix<std::complex<double> > tempfelix = rho_[num_time_-1]*MOs::Dagger(rho_[num_time_-1]);
		// 	tempfelix.SetOutputStyle(Matrix);
		// 	std::cout << tempfelix << std::endl;
		
		if(current_fidelity >= fidelity_)
		{	
			std::cout << "success\t";
			break;
		}
		current_delta_fidelity=current_fidelity-top_fidelity;
		
		if (count_%100==1){
			std::cout << "step " << count_ << ", " << failcount_ << " failures, which has an fidelity of " << current_fidelity << " and a delta fidelity of " << current_delta_fidelity << std::endl; 
		}		
								
		//the update in the controls
		if(current_delta_fidelity>=0){
			top_fidelity=current_fidelity;
			// the backward evolution
			lambda_[num_time_-1] = rho_desired_;
			for(size_t j = num_time_-1; j > 0; --j){
				//Udg...
				lambda_[j-1] = MOs::Dagger(Unitaries_[j])*lambda_[j];
			}
			/*	for(size_t j = 0; j < num_time_; ++j){
				for(size_t k =0; k < num_controls_; ++k){
					// std::cout << power_ << std::endl;
					gradient_[k][j]=(this->*gradphi)(j,k);
					controls_[k][j]+=pow(base_a_,power_)*gradient_[k][j];//-alpha_*controls_[k][j]);
				}
			}	*/
/*			if(pos_count_ ==2){
				power_++;
				pos_count_ =0;
			}
			double tempgrad, newgradient, newcontrol;
			for(size_t k =0; k < num_controls_; ++k){
				for(size_t j = 0; j < num_time_; j+=nsubpixels_[k]){
					tempgrad=0;
					for(size_t l = 0; l < nsubpixels_[k]; ++l){
						(this->*ptrSetHgradnonlin)(j+l);
						tempgrad+=(this->*gradphi)(j+l,k);
					}
					newgradient=tempgrad;
					newcontrol=controls_[k][j]+pow(base_a_,power_)*newgradient;
					for(size_t l = 0; l < nsubpixels_[k]; ++l)
					{	gradient_[k][j+l]=newgradient;
			//			std::cout << controls_[k][j+l] << "+" << pow(base_a_,power_) << '*'  << newgradient << " ; " ;
						controls_[k][j+l]=newcontrol;  
					}
				}
			} //std::cout << std::endl;//////
			++pos_count_;//if we get to correct directions we set power to be bigger		
		}
		else{
			power_--;
			pos_count_ =0;
			for(size_t j = 0; j < num_time_; ++j){
				for(size_t k =0; k < num_controls_; ++k){
		//			std::cout << controls_[k][j]-pow(base_a_,power_)*(base_a_)*gradient_[k][j] << "+" <<pow(base_a_,power_)*(1) << '*'  << gradient_[k][j] << " ; " ;
					controls_[k][j]+=pow(base_a_,power_)*(1-base_a_)*gradient_[k][j];
					// std::cout << power_ << std::endl;
				}
			}
			failcount_++;
		//	std::cout << std::endl;//////		
		}

		//if(count_>25) break;
	}
	std::cout << "step " << count_ << ", " << failcount_ << " failures, which has an fidelity of " << current_fidelity << " and a delta fidelity of " << current_delta_fidelity << std::endl;  
}
*/
inline double RobustGrape::Phi0(const size_t point) const{
	//the measure implemented, phi (PHI_0) is phi = trace[rho_desired* UN-1....U_0 rho_initial U_0\dg ... U_N-1\dg]
	std::complex<double> temp1=0.0;
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += conj(rho_desired_(p,q))*rho_[point][num_time_-1](p,q);
		}
	}
	return real(temp1);
}
inline double RobustGrape::GradPhi0(const size_t point, const size_t j, const size_t k) const{
	//phi_0
	std::complex<double> temp1=0;
	matrix<std::complex<double> > Rtemp;
	Rtemp=QOs::SuperHami(H_controls_tdep_[point][k],rho_[point][j]);
	//Fill in the upper triangle
	for (size_t p=0; p< dim_; ++p){
		for(size_t q=p+1; q < dim_; ++q){
			Rtemp(p,q)=std::conj(Rtemp(q,p));
		}
	}	
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 +=std::conj(lambda_[point][j](p,q))*Rtemp(p,q);	
		}
	}
	// std::cout << h_*QOs::Expectation(L,-i*H*R+i*R*H)- h_*temp1 << std::endl;
	// std::cout << temp1 << std::endl;
	// 	Rtemp.SetOutputStyle(Matrix);
	// std::cout << Rtemp << std::endl;		
	return epsilon_*real(temp1);
}
inline double RobustGrape::Phi3(const size_t point) const{
	//the measure implemented, phi (PHI_3) is phi = trace[U_desired* UN-1....U_0]/D
	std::complex<double> temp1=0.0;
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += std::conj(rho_desired_(p,q))*rho_[point][num_time_-1](p,q);
		}
	}
	return std::real(temp1)*one_on_dim_;
}
inline double RobustGrape::GradPhi3(const size_t point, const size_t j, const size_t k) const{
	//phi_3
	std::complex<double> temp1=0;
	matrix<std::complex<double> > Rtemp;
	Rtemp=H_controls_tdep_[point][k]*rho_[point][j];
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 +=std::conj(lambda_[point][j](p,q))*Rtemp(p,q);	
		}
	}	
	return epsilon_*std::imag(temp1);
}
inline double RobustGrape::Phi4(const size_t point) const{
	//the measure implemented, phi (PHI_4) is phi = |trace[U_desired* UN-1....U_0 ]|^2/D^2
	std::complex<double> temp1=0.0;
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += std::conj(rho_desired_(p,q))*rho_[point][num_time_-1](p,q);
		}
	}
	return real(temp1*std::conj(temp1))*one_on_dim_*one_on_dim_;
}
inline double RobustGrape::GradPhi4(const size_t point, const size_t j, const size_t k) const{
	//phi_4
	//2.0*h_*imag(QOs::Expectation( MOs::Dagger(L),H*R)*QOs::Expectation(MOs::Dagger(R),L));
	std::complex<double> temp1=0, temp2=0;
	matrix<std::complex<double> > Rtemp;
	Rtemp=H_controls_tdep_[point][k]*rho_[point][j];
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 +=std::conj(lambda_[point][j](p,q))*Rtemp(p,q);
			temp2 += std::conj(rho_[point][j](p,q))*lambda_[point][j](p,q);	
		}
	}
	return 2.0*epsilon_*std::imag(temp1*temp2);//*one_on_dim_*one_on_dim_;
}
inline double RobustGrape::Phi4Sub2(const size_t point) const{
	//the measure implemented, phi (PHI_4_sub) is phi = |sum_{i=0,1} <i|U_desired* UN-1....U_0|i>|^2/2^2
	std::complex<double> temp1=0.0;
	for(size_t q=0; q< 2; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += std::conj(rho_desired_(p,q))*rho_[point][num_time_-1](p,q);
		}
	}
	return std::real(temp1*std::conj(temp1))*0.25;
}
inline double RobustGrape::GradPhi4Sub2(const size_t point, const size_t j, const size_t k) const{
	//Phi_4_subsystem
	//2.0*h_*imag(QOs::Expectation( MOs::Dagger(L),H*R)*QOs::Expectation(MOs::Dagger(R),L));
	std::complex<double> temp1=0, temp2=0;
	matrix<std::complex<double> > Rtemp;
	Rtemp=H_controls_tdep_[point][k]*rho_[point][j];
	for(size_t q=0; q< 2; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 +=std::conj(lambda_[point][j](p,q))*Rtemp(p,q);
			temp2 += std::conj(rho_[point][j](p,q))*lambda_[point][j](p,q);	
		}
	}
	return 2.0*epsilon_*std::imag(temp1*temp2);//*one_on_dim_*one_on_dim_;
}

void RobustGrape::SetLab2Quadratures(const size_t point, size_t j)
{
	std::complex<double> w=H_drift_[point](1,1);
	for(size_t k = 0; k < num_controls_; ++k){	
			H_controls_tdep_[point][k]= H_controls_[point][k]*cos(((double)j)*w*h_ - (k%2)*3.14159265358/2);
	}
}

//pre: gradphiHelper set to appropriate gradphi function
//inline double Grape::SetGradLab2Quadratures(const size_t j, const size_t k) const{
//	std::complex<double> w=H_drift_(1,1);
//	return (this->*gradphiHelper_)(j, k)*cos(((double)j)*w*h_ - (k%2)*3.14159265358/2);
//}/



#endif /* Grape_h */

