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
#ifndef GrapeUnitary_h
#define GrapeUnitary_h
#include "MatrixExponential.hpp"
#include "QuantumOperations.hpp"
#include "ControlOptimizer.hpp"
class GrapeUnitary;

typedef  double (GrapeUnitary::* Fid)() const;
typedef double (GrapeUnitary::* GradFid)(const size_t j, const size_t k);
typedef  void (GrapeUnitary::* Hnonlin)(size_t j);

class GrapeUnitary : public ControlOptimizer{
	public:
		GrapeUnitary(size_t dim, size_t num_controls, size_t num_time);// number of controls, the number of time points
		virtual ~GrapeUnitary();
	
		void SetPhysicalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, double qubit_freq_, double anharmonicity_);
		
		void UnitaryTransfer(Fid phi, GradFid gradphi, Hnonlin ptrSetHnonlinearity, Hnonlin ptrSetHgradnonlin);
		void UnitaryTransfer(Fid phi, GradFid gradphi)
		{   UnitaryTransfer(phi, gradphi, &GrapeUnitary::SetLinear, &GrapeUnitary::SetGradLinear);
		}
		
		//getshortestcontrols, getcontrols, getrobustcontrols
		void sweeptimes(const char* outfile, Fid phi, GradFid gradphi); //should pass ptr to one-iteration function, then also add "Ascent()" fn
		void sweeptimesandcompare();
		void sweepenergies(const double min, const double max, const matrix<std::complex<double> >& Shift, const matrix<std::complex<double> >& Shift2, matrix<std::complex<double> >& Hsweep, matrix<std::complex<double> >& Hsweep2, const char* outfile, Fid phi, GradFid gradphi); //should pass ptr to one-iteration function, then also add "Ascent()" fn
		void averageconfigurations();
		
		double Phi3() const;
		double GradPhi3(const size_t j, const size_t k) const;
		double Phi4() const;
		double GradPhi4(const size_t j, const size_t k) ;
				
		void SetLinear(size_t j);
		void SetGradLinear(size_t j);
		
		void SetOscillation(size_t j);
		void SetGradOscillation(size_t j);

		size_t nconfigs_;
		Hnonlin ptrSetHnonlinearity; 
		Hnonlin ptrSetHgradnonlin;
		
	protected:
		
		//system specific parameters:
		double qubit_freq_;
		double anharmonicity_;
		double drive_freq_;
		double filter_norm_;
		double filter_rate_;
		double config_;		
		
};

inline GrapeUnitary::GrapeUnitary(size_t dim, size_t num_controls, size_t num_time) : ControlOptimizer(dim, num_controls, num_time){
	
	config_=0;
	drive_freq_=0;
	ptrSetHnonlinearity=&GrapeUnitary::SetLinear; 
	ptrSetHgradnonlin= &GrapeUnitary::SetGradLinear;
}

GrapeUnitary::~GrapeUnitary(){

}

inline void GrapeUnitary::SetPhysicalParameters(const double fidelity, const double h, const double base_a, const double epsilon, const double tolerance, size_t max_iter, 
											double qubit_freq, double anharmonicity)	
{	qubit_freq_=qubit_freq;
	anharmonicity_=anharmonicity;	
	nconfigs_=1;
	SetNumericalParameters(fidelity, h, base_a, epsilon, tolerance, max_iter);
}

void GrapeUnitary::UnitaryTransfer(Fid phi, GradFid gradphi, Hnonlin ptrSetHnonlinearity, Hnonlin ptrSetHgradnonlin){
	static matrix<std::complex<double> > ident(dim_,dim_); 
	MOs::Identity(ident);
	matrix<std::complex<double> >* lastrho = &ident; 
		
	//Some flags to make sure all parameters are set
	size_t test=0;
	for(size_t k = 0; k < num_controls_+2; ++k)
	{
		test+=controlsetflag_[k];
		if(controlsetflag_[k]) cout << "flag " << k << " not set\n";
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
		
		
		
	//	cout << nconfigs_ << "configs\n";
	//	for(phase_=0; phase_<1.99*M_PI; phase_+=2*M_PI/nphases) 
		for(config_=0; config_<nconfigs_; config_++) 
		{	
			// the propagators and the foraward evolution
			lastrho = &ident;
	//		cout << "+";

			for(size_t j = 0; j < num_time_; ++j){			
				(this->*ptrSetHnonlinearity)(j);// updates H_controls_tdep_ for time dependent controls where time dependence is define by ptrSetHnonlinearity (hard coded)
				Htemp_ = *H_drift_;
				for(size_t k = 0; k < num_controls_; ++k){				
					for(size_t q=0; q< dim_; ++q){
						for(size_t p=0; p<dim_; ++p){
							Htemp_(q,p) += H_controls_tdep_[k](q,p);
						}
					}
				}
				
				Unitaries_[j]= ExpM::EigenMethod(Htemp_,-i*h_);	

				rho_[j] = Unitaries_[j]*(*lastrho);
				lastrho=&(rho_[j]);
				//cout << Unitaries_[j] << endl;
			}	
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
	//		cout << "-\n";
			for(size_t k =0; k < num_controls_; ++k){
				for(size_t j = 0; j < num_time_; j+=nsubpixels_[k]){
					tempgrad=0;
					for(size_t l = 0; l < nsubpixels_[k]; ++l){
						(this->*ptrSetHgradnonlin)(j+l);
						tempgrad+=(this->*gradphi)(j+l,k);
					}
					//		std::cout << "wha\n";
					for(size_t l = 0; l < nsubpixels_[k]; ++l)
					{	tempgradient_[k][j+l]=tempgrad;
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
		if (count_%170==1)
		{	std::cout << "step " << count_ << ", " << failcount_ << " failures, with fidelity " << top_fidelity_ << " and delta fidelity " << delta_fidelity*delta_fidelity << " " << tolerance_ << std::endl; 			  
			ofstream dataout;
			UFs::OpenFile("ctrl.dat",dataout, 16);
			
			for(size_t j =0; j < num_time_; j++){
				dataout << h_* j;
				for(int k=0; k<num_controls_; k++)
				 dataout << "\t" << controls_[k][j];
				dataout <<std::endl;
			}
			dataout.close();
		}
   }
	//
	rho_[num_time_-1].SetOutputStyle(Matrix);	
	std::cout << rho_[num_time_-1] << std::endl;
	rho_desired_.SetOutputStyle(Matrix);
	std::cout << rho_desired_ << std::endl;	
	std::cout << "finished on step " << count_ << ", " << failcount_ << " failures, with fidelity " << top_fidelity_ << " and delta fidelity " << delta_fidelity << std::endl;  
	ofstream dataout;
	UFs::OpenFile("ctrl.dat",dataout, 16);
	ofstream popout;
	UFs::OpenFile("pops.dat",popout, 16);
			
	int initpop=1;
	std::complex<double> sign;
	for(size_t j =1; j < num_time_; j++){
				dataout << h_* j;
				for(int k=0; k<num_controls_; k++)
				 dataout << "\t" << controls_[k][j];
				dataout <<std::endl;
		sign = rho_[j](!initpop,initpop)/abs(rho_[j](!initpop,initpop));
		popout << h_* j << "\t" << pow(abs(rho_[j](initpop,initpop)),2) << "\t" << pow(real(rho_[j](initpop,!initpop)/sign),2) << "\t" << pow(abs(imag(rho_[j](initpop,!initpop)/sign)),2) << "\t" << pow(abs(rho_[j](initpop,2)),2) << endl;
		//popout << h_* j;
		//for(int d=0; d<dim_; d++)
		//	popout << '\t' << abs(rho_[j](d,d));
		//popout << endl;
	}
	dataout.close();
	popout.close();
}


//Diagonalize unwanted square driving
/*		matrix<complex<double> > Hd(H_drift_); 
		matrix<complex<double> > Hx(theControls_[0]->Hcontrol_); 
		matrix<complex<double> > Hx2(theControls_[0]->Hcontrol_); 
		
		Hx2 = rab/tgate_ * Hx2;
		Hx = rab/tgate_ * Hx;
		Hx(3,0) = Hx(0,3) = 0;
		Hx(1,4) = Hx(4,1) = 0;		
		Hd = Hd + Hx;
		Hd.SetOutputStyle(Matrix);
								
		double vals[dim_];
		Hx2 = Hx2 - Hx;
		QOs::MatrixElements(Hd, dim_, vals, Hx2);
		QOs::MatrixElements(Hd, dim_, vals, Hd);
		cout << Hx2 << endl;
		
		 Hx = Hx2;		
		 Hx(4,0) = Hx(0,4) = 0;	
			
		 Hd = Hd + Hx;		
		 Hx2 = Hx2 - Hx;
		 //cout << Hx2 << endl;
		 QOs::MatrixElements(Hd, dim_, vals, Hx2);
		 QOs::MatrixElements(Hd, dim_, vals, Hd);
		 //cout << Hx2 << endl;
		
		 
		double rab2 = M_PI/4/real(Hx2(0,4))/tgate_;
		cout << "rab2 " << rab2 << endl;
		
		 Hd = H_drift_;
		 Hx2 = sqrt(abs(rab2)) * rab/tgate_ * theControls_[0]->Hcontrol_;
		 Hx = Hx2;
		Hx(3,0) = Hx(0,3) = 0;
		Hx(1,4) = Hx(4,1) = 0;		
		Hd = Hd + Hx;
		Hx2 = Hx2 - Hx;
		QOs::MatrixElements(Hd, dim_, vals, Hx2);
	//	QOs::MatrixElements(Hd, dim_, vals, rho_desired_);
	//	rho_desired_.SetOutputStyle(Matrix);
	//	cout << rho_desired_ << endl;
		
		QOs::MatrixElements(Hd, dim_, vals, Hd);
		//cout << Hx2 << endl;
		
		 Hx = Hx2;		
		 Hx(4,0) = Hx(0,4) = 0;	
			
		 Hd = Hd + Hx;		
		 Hx2 = Hx2 - Hx;
	//	 cout << Hx2 << endl;
		 QOs::MatrixElements(Hd, dim_, vals, Hx2);
		 QOs::MatrixElements(Hd, dim_, vals, Hd);
		 //cout << Hx2 << endl;
			
		double stk = real(Hd(0,0)-H_drift_(0,0)-Hd(4,4)+H_drift_(4,4));
		
//		rab *= sqrt(abs(rab2*M_PI/4/real(Hx2(0,4))/tgate_));
//		rab_amp = rab*rab/tgate_;
		
		//rab = M_PI/4/real(Hx2(0,4))/tgate_;
//		stk = -real(Hd(0,0)-Hd(4,4));
//		
//		Hx = H_drift_;
//		H_drift_=Hd;
//		Hd = theControls_[0]->Hcontrol_;
//		theControls_[0]->Hcontrol_ = Hx2;
*/



void GrapeUnitary::sweeptimes(const char* outfile, Fid phi, GradFid gradphi){  //ofstream datafile
	
	ofstream dataout;
	UFs::OpenFile(outfile,dataout, 16);	
	
	//harmon
	size_t rwa_num_time=theControls_[0]->npixels_;
	//cout << theControls_[0]->Hcontrol_ << endl;
	
	ofstream datatry;
	UFs::OpenFile("tries.dat",datatry, 16);
	
	
	double lambda, lamb2;
	double rab = M_PI;
	double Delta, Delta2;
	double sysdata[4];

//0-1 driving
		sysdata[0] = Delta = anharmonicity_;
		sysdata[1] = lambda = 2*real(theControls_[0]->Hcontrol_(1,2));
		sysdata[2] = Delta2 = 1; //real((*H_drift_)(0,0));
		sysdata[3] = lamb2 = 0; //2*real(theControls_[0]->Hcontrol_(0,1));

//1-2 driving
	//	sysdata[0] = Delta = anharmonicity_;
//		sysdata[1] = lambda = 2*real(theControls_[0]->Hcontrol_(2,3));
//		sysdata[2] = Delta2 = real((*H_drift_)(0,0));
//		sysdata[3] = lamb2 = 2*real(theControls_[0]->Hcontrol_(0,1));
//	
	cout << "lamb2/lambda " << lamb2/lambda << endl;
	
	cout << "aaaa" << lamb2 << endl;
	for(size_t itime=0; itime<1; itime++, rwa_num_time+=40*theControls_[0]->nsubpixels_)
	{	
		double bestfidel=0, lastfidel=0;
		size_t bestry, bestry2, bestry3;
		double rwa_dt = tgate_/double(rwa_num_time);
		SetNumTimes(rwa_num_time);
		
		//5.0000000000000000e+00	22	89	23
		//9.0000000000000000e+00	35	72	28
		//1.3000000000000000e+01	42	56	24
		//1.7000000000000000e+01	49	46	24
		//2.1000000000000000e+01	54	39	23

//		for(size_t tries=0; tries<60; tries++)
//		{	
//			for(size_t tries2=30; tries2<90; tries2++)
//			{
//			//std::cout << tgate_ << " " << tries<< " " << tries2<< " " <<  << "\n";
//			for(size_t tries3=10; tries3<35; tries3++)
//			{
		for(size_t tries=53; tries<54; tries++)
		{	
			for(size_t tries2=42; tries2<43; tries2++)
			{
			//std::cout << tgate_ << " " << tries<< " " << tries2<< " " <<  << "\n";
			for(size_t tries3=0; tries3<1; tries3++)
			{
					
		//	if(tries==62) iscrap=true;
			std::cout << tgate_ << " " << tries<< " " << tries2<< " " << tries3 << "\n";
			//allocate controls (..dynamically)
			
			Control refcontrol(*(theControls_[0]), theControls_[0]->nsubpixels_);
			refcontrol.ShiftedGaussian( rab, tgate_/1.5, NULL, NULL);
			refcontrol.Normalize(rab);
			refcontrol.writefile("ref.dat");
			
			for(int k=0; k<4; k++)
				theControls_[k]->init(&refcontrol, sysdata); 
			
			theControls_[0]->writefile("gaussiantrig.dat");
			theControls_[1]->writefile("gaussiantrig1.dat");
			theControls_[2]->writefile("gaussiantrig2.dat");
			theControls_[3]->writefile("gaussiantrig3.dat");				

				//ShiftedGaussian(rab, tgate_/1.5, theControls_[0], NULL);
			//theControls_[0]->TanhPulse(0, tgate_,  rab/tgate_, tgate_/2.5);		
			//range: (100-180)/250	
			//theControls_[0]->SquarePulse( 0, tgate_, rab/tgate_);
			//theControls_[0]->writefile("cont0.dat");
			//theControls_[0]->Pixelate();
			//theControls_[0]->writefile("pix0.dat");
			//double norm1 = theControls_[0]->SquareNormalize(rab_amp);
			//double norm1 = theControls_[0]->Normalize(rab);
			//double norm1 = theControls_[0]->SquareNormalize(rab_amp);
			//theControls_[1]->Pixelate();

			//cout << "norm " << norm1 << endl;
			
			//theControls_[0]->writefile("one0.dat");			
							
//			cout << "x = " << theControls_[0]->u_[num_time_/2] << endl;
//			cout << "x3 = " << rab_amp2c << endl; //sqrt(abs( theControls_[0]->u_[num_time_/2]*theControls_[0]->u_[num_time_/2]*theControls_[0]->u_[num_time_/2]*theControls_[0]->u_[num_time_/2]*shift2a))<< endl;
														
			//	theControls_[1]->Derivative(-1/Delta, theControls_[0]);
//			theControls_[1]->Null();
//			theControls_[1]->init(theControls_[0], sysdata); 

			//
			//if(tries==2)
			double integ=0, avg=0;
			for(size_t j = 0; j < num_time_; ++j)
			{	double x = theControls_[0]->u_[j];
				//avg+=h_ * (lambda*lambda/Delta-lamb2*lamb2/Delta2)*(x*x)/4;
				avg+=h_ *  theControls_[2]->u_[j];
			}
			avg/=tgate_;
			for(size_t j = 0; j < num_time_; ++j)
			{	//double x = theControls_[0]->u_[j];
				double x = refcontrol.u_[j];
//				//1ch
//				theControls_[2]->u_[j]=(lambda*lambda-4)*(x*x)/Delta/4;
//								-((27 -16 *lambda*lambda + 2*lambda*lambda*lambda*lambda) *x*x*x*x )/(Delta*Delta*Delta)/16;
//				theControls_[2]->u_[j]-=x*x*x*x*(12-7*lambda*lambda+lambda*lambda*lambda*lambda)/(Delta*Delta*Delta)/16;
//				theControls_[3]->u_[j]=-lambda*theControls_[1]->u_[j]*theControls_[0]->u_[j]/4/Delta;
//				theControls_[0]->u_[j]-=x*x*x*x*x*(112-76*lambda*lambda+13*lambda*lambda*lambda*lambda)/(Delta*Delta*Delta*Delta)/128;
//				theControls_[0]->u_[j]+=(lambda*lambda-4)*(x*x*x)/Delta/Delta/8;
	//			theControls_[0]->u_[j]+=(2-4)*(x*x*x)/2.5/2.5/8;
//						-(x*x*x*x*x* (104 - 84 *lambda*lambda + 13*lambda*lambda*lambda*lambda))/(8*Delta*Delta*Delta*Delta)/16;
//				theControls_[1]->u_[j]*=(1+3*x*x*(65-28*lambda*lambda)/(12*Delta*Delta*Delta)/4);
				//2ch  Adi
				//theControls_[2]->u_[j]=(lambda*lambda/Delta-lamb2*lamb2/Delta2)*(x*x)/4;
	//			theControls_[0]->u_[j]+=(lambda*lambda/Delta/Delta+lamb2*lamb2/Delta2/Delta2)*(x*x*x)/8;
//				integ += h_ * ((lambda*lambda/Delta-lamb2*lamb2/Delta2)*(x*x)/4);
				//integ += h_ * (theControls_[2]->u_[j] - avg);
//				integ += h_ * (theControls_[2]->u_[j]);
//				theControls_[2]->u_[j]=integ/(h_*(j+0.5));//avg;
				//theControls_[0]->u_[j]=x*cos(integ) - theControls_[1]->u_[j]*sin(integ);
				//theControls_[1]->u_[j]=theControls_[1]->u_[j]*cos(integ) + x*sin(integ);
//				theControls_[2]->u_[j]-=(lambda*lambda-0.5)*(lambda*lambda+1.5)*(x*x*x*x)/anharmonicity_/anharmonicity_/anharmonicity_/8;
//				theControls_[2]->u_[j]-=(lambda*lambda-0.5)*theControls_[1]->u_[j]*theControls_[1]->u_[j]/anharmonicity_/anharmonicity_/anharmonicity_;
				//theControls_[1]->u_[j]=theControls_[0]->u_[j]*cos(
//				theControls_[0]->u_[j]-=x*(0.0625*x*x*x*x*(0.625/4 - 0.5*lambda*lambda + 0.625*lambda*lambda*lambda*lambda + 0.5* (-0.5 - 1.25*lambda*lambda)) - 0.5* (0.5 + lambda*lambda)* (theControls_[1]->u_[j]*theControls_[1]->u_[j])/4)/anharmonicity_/anharmonicity_/anharmonicity_/anharmonicity_;
//				theControls_[1]->u_[j]=0;
				
				//2Ch Drg
	//			theControls_[2]->u_[j]=(x*x)*( -lamb2*lamb2/Delta2 +(lambda*lambda - 4)/Delta)/4;
	//			theControls_[0]->u_[j]-=(x*x*x)*( (4-lambda*lambda)/Delta/Delta - lamb2*lamb2/Delta2/Delta2 )/8;
	//			theControls_[0]->u_[j]-=(x*x*x)*( (4-3/2)/2.5/2.5 - 1/2/2.5/2.5 )/8;
				//theControls_[0]->u_[j]-=(x*x*x)*( (4-lambda*lambda)/Delta/Delta + lamb2*lamb2*(1/Delta2+2/Delta)/Delta2)/8;

				//2Ch Drg/2
				//0: theControls_[2]->u_[j]=(x*x)*( -lamb2*lamb2/Delta2 +(lambda*lambda - 4/4)/Delta)/4;
	//			theControls_[0]->u_[j]-=(x*x*x)*( (4/4 -lambda*lambda)/Delta/Delta - lamb2*lamb2/Delta2/Delta2 )/8;
				//theControls_[0]->u_[j]+=(x*x*x)*( lamb2*lamb2/Delta2/Delta2 )/4;
				

			}
			
		theControls_[0]->writefile("one0.dat");
		theControls_[1]->writefile("one1.dat");
		theControls_[2]->writefile("one2.dat");
	//	theControls_[3]->writefile("one3.dat");
			
		//	for(int k=0; k<3; k++)
		//		theControls_[k]->Pixelate();
			//else 
			if(tries>1)
			{
	//			theControls_[0]->RandomControl(0, M_PI/tgate_*2.0);
	//			cout << theControls_[0]->Normalize(rab)  << " norm\n";
	//			if(tries%2)
	//					theControls_[1]->RandomControl(-M_PI/anharmonicity_/tgate_/3, M_PI/anharmonicity_/tgate_/3);
			}
	//		theControls_[0]->writefile("one0.dat");
	//	theControls_[1]->writefile("one1.dat");
	//	theControls_[2]->writefile("one2.dat");
	//cout << "unitary\n";
	
	
		UnitaryTransfer(phi, gradphi, ptrSetHnonlinearity, ptrSetHgradnonlin); 
	
		//	H_drift_ = Hx;
		//	theControls_[0]->Hcontrol_ = Hd;
		//	cout << tgate_ <<"\t" << 1-top_fidelity_ << "\n";
			if(top_fidelity_>bestfidel){ bestfidel = top_fidelity_; bestry = tries; bestry2 = tries2;bestry3 = tries3;}
			//if(top_fidelity_>fidelity_) break;
				//if(top_fidelity_<lastfidel) break;
				lastfidel=top_fidelity_;
			}
			}
		}
		
		//cout << 1-abs(rho_[num_time_-1](1,2)) << endl;;
		dataout <<  tgate_ <<"\t" << 1-bestfidel << "\n";
		datatry << tgate_ <<"\t" <<bestry<<"\t" <<bestry2<<"\t" <<bestry3 << "\n";
		cout << "best for " << tgate_ <<" was " << 1-bestfidel << "\n";
		//if (top_fidelity_==fidelity_)
		

	}
	datatry.close();
	dataout.close();
}

void GrapeUnitary::sweepenergies(const double min, const double max, const matrix<std::complex<double> >& Shift, const matrix<std::complex<double> >& Shift2, matrix<std::complex<double> >& Hsweep, matrix<std::complex<double> >& Hsweep2, const char* outfile, Fid phi, GradFid gradphi)
{  //datafile
	
	ofstream dataout;
	UFs::OpenFile(outfile,dataout, 16);
	
	size_t rwa_num_time=theControls_[0]->npixels_;
	matrix<std::complex<double> > H_bckup (Hsweep);
	matrix<std::complex<double> > H_bckup2 (Hsweep2);
	
	double lambda, lamb2;
	double rab = M_PI;
	double Delta = anharmonicity_;
	double Delta2 = anharmonicity_;
	
		cout << min << " " << max << endl;;
	for(double shift=min; shift<max; shift +=  (max-min)/80 )
	{	
		Hsweep = H_bckup + shift*Shift;
		Hsweep2 = H_bckup2 + shift*Shift2;
	//	if(abs(anharmonicity_)>abs(shift)) {
	//		Delta = shift;
	//		Delta2 = anharmonicity_;
	//	} else {
			Delta = anharmonicity_;
			Delta2 = shift;
	//	}
		
		double sysdata[4];
		sysdata[0] = Delta = anharmonicity_;
		sysdata[1] = lambda = 2*real(theControls_[0]->Hcontrol_(2,3));
		sysdata[2] = Delta2 = real((*H_drift_)(0,0));
		sysdata[3] = lamb2 = 2*real(theControls_[0]->Hcontrol_(0,1));
		

		
		Control refcontrol(*(theControls_[0]), theControls_[0]->nsubpixels_);
			refcontrol.ShiftedGaussian( rab, tgate_/3, theControls_[0], NULL);
			refcontrol.Normalize(rab);
			refcontrol.writefile("ref.dat");
			
		double bestfidel=0;
		//double rwa_dt = tgate_/double(rwa_num_time);
		//SetNumTimes(rwa_num_time);
		for(int k=0; k<3; k++)
				theControls_[k]->init(&refcontrol, sysdata);
			std::cout << tgate_ << " " << lamb2<< "\n";
				
		for(size_t tries=0; tries<1; tries++)
		{	
		//	if(tries==62) iscrap=true;
			std::cout << tgate_ << " " << tries<< "\n";
			
			//allocate controls (..dynamically)
		//	theControls_[0]->ShiftedGaussian(rab, tgate_/3, theControls_[0], NULL);
			//theControls_[0]->SquarePulse( 0, tgate_, M_PI/tgate_);
			//theControls_[0]->writefile("one0.dat");
		//	double norm = theControls_[0]->Normalize(rab);
			//if(tries>0) theControls_[0]->SquarePulse( 0, tgate_, M_PI/tgate_);
		//	double mindel = abs(Delta/lambda/lambda)<abs(Delta2/lamb2/lamb2) ? Delta : -Delta2;
		//	theControls_[1]->Derivative(-1/mindel, 1, theControls_[0], NULL);
		//	theControls_[1]->Derivative(-1/anharmonicity_, 1, theControls_[0], NULL);
		//	theControls_[1]->Derivative(-(lambda*lambda/Delta-lamb2*lamb2/Delta2)/4, 1, theControls_[0], NULL);
			
			//
			//if(tries==2)
			for(size_t j = 0; j < num_time_; ++j)
			{	double x = theControls_[0]->u_[j];
//				//1ch
		//		theControls_[2]->u_[j]=(lambda*lambda-4)*(x*x)/anharmonicity_/4;
		//						-((27 -16 *lambda*lambda + 2*lambda*lambda*lambda*lambda) *x*x*x*x )/(Delta*Delta*Delta)/16;
		//		theControls_[3]->u_[j]=-lambda*theControls_[1]->u_[j]*x/2/Delta;
		//		theControls_[0]->u_[j]+=(lambda*lambda-4)*(x*x*x)/anharmonicity_/anharmonicity_/8;
		//				-(x*x*x*x*x* (104 - 84 *lambda*lambda + 13*lambda*lambda*lambda*lambda))/(8*Delta*Delta*Delta*Delta)/16;
					
		//		theControls_[1]->u_[j]*=(1+3*x*x*(65-28*lambda*lambda)/(12*Delta*Delta*Delta)/4);
				//2ch  Adi
		//		theControls_[2]->u_[j]=(lambda*lambda/Delta-lamb2*lamb2/Delta2)*(x*x)/4;

//				theControls_[2]->u_[j]-=(lambda*lambda-0.5)*(lambda*lambda+1.5)*(x*x*x*x)/anharmonicity_/anharmonicity_/anharmonicity_/8;
//				theControls_[2]->u_[j]-=(lambda*lambda-0.5)*theControls_[1]->u_[j]*theControls_[1]->u_[j]/anharmonicity_/anharmonicity_/anharmonicity_;
		//		theControls_[0]->u_[j]+=(lambda*lambda/Delta/Delta+lamb2*lamb2/Delta2/Delta2)*(x*x*x)/8;

//				theControls_[0]->u_[j]-=x*(0.0625*x*x*x*x*(0.625/4 - 0.5*lambda*lambda + 0.625*lambda*lambda*lambda*lambda + 0.5* (-0.5 - 1.25*lambda*lambda)) - 0.5* (0.5 + lambda*lambda)* (theControls_[1]->u_[j]*theControls_[1]->u_[j])/4)/anharmonicity_/anharmonicity_/anharmonicity_/anharmonicity_;
//				theControls_[1]->u_[j]=0;
				
				//2Ch Drg
		//		theControls_[2]->u_[j]=(x*x)*( -lamb2*lamb2/Delta2 +(lambda*lambda - 4)/Delta)/4;
		//		theControls_[0]->u_[j]-=(x*x*x)*( (4-lambda*lambda)/Delta/Delta - lamb2*lamb2/Delta2/Delta2 )/8;	
	
				//Drg/2
		//		theControls_[0]->u_[j]+=(x*x*x)*( lamb2*lamb2/Delta2/Delta2 )/4;
										
			}
			
			
			//else 
			if(tries>1)
			{
				theControls_[0]->RandomControl(0, M_PI/tgate_*2.0, theControls_[0], NULL );
				
				cout << theControls_[0]->Normalize(rab)  << " norm\n";
				
				if(tries%2)
						theControls_[1]->RandomControl( -M_PI/anharmonicity_/tgate_/3, M_PI/anharmonicity_/tgate_/3, theControls_[0], NULL);
			}
			
			
		theControls_[0]->Hcontrol_.SetOutputStyle(Matrix);
		cout << theControls_[0]->Hcontrol_ << endl;
		H_drift_->SetOutputStyle(Matrix);
		cout << Hsweep << endl;
		cout << *H_drift_ << endl;
		cout << lambda << " " << lamb2<< endl;

			
			UnitaryTransfer(phi, gradphi); 
			cout << tgate_ <<"\t" << 1-top_fidelity_ << "\n";
			if(top_fidelity_>bestfidel) bestfidel = top_fidelity_;
			//if(top_fidelity_>fidelity_) break;
		}
		
		//cout << 1-abs(rho_[num_time_-1](1,2)) << endl;;
	//	dataout <<  shift/anharmonicity_ <<"\t" << 1-bestfidel << "\n";
		dataout <<  shift/lambda <<"\t" << 1-bestfidel << "\n";
	
			cout << "best for " << shift/anharmonicity_ <<" Deltas was " << 1-bestfidel << "\n";
		//if (top_fidelity_==fidelity_)
		theControls_[0]->writefile("one0.dat");
		theControls_[1]->writefile("one1.dat");
		theControls_[2]->writefile("one2.dat");
		}
	dataout.close();
}
		

inline double GrapeUnitary::Phi3() const{
	//the measure implemented, phi (PHI_3) is phi = trace[U_desired* UN-1....U_0]/D
	std::complex<double> temp1=0.0;
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += std::conj(rho_desired_(p,q))*rho_[num_time_-1](p,q);
		}
	}
	return std::real(temp1)*one_on_dim_;
}
inline double GrapeUnitary::GradPhi3(const size_t j, const size_t k) const{
	//phi_3
	std::complex<double> temp1=0;
	static matrix<std::complex<double> > Rtemp;
	Rtemp=H_controls_tdep_[k]*rho_[j];
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 +=std::conj(lambda_[j](p,q))*Rtemp(p,q);	
		}
	}	
	return epsilon_*std::imag(temp1);
}
inline double GrapeUnitary::Phi4() const{
	//the measure implemented, phi (PHI_4) is phi = |trace[U_desired* UN-1....U_0 ]|^2/D^2
	std::complex<double> temp1=0.0;
	for(size_t q=0; q< dim_; ++q){
		for(size_t p=0; p<dim_; ++p){
			temp1 += std::conj(rho_desired_(p,q))*rho_[num_time_-1](p,q);
		}
	}
	return real(temp1*std::conj(temp1))*one_on_dim_*one_on_dim_;
}
inline double GrapeUnitary::GradPhi4(const size_t j, const size_t k) {
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

void GrapeUnitary::SetLinear(size_t j)
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
		H_controls_tdep_[k]=controls_[k][j]*(*H_controls_[k]);
		
	}
//	H_controls_[0](1,2)=H_controls_[0](2,1)=lambda;
//	H_controls_[1](1,2)=lambdai;
//	H_controls_[1](2,1)=-H_controls_[1](1,2);

}
void GrapeUnitary::SetGradLinear(size_t j)
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
		H_controls_tdep_[k]=*H_controls_[k];
	}
//	H_controls_[0](1,2)=H_controls_[0](2,1)=lambda;
//	H_controls_[1](1,2)=lambdai;
//	H_controls_[1](2,1)=-H_controls_[1](1,2);
}

void GrapeUnitary::SetOscillation(size_t j)
{
	double phase = config_*2*M_PI/nconfigs_; 

	for(size_t l=0; l<dim_; l++) {
		Udiag_[0][0][l]=exp(-std::complex<double>(0.0,-1.0*(l%2)*(phase +(double)j*anharmonicity_*h_)));
		Udiag_[0][1][l]=exp(std::complex<double>(0.0, -1.0*(l%2)*(phase +(double)j*anharmonicity_*h_)));
		
		Udiag_[1][0][l]=exp(-std::complex<double>(0.0,-1.0*(l%2)*(phase +(double)j*anharmonicity_*h_)));
		Udiag_[1][1][l]=exp(std::complex<double>(0.0,-1.0*(l%2)*(phase +(double)j*anharmonicity_*h_)));						
		
		Udiag_[2][0][l]=exp(-std::complex<double>(0.0,1.0*(l%2+2*(l>1))*(phase +(double)j*anharmonicity_*h_)));
		Udiag_[2][1][l]=exp(std::complex<double>(0.0,1.0*(l%2+2*(l>1))*(phase +(double)j*anharmonicity_*h_)));		
		
		Udiag_[3][0][l]=exp(-std::complex<double>(0.0,1.0*(l%2+2*(l>1))*(phase +(double)j*anharmonicity_*h_)));
		Udiag_[3][1][l]=exp(std::complex<double>(0.0,1.0*(l%2+2*(l>1))*(phase +(double)j*anharmonicity_*h_)));				
	}		
							
	for(size_t k = 0; k < num_controls_; ++k) {	
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				H_controls_tdep_[k](p,q)= Udiag_[k][0][q]*controls_[k][j]*(*H_controls_[k])(p,q)*Udiag_[k][1][p];
			}
		}
	}
	
}

void GrapeUnitary::SetGradOscillation(size_t j)
{
	double phase = config_*2*M_PI/nconfigs_; 

	for(size_t l=0; l<dim_; l++) {
	Udiag_[0][0][l]=exp(-std::complex<double>(0.0,-1.0*(l%2)*(phase +(double)j*anharmonicity_*h_)));
		Udiag_[0][1][l]=exp(std::complex<double>(0.0, -1.0*(l%2)*(phase +(double)j*anharmonicity_*h_)));
		
		Udiag_[1][0][l]=exp(-std::complex<double>(0.0,-1.0*(l%2)*(phase +(double)j*anharmonicity_*h_)));
		Udiag_[1][1][l]=exp(std::complex<double>(0.0,-1.0*(l%2)*(phase +(double)j*anharmonicity_*h_)));						
		
		Udiag_[2][0][l]=exp(-std::complex<double>(0.0,1.0*(l%2+2*(l>1))*(phase +(double)j*anharmonicity_*h_)));
		Udiag_[2][1][l]=exp(std::complex<double>(0.0,1.0*(l%2+2*(l>1))*(phase +(double)j*anharmonicity_*h_)));		
		
		Udiag_[3][0][l]=exp(-std::complex<double>(0.0,1.0*(l%2+2*(l>1))*(phase +(double)j*anharmonicity_*h_)));
		Udiag_[3][1][l]=exp(std::complex<double>(0.0,1.0*(l%2+2*(l>1))*(phase +(double)j*anharmonicity_*h_)));						
	}		
							
	for(size_t k = 0; k < num_controls_; ++k) {	
		for(size_t q=0; q< dim_; ++q){
			for(size_t p=0; p<dim_; ++p){
				H_controls_tdep_[k](p,q) = Udiag_[k][0][q]*(*H_controls_[k])(p,q)*Udiag_[k][1][p];
			}
		}
	}
	
}


#endif /* GrapeUnitary_h */

