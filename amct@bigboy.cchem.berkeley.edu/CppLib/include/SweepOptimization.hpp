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
#include "OptimizeEvolution.hpp"				
				
#ifndef SweepOptimization_h
#define SweepOptimization_h

class SweepOptimization : public OptimizeEvolution {
	public:
		SweepOptimization(size_t dim, size_t num_controls, size_t num_time, const double dt, char* filename)
		:OptimizeEvolution( dim, num_controls, num_time, dt, filename)
		{top_fidelity_=0;}
		virtual ~SweepOptimization(){} 

		void sweepinitialconditions();
		void sweeptimes(ptrPropagate pSearchMethod);
		void sweepenergies(const double min, const double max, const matrix<std::complex<double> >& Shift, const matrix<std::complex<double> >& Shift2, 
										matrix<std::complex<double> >& Hsweep, matrix<std::complex<double> >& Hsweep2, ptrPropagate pSearchMethod);
		void sweepfinegraining(ptrPropagate pSearchMethod);
		void sweeptimesandcompare(ptrPropagate pSearchMethod);
		void sweepfinegrainingandcompare(ptrPropagate pSearchMethod);
		void sweepdimandcompare(ptrPropagate pSearchMethod, const matrix<std::complex<double> >& RefMat, const matrix<std::complex<double> >& RefDesired,Control** refcons);

		size_t gatesteptime;				//when sweeping gate times increment gate time by this amount at each iteration
		size_t ngatetimes;							//number of gate times to sweep through
		double startfid;						//starting fidelity before sweep				
};							


void SweepOptimization::sweepinitialconditions(){
		double bestfidel=0, lastfidel=0;
		size_t bestry, bestry2, bestry3;
		Control *bestControls;
		double bestweight=1000, newweight;

		int numtries = 3;		
		for(size_t tries=0; tries<pow((float)numtries,(int)num_controls_); tries+=1)
		{	
			for(int k=0; k<num_controls_; k++)
			{	if(!(tries % (int)(pow((float)numtries,k)) ))
					((AnalyticControl*)theControls_[k])->scaling = ( (tries/(int)pow((float)numtries,k)) % numtries) / (double)numtries;
				cout << k << " scaling at " << ((AnalyticControl*)theControls_[k])->scaling << endl;
				((AnalyticControl*)theControls_[k])->init(); 
			}	
			cout<< "try " << tries << ", timeslices " <<  endl;
			UnitaryTransfer();
			
			//newweight = theControls_[0]->u_[0]*theControls_[0]->u_[0] + theControls_[1]->u_[0]*theControls_[1]->u_[0] + theControls_[2]->u_[0]*theControls_[2]->u_[0];
			if(1-top_fidelity_<(1-bestfidel)) {//*0.85 || ( 0.85*(1-top_fidelity_)<1-bestfidel && newweight<bestweight) ){ 
					bestfidel = top_fidelity_; bestry = tries; 
					bestweight = newweight;
					for(int k=0; k<num_controls_; k++)
						theControls_[k]->setasbest();
					writecontrols("controls_best");	
					writepopulations("popul_best");
			}
			if(bestfidel>fidelity_)// && newweight>bestweight && tries>pow((float)numtries,(int)num_controls_)/1.8) 
				break;
			lastfidel=top_fidelity_;
		}
		
		cout << "best for " << tgate_ <<": try "<<bestry <<" err " << 1-bestfidel << "\n";	
		top_fidelity_ = bestfidel; //warning: controls not tracked
}

void SweepOptimization::sweeptimes(ptrPropagate pSearchMethod=NULL){ 
	
	double rab=M_PI;
	size_t num_time=theControls_[0]->npixels_;
	ofstream dataout, ctrlout;
	char tempfilename[80];
	
	size_t nsub = theControls_[0]->npixels_/tgate_;
	//size_t nsub = theControls_[0]->nsubpixels_;
	
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_sweeptime.dat"),dataout, 16);
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_timectrl.dat"),ctrlout, 16);
		
	for(size_t itime=0; itime<=10; itime++, num_time+= nsub )
	{	
			SetNumTimes(num_time);
			cout << "sweeping time"<< tgate_<<  ", number of timeslices: " << num_time << endl;
			
			for(int k=0; k<num_controls_; k++)
				((AnalyticControl*)theControls_[k])->init(); 
				
		if(pSearchMethod!=NULL) (this->*pSearchMethod)();
		else{ (this->*pPropagate)();
			top_fidelity_=(this->*Phi)();
		}
		dataout <<  tgate_ <<"\t" << 1-top_fidelity_ << "\n";
		ctrlout << tgate_ <<"\t";
		for(int k=0; k<num_controls_; k++)
			ctrlout << theControls_[k]->u_best[num_time/2] << '\t';  //maybe average is better.
		ctrlout << endl;
		cout <<  "time " << tgate_ <<", error \t" << 1-top_fidelity_ << "\n";
	}
	dataout.close();
	ctrlout.close();
}

void SweepOptimization::sweepenergies(const double min, const double max, const matrix<std::complex<double> >& Shift, const matrix<std::complex<double> >& Shift2, 
										matrix<std::complex<double> >& Hsweep, matrix<std::complex<double> >& Hsweep2, ptrPropagate pSearchMethod=NULL)
{  //datafile
	
	double rab=M_PI;
	size_t num_time=theControls_[0]->npixels_;
	ofstream dataout, ctrlout;
	char tempfilename[80];
		
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_sweepnrg.dat"),dataout, 16);
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_nrgctrl.dat"),ctrlout, 16);
	
	size_t rwa_num_time=theControls_[0]->npixels_;
	matrix<std::complex<double> > H_bckup (Hsweep);
	matrix<std::complex<double> > H_bckup2 (Hsweep2);
	double deltabck = ((AnalyticControl*)theControls_[1])->drivefreq_;
	
	for(double shift=min; shift<max; shift += (max-min)/20 )
	{	
		Hsweep = H_bckup + shift*Shift;
		Hsweep2 = H_bckup2 + shift*Shift2;
		theControls_[1]->drivefreq_ = deltabck*(1+shift);
		
		SetNumTimes(num_time);
		cout << "sweeping anharmonicity"<< theControls_[1]->drivefreq_<<  ", number of timeslices: " << num_time << endl;
			
		for(int k=0; k<num_controls_; k++)
			((AnalyticControl*)theControls_[k])->init(); 
						
		if(pSearchMethod!=NULL) (this->*pSearchMethod)();
		else{ (this->*pPropagate)();
			top_fidelity_=(this->*Phi)();
		}
	
		dataout <<  deltabck*(1+shift) <<"\t" << 1-top_fidelity_ << "\n";
		ctrlout << deltabck*(1+shift) <<"\t";
		for(int k=0; k<num_controls_; k++)
			ctrlout << theControls_[k]->u_best[num_time/2] << '\t';  //maybe average is better.
		ctrlout << endl;
		cout <<  "time " << deltabck*(1+shift) <<", error \t" << 1-top_fidelity_ << "\n";
	}
	
	dataout.close();
	ctrlout.close();
}


inline void SweepOptimization::sweepfinegraining(ptrPropagate pSearchMethod=NULL){ 
	
	double rab=M_PI;
	size_t nsubpix=theControls_[0]->nsubpixels_;
	ofstream dataout;
	char tempfilename[80];
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_sweepgrain.dat"),dataout, 16);
	for(size_t itime=0; itime<=50; itime++, nsubpix+= 1)
	{	
			cout << "sweeping number of sub-timeslices: " << nsubpix << endl;
			h_*=((double)theControls_[0]->nsubpixels_)/nsubpix;
			SetNumTimes(tgate_/h_);
			
			for(int k=0; k<num_controls_; k++)
			{	theControls_[k]->nsubpixels_=nsubpix;
				((AnalyticControl*)theControls_[k])->init(); 
			}
		
		if(pSearchMethod!=NULL) (this->*pSearchMethod)();
		else{ (this->*pPropagate)();
			top_fidelity_=(this->*Phi)();
		}
	
		dataout <<  tgate_ <<"\t" << 1-top_fidelity_ << "\n";
		cout <<  "time " << tgate_ <<", error \t" << 1-top_fidelity_ << "\n";
	}
		dataout.close();
}

inline void SweepOptimization::sweeptimesandcompare(ptrPropagate pSearchMethod=NULL) 
{	
	srand ( time(NULL) );
	ofstream dataout, dataout2, timeout;
	char tempfilename[80];
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_sweeptimes.dat"),dataout, 16);
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_sweepconfigs.dat"),dataout2, 16);
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_runtimes.dat"),timeout, 16);
	size_t num_time=theControls_[0]->npixels_;
	long times0;
	//double avgerr;
	//double div =4;
	//size_t nsub = theControls_[0]->npixels_/tgate_/div;
	size_t nsub = (theControls_[0]->npixels_/tgate_)*gatesteptime;
	//size_t nsub2 = evolutions[0]->theControls_[0]->npixels_/tgate_/div;	
	size_t nsub2 = (evolutions[0]->theControls_[0]->npixels_/tgate_)*gatesteptime;	
	size_t num_time2 = evolutions[0]->theControls_[0]->npixels_;
	double targetfid = fidelity_; 
									
	for(size_t itime=0; itime<ngatetimes; itime++, num_time+= nsub, num_time2+=nsub2)
	{	
			double oldtgate = tgate_;
			SetNumTimes(num_time);
			for(int k=0; k<num_controls_; k++)
			{	//theControls_[k]->nsubpixels_ = nsub;
				((AnalyticControl*)(theControls_[k]))->init(); 	
				theControls_[k]->Interpolate();
			}
			for (size_t s=0; s<num_evol; s++ )
			{
				cout << "evolution"<< s<<"\n";
				evolutions[s]->SetNumTimes((evolutions[s]->theControls_[0]->npixels_/oldtgate)*tgate_);
				//evolutions[s]->SetNumTimes(num_time2);
				for(int k=0; k<num_controls_; k++)
				{	
					evolutions[s]->theControls_[k]->Replicate(theControls_[k]);
					evolutions[s]->theControls_[k]->Interpolate();
				}
				//avgerr=0;
				//evolutions[s]->SetRhoDesired(evolutions->U_[theControls_[0]->npixels_-1]);
				evolutions[s]->SetRhoDesired(rho_desired_);
				evolutions[s]->top_fidelity_=0;
			}	
			cout << "\n\ntime " << num_time_*h_<< " gradient ascent\n";
			fidelity_ = startfid;
			if(pSearchMethod!=NULL) (this->*pSearchMethod)();
			else{   (this->*pPropagate)();
					top_fidelity_=(this->*Phi)();
					cout << "inital error: " << 1 - top_fidelity_ << endl;
			}
			timeout << num_time*h_ <<"\t";
			dataout <<  num_time*h_ <<"\t" <<  1-top_fidelity_ <<"\t";
			dataout2 <<  num_time*h_ <<"\t" <<  1-top_fidelity_ <<"\t";

			for (size_t s=0; s<num_evol; s++ )
			{			
				times0=-clock();
				if(pSearchMethod!=NULL) (evolutions[s]->*pSearchMethod)();
				else{   (evolutions[s]->*pPropagate)();
						evolutions[s]->top_fidelity_=(evolutions[s]->*Phi)();
				}	
				for(int k=0; k<num_controls_; k++)
				{	
					//theControls_[k]->nsubpixels_ = theControls_[k]->npixels_ / evolutions[s]->theControls_[k]->npixels;
					theControls_[k]->Replicate(evolutions[s]->theControls_[k]);					
					theControls_[k]->Interpolate();					
					if(!s)theControls_[k]->writefilewithsuffix("_spli");
					else theControls_[k]->writefilewithsuffix("_pix");					
				}
				(this->*pPropagate)();
				top_fidelity_=(this->*Phi)();
				cout << "actual error: " << 1 - top_fidelity_ << endl;
				times0+=clock();
				//avgerr += 1-evolutions[s]->top_fidelity_;
				timeout << times0/1000.0 << '\t';
				dataout <<   1-top_fidelity_ <<"\t";	
				dataout2 <<  1-evolutions[s]->top_fidelity_ <<"\t";	
			}	
			fidelity_ = targetfid;
			if(pSearchMethod!=NULL) (this->*pSearchMethod)();
			else{   (this->*pPropagate)();
					top_fidelity_=(this->*Phi)();
			}
			cout << "final error: " << 1 - top_fidelity_ << endl;
			dataout <<   1-top_fidelity_ <<"\t";
			timeout << endl;
			dataout << endl;
			dataout2 << endl;

	}
	dataout.close();
	dataout2.close();
	timeout.close();
}

inline void SweepOptimization::sweepfinegrainingandcompare(ptrPropagate pSearchMethod=NULL){
	
	srand ( time(NULL) );
	size_t nsubpix=evolutions[0]->theControls_[0]->nsubpixels_;
	ofstream dataout, dataout2, timeout;
	char tempfilename[80];
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_sweepgrains.dat"),dataout, 16);
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_grainconfigs.dat"),dataout2, 16);
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_grainruntimes.dat"),timeout, 16);
	double nsuperpixels, avgerr;
	cout << "npixels: " << evolutions[0]->theControls_[0]->npixels_ << endl;
	
	nsuperpixels = evolutions[0]->theControls_[0]->npixels_/ evolutions[0]->theControls_[0]->nsubpixels_;
	long times0;
	for(int k=0; k<num_controls_; k++)
				((AnalyticControl*)(theControls_[k]))->init(); 
	double targetfid = fidelity_;
	fidelity_ = startfid;
	UnitaryTransfer();
	fidelity_=targetfid;
	
	double backupcontrols[num_controls_][num_time_];
	
	for(int k=0; k<num_controls_; k++)
		for(int j=0; j<num_time_; j++)	
			backupcontrols[k][j] = theControls_[k]->u_[j];
	
	size_t pixinc = 9;
	//size_t pixinc = 1;
	
	cout << "ngatetimes: " << ngatetimes << endl;
	//for(size_t it=0; it<=40; it++, nsubpix+= pixinc)
	for(size_t it=0; it<=ngatetimes; it++, nsubpix+= pixinc)
	{	
			for(int k=0; k<num_controls_; k++)
				for(int j=0; j<num_time_; j++)	
					 theControls_[k]->u_[j] = backupcontrols[k][j];
			size_t totaln = nsubpix*nsuperpixels;
			
			timeout << nsubpix <<"\t";
			dataout <<  nsubpix <<"\t";
			dataout2 <<  nsubpix <<"\t";
			for (size_t s=0; s<num_evol; s++)
			{
				nsubpix = evolutions[s]->theControls_[0]->nsubpixels_;
				totaln = nsubpix*nsuperpixels;
				evolutions[s]->h_= tgate_ / totaln;
				evolutions[s]->SetNumTimes(totaln);
				for(int k=0; k<num_controls_; k++)
				{	
					//evolutions[s]->theControls_[k]->nsubpixels_+=1;
					evolutions[s]->theControls_[k]->dt_=evolutions[s]->h_;
					evolutions[s]->theControls_[k]->Replicate(theControls_[k]);
					evolutions[s]->theControls_[k]->Interpolate();
				}
				avgerr=0;
				//evolutions[s]->SetRhoDesired(evolutions->U_[theControls_[0]->npixels_-1]);
				evolutions[s]->SetRhoDesired(rho_desired_);
				evolutions[s]->top_fidelity_=0;
			}	
			cout << "\n\npixels " << totaln<< " gradient ascent\n";	
			for (size_t s=0; s<num_evol; s++)
			{	
				times0=-clock();
				if(pSearchMethod!=NULL) (evolutions[s]->*pSearchMethod)();
				else{   (evolutions[s]->*pPropagate)();
						cout << "\n\n1pixels " << totaln<< " gradient ascent\n";
						evolutions[s]->top_fidelity_=(evolutions[s]->*Phi)();
				}	
					
				for(int k=0; k<num_controls_; k++)
				{	
					//theControls_[k]->nsubpixels_ = theControls_[k]->npixels_ / evolutions[s]->theControls_[k]->npixels;
					theControls_[k]->Replicate(evolutions[s]->theControls_[k]);					
					theControls_[k]->Interpolate();					
					if(!s)theControls_[k]->writefilewithsuffix("_fromspli");
					else if(s==2) theControls_[k]->writefilewithsuffix("_fromfilt");		
					else theControls_[k]->writefilewithsuffix("_frompix");
				}			
				
				(this->*pPropagate)();
				top_fidelity_=(this->*Phi)();
				cout << "actual error: " << 1 - top_fidelity_ << endl;
				times0+=clock();
				//avgerr += 1-evolutions[s]->top_fidelity_;
				timeout << times0/1000.0 << '\t';
				dataout2 <<   1-top_fidelity_ <<"\t";	
				dataout <<  1-evolutions[s]->top_fidelity_ <<"\t";	
				for(int k=0; k<num_controls_; k++)
				{	
					evolutions[s]->theControls_[k]->nsubpixels_+=pixinc;
				}
			}	
			timeout << endl;
			dataout << endl;
			
			dataout2 <<   1-top_fidelity_ <<endl;
			if(!it) pixinc++;
	}
	dataout.close();
	timeout.close();
	dataout2.close();
}


inline void SweepOptimization::sweepdimandcompare(ptrPropagate pSearchMethod, const matrix<std::complex<double> >& RefMat,const matrix<std::complex<double> >& RefDesired, Control** refcons){  //ofstream datafile
	
	srand ( time(NULL) );
	size_t nsubpix=evolutions[0]->theControls_[0]->nsubpixels_;
	ofstream dataout, timeout;
	char tempfilename[80];
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_sweepdims.dat"),dataout, 16);
	UFs::OpenFile(strcat(strcpy(tempfilename,filename_), "_dimruntimes.dat"),timeout, 16);
	double avgerr;
	int dim = dim_;
	
	cout << "Ref dim " << RefMat.GetRows() << endl;
	
	for(int k=0; k<num_controls_; k++)
		((AnalyticControl*)(theControls_[k]))->init(); 
	
	
	for(size_t it=0; it<4; it++, dim*=2)
	{	
			cout << "dim " << dim << endl;
			*H_drift_ = RefMat;
			dim_ = dim;
			H_drift_->resize(dim, dim);
			cout << *H_drift_ << endl;
			
			for(int k=0; k<num_controls_; k++)
			{	theControls_[k]->Hcontrol_= refcons[k]->Hcontrol_;
				theControls_[k]->dim_ = dim;	
				theControls_[k]->Hcontrol_.resize(dim,dim);
				theControls_[k]->Z.resize(dim,dim);
				//cout << theControls_[k]->Hcontrol_ << endl;
				delete [] theControls_[k]->W;
				theControls_[k]->W = new double[dim];
				ExpM::EigenMethod(theControls_[k]->Hcontrol_,1.0, &(theControls_[k]->Z), theControls_[k]->W);
			} 
			for(int j=0; j<num_time_; j++)
			{	delete [] W[j];
				W[j] = new double[dim];
				Z[j].resize(dim,dim);
			}
			rho_desired_ = RefDesired;
			rho_desired_.resize(dim,dim);
			rho_desired_.SetOutputStyle(Matrix);
			cout << rho_desired_ << endl;
			
			//Z[0].SetOutputStyle(Matrix);
			//cout << Z[0] << endl;
			for(size_t j=0; j<num_time_; j+=10)
				cout << theControls_[1]->u_[j] << " ";
			cout << "\n";	
			
			top_fidelity_=0.0;
			if(pSearchMethod!=NULL) (this->*pSearchMethod)();
			else{ (this->*pPropagate)();
				top_fidelity_=(this->*Phi)();
			}
			
			//for(size_t j=0; j<num_time_; j+=10)
			//	cout << theControls_[1]->u_[j] << " ";
			//cout << "\n";
			
			timeout << dim <<"\t";
			dataout <<  dim <<"\t";
			continue;
			for (size_t s=0; s<2; s++)
			{
				cout << s << " evolution\n";
				evolutions[s]->dim_= dim;
				evolutions[s]->fidelity_ = 0.00001;
				for(int k=0; k<num_controls_; k++)
				{	evolutions[s]->theControls_[k]->Hcontrol_ = refcons[k]->Hcontrol_;
					evolutions[s]->theControls_[k]->Hcontrol_.resize(dim,dim);
					delete [] evolutions[s]->theControls_[k]->W;
					evolutions[s]->theControls_[k]->W = new double[dim];
					evolutions[s]->theControls_[k]->Z.resize(dim,dim);
					ExpM::EigenMethod(evolutions[s]->theControls_[k]->Hcontrol_,1.0, &(evolutions[s]->theControls_[k]->Z), evolutions[s]->theControls_[k]->W);
				}
				for(int j=0; j<evolutions[s]->num_time_; j++)
				{	delete [] evolutions[s]->W[j];
					evolutions[s]->W[j] = new double[dim];					
					evolutions[s]->Z[j].resize(dim,dim);
				}
				evolutions[s]->SetHdrift(*H_drift_);					
				evolutions[s]->rho_desired_.resize(dim,dim);
				//evolutions[s]->SetRhoDesired(evolutions->U_[theControls_[0]->npixels_-1]);
				evolutions[s]->SetRhoDesired(rho_desired_);
				avgerr=1;
				size_t nsubpix=evolutions[s]->theControls_[0]->nsubpixels_;
				size_t nsuperpixels = evolutions[s]->num_time_/ nsubpix;
				size_t totaln = nsubpix*nsuperpixels;
				while(avgerr>(1-fidelity_)*1.2)
				{
					
					for(int k=0; k<num_controls_; k++)
					{	
						((AnalyticControl*)(evolutions[s]->theControls_[k]))->Replicate(0,0,NULL);
					//	evolutions[s]->theControls_[k]->Interpolate();
					}
					avgerr=0;
					evolutions[s]->top_fidelity_=0;
					evolutions[s]->UnitaryTransfer();
					avgerr += 1-evolutions[s]->top_fidelity_;
					if(avgerr>(1-fidelity_)*1.2)
					{
						nsubpix+=2;
						totaln = nsubpix*nsuperpixels;
						evolutions[s]->h_= tgate_ / totaln;
						evolutions[s]->SetNumTimes(totaln);
						for(int k=0; k<num_controls_; k++)
						{	
							evolutions[s]->theControls_[k]->nsubpixels_=nsubpix;
							evolutions[s]->theControls_[k]->dt_=evolutions[s]->h_;
						}
					}
				}
				evolutions[s]->fidelity_ = 0.999;
				for(int k=0; k<num_controls_; k++)
				{
					if(k==1) ((AnalyticControl*)(evolutions[s]->theControls_[k]))->SquarePulse(freqs_[3]); 
					else ((AnalyticControl*)(evolutions[s]->theControls_[k]))->Null();
					evolutions[s]->theControls_[k]->Interpolate();
				}
				cout << "speed test\n";
				evolutions[s]->times0=evolutions[s]->times1=0;
				evolutions[s]->UnitaryTransfer();
				cout << "tim " << evolutions[s]->times0 << endl;
									
				timeout << evolutions[s]->times0/1000.0 << '\t' << evolutions[s]->times1/1000.0 << '\t';
				dataout <<  evolutions[s]->theControls_[0]->nsubpixels_ <<"\t";	
			}	
			timeout << ((double)evolutions[0]->times1)/evolutions[1]->times1 << endl;
			dataout << endl;			
	}
	dataout.close();
	timeout.close();
}



#endif /* SweepOptimization_h */

