/*
Name: Pulses.cpp
Author: Jay Gambetta

Dependences: UsefullFunctions.hpp
Brief Discription: test the pulses defined in the UsefullFunctions.
Limitations: None

Version History
	v0: February 28, 2009.

*/
#include <UsefullFunctions.hpp>
using namespace std;

int main (int argc, char const *argv[])
{
	size_t num_time=1000;
	double tgate=10.0, dt, amp=1.0, sigma = 1.0, t;
	dt=tgate/double(num_time);
	string const outfile = "Pulses.dat";
	ofstream dataout;
	UFs::OpenFile(outfile,dataout, 16);
	for(size_t j =0; j < num_time; j++){
		t =j*dt;
		dataout << t << '\t' << USs::SquarePulse(t,0.0,tgate,amp) << '\t' << USs::TanhPulse(t,0.0,tgate,amp,sigma) << '\t';
		dataout << USs::Gaussian(t, 0.0, tgate, amp,4*sigma) << '\t' << USs::TruncatedGaussian(t,0.0,tgate,amp,4*sigma) << endl;
	}
	dataout.close();
	
	return 0;
}