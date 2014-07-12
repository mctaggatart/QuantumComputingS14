/*
Name: TestErrorFunctions.cpp
Author: Jay Gambetta

Dependences: SpecialFunctions.hpp
Brief Discription: This program tests the special functions Erf, Erfc, Erffast, Erfcfast, and implicity LogGamma, 
	GammaP, GammaQConFrac, GammaPSeries
Limitations: None

Version History
	v0: February 12, 2008.

*/ 
#include <SpecialFunctions.hpp>
using namespace std;

int main(void){

	cout.setf(ios::scientific, ios::floatfield); 
	cout.precision(30);


	cout << " -------Error Function erf(z) ------- "<< endl;
	double a, b;
	cout << "Please enter the range [a, b] you want the Error Function worked out to: " << endl;
	cin >> a >> b;

	
	for(double z=a; z <=b; z=z+0.1){
	cout << "z = " << z << " Erf(z) " << SFs::Erf(z) << " ErfFast(z) " << SFs::ErfFast(z) << endl;
	}
	
	cout << " ------ Complement Error function erfc(z) --------" << endl;

	for(double z=a; z <=b; z=z+0.1){
	cout << "z = " << z << " Erfc(z) " << SFs::Erfc(z) << " ErfcFast(z) " << SFs::ErfcFast(z) << endl;
	}
	
}

