/*
Name: TestSpecialFunctions.cpp
Author: Jay Gambetta

Dependences: SpecialFunctions.hpp
Brief Discription: This program tests the special functions LogGamma, Factorial, GammaP, GammaQ, NChooseR, Beta
Limitations: None

Version History
	v0: February 12, 2008.

*/ 
#include <SpecialFunctions.hpp>
using namespace std;

int main(void){
	cout.setf(ios::scientific, ios::floatfield); 
	cout.precision(30);


	cout << " -------Gamma(z) test ------- "<< endl;
	double z_final;
	cout << "Please enter the final number you want the Gamma Function worked out to: " << endl;
	cin >> z_final;

	
	for(double z=1.00; z <=z_final; z=z+0.1){
	cout << "z = " << z << " Gamma(z) + " << exp(SFs::LogGamma(z)) << endl;
	}
	
	cout << " ------ n! test --------" << endl;

    int n;
	
	cout << "Please enter the final number you want the factorial worked out to: " << endl;
	cin >> n;
	
	for(int i=0; i <= n; i++){
		cout << i << "! = " << SFs::Factorial(i) << endl;
	}

	
	cout << " ------ nCr! test --------" << endl;

    int N;
	
	
	cout << "Please enter the number n in nCr: " << endl;
	cin >> N;
	
	for(int R=0; R <= N; R++){
		cout << N << "C" << R << " = " << SFs::NChooseR(N,R) << endl;
	}

	
	cout << " -------Beta(z,w) test ------- "<< endl;
	double w_final;
	cout << "Please enter the final number you want the Beta Function worked out to: " << endl;
	cin >> w_final;

	
	for(double w=1.00; w <=w_final; w=w+0.1){
	cout << "z = " << z_final << " w = " << w  << " Beta(z,w) + " << SFs::Beta(z_final,w) << endl;
	}
	
	cout << " -------GammaP(z,w) test ------- "<< endl;
	cout << "Please enter the final number you want the GammaP Function worked out to: " << endl;
	cin >> w_final;

	
	for(double w=1.00; w <=w_final; w=w+0.1){
	cout << "z = " << z_final << " w = " << w  << " Beta(z,w) + " << SFs::GammaP(z_final,w) << endl;
	}
	

	cout << " -------GammaQ(z,w) test ------- "<< endl;
	cout << "Please enter the final number you want the GammaQ Function worked out to: " << endl;
	cin >> w_final;

	
	for(double w=1.00; w <=w_final; w=w+0.1){
	cout << "z = " << z_final << " w = " << w  << " Beta(z,w) + " << SFs::GammaQ(z_final,w) << endl;
	}

	//cout << "The Gamma  of -" << z_final << " is " << endl;
	//cout << exp(MyFunctions::LogGamma(-z_final)) << endl;

	cout << "The Factorial of -" << n << " is " << endl;
	cout << SFs::Factorial(-n) << endl;
	return 0;
}

