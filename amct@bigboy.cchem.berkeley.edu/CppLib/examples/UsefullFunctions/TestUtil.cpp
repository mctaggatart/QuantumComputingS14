/*
Name: TestUtil.cpp
Author: Jay Gambetta

Dependences: UsefullFunctions.hpp
Brief Discription: This program test the following functions Max, Min, Sign, Swap in UsefullFunctions.hpp
Limitations: None

Version History
	v0: February 12, 2008.
	v1: Janurary 4th, 2009. changed to suport useful functions
	
*/
#include <UsefullFunctions.hpp>
using namespace std;

int main(void){
	

	cout << "Please enter two numbers " << endl;
	double a, b;
	cin >> a >> b;
	
	cout << "The first number is " << a << endl;
	cout << "The second number is " << b << endl;
	cout << "The maximum of them is " << USs::Max<double>(a,b) << endl;
	cout << "The minimum of them is " << USs::Min<double>(a,b) << endl;
	cout << "The sign of the second number is " << USs::Sign<double>(b) << endl;
	
	cout << "We Swap them " << endl;
	USs::Swap<double>(a,b);
	
	cout << "The first number is now " << a << endl;
	cout << "The second number is now " << b << endl;
	cout << "The sign of the second number is " << USs::Sign<double>(b) << endl;
	
	cout << "Lets swap them back" << endl;
	USs::Swap<double>(a,b);

	cout << "The first number is now " << a << endl;
	cout << "The second number is now " << b << endl;
	cout << "The sign of the second number is " << USs::Sign<double>(b) << endl;
	
	return 0;
}

