/*
Name: TestPolar.cpp
Author: Jay Gambetta

Dependences: UsefullFunctions.hpp
Brief Discription: test the polar to cartesian conversion program in UFs.h
	This being: 
		Cart2Polar(x,y)
		Polar2Cart(amp,theta)
	it takes in a vector in one format an prints it in a new format.
Limitations: None

Version History
	v0: February 12, 2008.
	v1: Janurary 4th, 2009. changed to suport useful functions


*/
#include <UsefullFunctions.hpp>
using namespace std;

int main(void){
	double a, b;
	vector<double> polar(2), cart(2);
	cout << "Please enter x and y: ";
	cin >> a >> b;
	cart[0] = a;
	cart[1] = b;
	cout << "x and y are: " << cart[0] << " " << cart[1] << endl;
	
	polar = USs::Cart2Polar(cart);
	cout << "In polar the magnitude and the angle are: " << polar[0] << " " << polar[1] << endl;
	cout << "Lets convert them back to cartisian \n" << endl;
	cart = USs::Polar2Cart(polar);
	cout << "x and y are: " << cart[0] << " " << cart[1] << endl;
	
	return 0;
}