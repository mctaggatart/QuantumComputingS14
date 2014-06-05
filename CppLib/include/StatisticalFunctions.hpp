/*
Name: StatisticalFunctions.hpp
Author: Jay Gambetta

Dependences: UtilityFunctions.hpp  
Brief Discription: a namespace containing some functions used to do statistics on some data (vector)
Limitations: None

Version History
	v3: january 3rd, 2009
	v4: september 3rd, 2009
		Added some new functions, added capability to use arrays as alternative
		
*/
#ifndef StatisticalFunctions_h
#define StatisticalFunctions_h
#include "UtilityFunctions.hpp"

namespace STs {//Utility functions 
	
	double Mean(const std::vector<double>& A);
	double Mean(const double* A, const double n);
	double STD(const std::vector<double>& A);
	double STD(const double* A, const double n);
	double Sum(const std::vector<double>& A);
	double Sum(const double* A, const double n);
	double MaxElement(const std::vector<double>& A);
	std::vector<double> Histogram(const std::vector<double>& A, const double bin_min, const double bin_max, const size_t bin_no);
	
}
//------------Member funtions------------------
inline double STs::Mean(const std::vector<double>& A){
	double temp = 0.0;
	for (size_t i=0; i < A.size(); i++){
		temp += A[i];
	}
	return temp/double(A.size());
}
inline double STs::Mean(const double* A, const double n){
	double temp = 0.0;
	for (size_t i=0; i < n; i++){
		temp += A[i];
	}
	return temp/n;
}
inline double STs::STD(const std::vector<double>& A){
	double temp = 0.0;
	double mean = Mean(A);
	for (size_t i=0; i < A.size(); i++){
		temp += pow((A[i]-mean),2);
	}
	return temp/double(A.size());
}
inline double STs::STD(const double* A, const double n){
	double temp = 0.0;
	double mean = Mean(A, n);
	for (size_t i=0; i < n; i++){
		temp += pow((A[i]-mean),2);
	}
	return temp/n;
}
inline double STs::Sum(const std::vector<double>& A){
	double temp = 0.0;
	for (size_t i=0; i < A.size(); i++){
		temp += A[i];
	}
	return temp;
}
inline double STs::Sum(const double* A, const double n){
	double temp = 0.0;
	for (size_t i=0; i < n; i++){
		temp += A[i];
	}
	return temp;
}
inline double STs::MaxElement(const std::vector<double>& A){
	double temp = A[0];
	for (unsigned i=1; i < A.size(); i++){
		if( A[i]> temp )
			temp=A[i];
	}
	return temp;
}
inline std::vector<double> STs::Histogram(const std::vector<double>& A, const double bin_min, const double bin_max, const size_t number_of_bins){
	std::vector<double> hist(number_of_bins+1);//Need to understand the +1
	double bin_step = (bin_max-bin_min)/double(number_of_bins);
	std::cout << "The bin size is: " << bin_step << std::endl;
	size_t bin;
	for (size_t i=0; i < A.size(); i++){
		bin = int(floor((A[i] -bin_min)/bin_step));
		hist[bin] +=1;
	}
	return hist;
}

#endif /* StatisticalFunctions_h */

