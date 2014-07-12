/*
Name: UtilityFunctions.hpp
Author: Jay Gambetta

Dependences: the std library  
Brief Discription: a namespace containing my ulitily functions -- error, open file etc.
Limitations: None

Version History
	v0: April 17, 2006
	v1:	February 12, 2008. Change the Polar2Cart and Cart2Polar to support vectors
	v2: November 21, 2008. Starting to remove the std namespace so less confusion can arise
	v3: Janurary 3rd, 2009. Split into three, utility functions, statistical functions and usefull functions
	v4: August 15, 2009.  Added new string to integer/double functions

*/
#ifndef UtilityFunctions_h
#define UtilityFunctions_h
#include <iostream> 
#include <fstream> 
#include <sstream>
#include <string> //std string class
#include <cmath> //cmath declares a set of functions to compute common mathematical operations and transformations
#include <vector> //std vector class
#include <complex> //std complex class
#include <limits> // the std limit libary -- contains highest double, smallest double etc
#include <cstdlib> // would like to remove these programs
#include <cstring>

enum Verbose {yes, no};
enum Verbose verbose;

namespace UFs {//Utility functions `
	
	void MyError(const std::string error_text);
	void FileCheck(std::ifstream&, const std::string);
	void FileCheck(std::ofstream&, const std::string);
	void OpenFile(const std::string filename, std::ofstream& dataout, const std::streamsize prec);
	void AppendFile(const std::string filename, std::ofstream& dataout, const std::streamsize prec);
	std::ifstream& FindInFile(std::ifstream&, const char*);
	std::string dtos(const double, const int);
	std::string itos(const int number);
	double stod(const std::string text);
	int stoi(const std::string text);
	
		
}
//------------Member funtions------------------
inline void UFs::MyError(const std::string error_text){
	//returns an error message that contains the string from which it was called.
	std::cerr << "--------------------------------------------------------------- " << std::endl;
	std::cerr << std::endl;
	std::cerr << "********    You stupid idiot you have made an error    ******** " << std::endl;
	std::cerr << std::endl;
	std::cerr << "--------------------------------------------------------------- " << std::endl;
	std::cerr << "The following error has occurred... " << std::endl;
	std::cerr << std::endl;
	std::cerr << error_text << std::endl;
	std::cerr << std::endl;
	std::cerr << "Exiting ..." << std::endl;
	exit(1);
}
inline void UFs::FileCheck(std::ifstream& entry, const std::string paramfile){
	//Checks if file can be opened
	if (!entry){
		std::cerr << "can not open file: " << paramfile << " for read " << std::endl;
		std::cerr << "Exiting..." << std::endl;
		exit(1);
	}
	if(verbose==yes)
		std::cout << "Data File: " << paramfile << " was loaded correctly " << std::endl;  
}
inline void UFs::FileCheck(std::ofstream& outfile, const std::string paramfile){
	//Checks if file can be opened
	if (!outfile){
		std::cerr << "can not open file: " << paramfile << " for ouputing " << std::endl;
		std::cerr << "Exiting..." << std::endl;
		exit(1);
	}
	if(verbose==yes)
		std::cout << "Data File: " << paramfile << " was loaded correctly " << std::endl;
}
inline void UFs::OpenFile(const std::string filename, std::ofstream& dataout, const std::streamsize prec){
	//opens the file for output
	dataout.open(filename.c_str(), std::ios::out);
	dataout.setf(std::ios::scientific, std::ios::floatfield);
	dataout.precision(prec);
	UFs::FileCheck(dataout,filename.c_str());
}
inline void UFs::AppendFile(const std::string filename, std::ofstream& dataout, const std::streamsize prec){
	//opens the file to append data to it
	dataout.open(filename.c_str(), std::ios::app);
	dataout.setf(std::ios::scientific, std::ios::floatfield);
	dataout.precision(prec);
	UFs::FileCheck(dataout,filename.c_str());
}

//Bellow no examples have been writtent yet....
std::ifstream& UFs::FindInFile(std::ifstream& entry, const char* parameter){
	//Finds a string in a file
	char buff[200]; // declares a string buff of 200 characters (Note this is a limitation of the program 
    entry.seekg(0); // Not sure what this does 
    entry.clear(); // Note sure what this does 
 
    while (!entry.eof()){
		entry >> buff; // assigns entry to buff
	    if(strcmp(parameter, buff) == 0) break; // compares the buff to the search parameter - if 0 exits, if -1 goes on
	}
	
	if(entry.eof()){ // if we cant find it by the end on the file terminate the program
		std::cerr << "Routine UFs::FindInFile: Warning: " << parameter  << " not found." << std::endl;
		std::cerr << "Exiting..." << std::endl;
		exit(1);
	}
	return entry;
}

inline std::string UFs::dtos(const double number, const int number_of_dec){
	std::stringstream strout;
	std::string run;
	
	strout.setf(std::ios::fixed, std::ios::floatfield);
	strout.precision(number_of_dec);
	strout.str("");
	strout << number;
	run = strout.str();
	return run;
}
inline std::string UFs::itos(const int number){
	std::stringstream strout;
	std::string run;
	
	strout.str("");
	strout << number;
	run = strout.str();
	return run;
}
inline double UFs::stod(const std::string text){
	double number;
	std::istringstream stuff(text);
	stuff>>number;
	return number;
}
inline int UFs::stoi(const std::string text){
	int number;
	std::istringstream stuff(text);
	stuff>>number;
	return number;
}


#endif /* UtilityFunctions_h */

