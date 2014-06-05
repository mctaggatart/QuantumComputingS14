/*
 Name:  Tictoc.hpp
 Author: Botan Khani

 Dependences: none
 Brief Description: Keep track of time using this timer function!  
 					Use tic() to start the clock, toc() to end the clock and 
 					return the time in seconds.
 Limitations: None

 Version History
 v0: June 17, 2010
 */
#ifndef Tictoc_h
#define Tictoc_h

#include <ctime>
//#include <time>
//using namespace std;
namespace std
{
	double CLOCK_START, CLOCK_END, TIMER;
	inline void Tic()
	{
		CLOCK_START=clock();
	}
	inline double Toc()
	{
		CLOCK_END = clock();
		TIMER     = (CLOCK_END-CLOCK_START) / CLOCKS_PER_SEC;

		CLOCK_START=CLOCK_END;
		return TIMER;
	}
};
//------------Member funtions------------------

#endif

