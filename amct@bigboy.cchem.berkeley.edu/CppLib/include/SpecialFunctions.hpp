/*
Name: SpecialFunctions.hpp
Author: Jay Gambetta

Dependences: UFs.h
Brief Discription: a namespace containing all the special functions i have written
Limitations: some of these special functions exists in the std library so i would use them 
	over these as they may contain higher order approximations then contian in the 
	simple examples I use here .

Version History
	v0: April 17, 2006
	v1:	February 12, 2008. Rewrote into a namespace from indivual c code. 
	v2: Janurary 4th 2009. Change the name to SpecialFunctions.hpp
*/
#ifndef SpecialFunctions_h
#define SpecialFunctions_h
#include "UtilityFunctions.hpp"

namespace SFs {// Special Functions
	double LogGamma(const double z);
	double Factorial(const int n);
	double NChooseR(const int n, const int r);
	double LogFactorial(const int n);
	double Beta(const double z, const double w);
	double GammaP(const double a, const double x);
	double GammaQ(const double a, const double x);
	int GammaPSeries(double& gamma_p_ser, const double a, const double x, double& log_gamma);
	int GammaQConFrac(double& gamma_q_cf, const double a, const double x, double& log_gamma);
	double Erf(const double x);
	double Erfc(const double x);
	double ErfFast(const double x);
	double ErfcFast(const double x);
}
//------------Member funtions------------------
inline double SFs::LogGamma(const double z){
	/*	Info:
	 		This function returns log|Gamma(z)| for Re z > 0
	 		We use the Lanczos approximation 
	 			Gamma(z+1) = (z+gamma+0.5)^{z+0.5}exp(-(z+gamma+0.5))*sqrt{2 pi} [c_0+\sum_i=1 c_i/(z+i)+eps]
	 			-> log(Gamma) 	= log[Gamma(z+1)/z] 
	  							= (z+0.5)*log(z+gamma+0.5) -(z+gamma+0.5) + log(sqrt{2 pi}[c_0+\sum_i=1 c_i/(z+1)])
	 		At the sixth order the parameters are
	 			gamma = 5, N = 6, 
	 			c = {1.000000000190015, 76.18009172947146, -86.50532032941677, 24.01409824083091, 
	 				-1.231739572450155, 1.2108650973866179e-3, -5.395239384953e-6 }
	 		At this order the error |eps| < 2 x 10^-10 
			Also works with complex numbers and has this bound 
	 	Limitations: 
	 		1. It is a approximation with error 2 x 10^-10
	*/
	double z_temp, temp, sum;
	static const double coef[7]={1.000000000190015, 76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155, 1.2108650973866179e-3, -5.395239384953e-6 };
	
	if (z<0) UFs::MyError("Routine SFs::LogGamma: Currently this implementation of the Gamma Function only supports positive values");
	z_temp = z;
	temp = z+5.5;  // this is the first part of the approximation:
	temp -= (z+0.5)*log(temp); 
	sum = coef[0]; // this is the sum part
	for(int j=1; j<7;j++) sum += coef[j]/++z_temp;
	return -temp+log(2.5066282746310005*sum/z);
}
inline double SFs::Factorial(const int n){
	/*	Dependences: SFs::LogGamma
	 	Info:
	 		This function returns n! as a double, It uses the static storage so that once it has 
	 		been called once the computation is not repeated.
	 	Limitations:
	 		1. It has roundoff error since we use double. 
	 		2. It is recursive up to 32! so the roundoff error will multiple
	 		3. After 32 the accuracy is limited by the Lanczros approximation for a Gamma function 
	 			-see SFs::LogGamma
	*/
	static int ntop=4; // declares ntop as an integer an set it for the first coll to be 4, for the second call it will have the value the function exit with last time.
	static double a[33] = {1.0, 1.0, 2.0, 6.0, 24.0}; // the other 28 are zero

	int j;
	if (n<0) 
		UFs::MyError("Routine SFs::Factorial: The factorial of a negative value is complex infinity"); 
	else if (n>32) 
		return exp(SFs::LogGamma(double(n+1.0))); // N! = Gamma(n+1) and for large values we will assume that the error in Gamma(n+1) in the Lanczos approximation is less than recursive numerical error
	while (ntop < n){

		j = ntop++;  //for the first call j = 4 ntop = 5 so we assign ntop to j then add one to ntop 
		//cout << j << " " << ntop++ << endl;

		a[ntop] = a[j]*ntop;
	}

	return a[n];
}
inline double SFs::NChooseR(const int n, const int r){
	/*	Dependences: SFs::LogFactorial	
		Info:
			This function returns the binomial coefficient nCr 
	 			nCr = \frac{n!}{r! (n-r)!}   where n >= r >= 0
	 	
			We don't use Factorial.cpp as this overflow (at n = 170 on my mac) whereas 
			the binomial coefficient still has a machine size number at these large n
		
			We use instead nCr = exp(log(n!) - log(r!) - log((n-r)!)) and rather then just 
			call LagGamma(n+1) every time for log(n!) we make up a look up table using 
			SFs::LogFactorial	
		Limitations:
	 	1. It is limited by the 6th order Lanczros approximation for a Gamma Function - see SFs::LogGamma 
	*/
	
	// The floor cleans up roundoff errors at small values of n and r
	return floor(0.5 + exp(SFs::LogFactorial(n)-SFs::LogFactorial(r)-SFs::LogFactorial(n-r))); 
}
inline double SFs::LogFactorial(const int n){
	/*	Dependences: SFs::LogGamma
		Info: 
	 		This routine creates a look up table of the log(n!) using the Lanczos approximation 
	 		for the gamma function - see LogGamma.cpp 
	 		The Gamma function is related to the factorial by n! = Gamma(n+1)
	 	Limitations:	
	 		1. Uses the Lanczos approximation for Gamma(n+1) 				
	*/
	static double a[101]; //a static array that remember previous calls to speed up tho calculation - all values are set to zero to start

	if (n<0) 
		UFs::MyError("Routine SFs::LogFactorial: The factorial of a negative value is complex infinity");
	else if(n<=1) 
		return 0.0;
	else if (n <= 100) // in the range of the table 
		return (a[n]!=0.0 ? a[n] :(a[n]=SFs::LogGamma(n+1.0))); // if it has been cal before return the result if not do the calculation
	else 
		return SFs::LogGamma(n+1.0); //
}
inline double SFs::Beta(const double z, const double w){
	/*	Dependences: LogGamma.cpp
	 	Info:
	 		This function returns the beta function using the Lanczos approximation for 
			the Gamma function
	 			Beta(z,w) = \frac{Gamma(z)*Gamma(w)}{Gamma(z+w)}
	 			-> Beta = exp(log(Gamma(z))+ log(Gamma(w))-log(Gamma(z+w)))
	 	Limitations:
	 		1. It is limited to the accuracy of the Lanczos approximation - see SFs::LogGamma
	*/
	return exp(SFs::LogGamma(z) + SFs::LogGamma(w) - SFs::LogGamma(z+w));
}
inline double SFs::GammaP(const double a, const double x){
	/*	Dependences: GammaPSeries.cpp GammaQConFrac.cpp 
		Info:
			This routine returns the incomplete Gamma function P(a,x)   
			P(a,x) 	= \gamma(a,x)/\Gamma(a)  
					= \int_0^x e^{-t} t^{a-1} dt     (a>0)
			using either the series or continued fraction method depending on 
	 			x < a + 1 -> series 
	 			x > a + 1 -> continued fraction method 
		Limitations:
			1. Uses iterative methods - which does at most 100 iterations
			2. Uses the Lanzcros approximation for the Gamma function - see SFs::LogGamma
	*/ 
	
	double gamma_p=0.0, gamma_q, log_gamma;
	
	if (x<0.0 || a < 0.0)
		UFs::MyError("Routine SFs::GammaP: input parameters must be above zero");
	if (x < a+ 1.0){
		GammaPSeries(gamma_p, a, x, log_gamma);
		return gamma_p;
	}
	else {
		GammaQConFrac(gamma_q, a, x, log_gamma);
		return 1.0 - gamma_q;
	}
}
inline double SFs::GammaQ(const double a, const double x){
	/*	Dependences: SFs::LogGamma SFs::GammaPSeries SFs::GammaQConFrac
		Info:
	 		This routine returns the incomplete Gamma function Q(a,x)   
				Q(a,x) 	= \Gamma(a,x)/\Gamma(a)  
	 					= \int_x^infty e^{-t} t^{a-1} dt     (a>0)
	 		using either the series or continued fraction method depending on 
	 			x < a + 1 -> series 
	 			x > a + 1 -> continued fraction method 
		Limitations:
	 		1. Uses iterative methods - which does at most 100 iterations
	 		2. Uses the Lanzcros approximation for the Gamma function - see SFs::LogGamma
	*/
	double gamma_p=0.0, gamma_q, log_gamma;
	
	if (x<0.0 || a < 0.0)
		UFs::MyError("Routine SFs::GammaP: input parameters must be above zero");
	if (x < a+ 1.0){
		GammaPSeries(gamma_p, a, x, log_gamma);
		return 1.0-gamma_p;
	}
	else {
		GammaQConFrac(gamma_q, a, x, log_gamma);
		return gamma_q;
	}
}
inline int SFs::GammaPSeries(double& gamma_p_ser, const double a, const double x, double& log_gamma){
	/*	Dependences: LogGamma.cpp 
	 	Info:
	 		This routine returns the incomplete Gamma function P(a,x) evaluated by its 
			series representation 
	 			P(a,x) = \gamma(a,x)/\Gamma(a)  
	 		where 
	 			\gamma(a,x) = exp(-x)x^a\sum_0^infty \frac{\Gamma(a)}{\Gamma(a+1+n)} x^n 
	 			 = exp(-x +a*log(x)-log(Gamma(a)))*\sum\frac{x^n}{\Prod_{m=0}^n (a+n-m) }
	 	Limitations:
	 		1. Uses a series approximation - which does at most 100 iterations
	 		2. Uses the Lanzcros approximation for the Gamma function - see LogGamma.cpp
	*/
	const int ITMAX= 100; // maximum number of iterations
	const double EPS = std::numeric_limits<double>::epsilon(); //the smallest possible double that con be added to 1 on this machine
	double sum, del, ap; 

	log_gamma = SFs::LogGamma(a);
	if (x<= 0.0){ 
		if (x < 0.0) 
			UFs::MyError("Routine SFs::GammaPSeries: x is less than 0"); // negative values cant be handled by this routine 
		gamma_p_ser = 0.0;
		return 0;
	}
	else {
		ap = a;
		del=sum=1.0/a;  
		for (int n=0; n < ITMAX; n++){
			++ap; //adds one first
			del *= x/ap; // this get the current x^n/(a*(a+1) ..(a+n)
			sum += del; 
			if (fabs(del) < fabs(sum)*EPS){ // not sure why this is guaranteed to converge
				gamma_p_ser = sum*exp(-x+a*log(x)-log_gamma);
				return 0; 
			}
		}
		UFs::MyError("Routine SFs::GammaPSeries: a is too large, ITMAX to small");
		return 0;
	}
}
inline int SFs::GammaQConFrac(double& gamma_q_cf, const double a, const double x, double& log_gamma){
	/*	Dependences: SFs::LogGamma
		Info:
			This routine returns the incomplete Gamma function Q(a,x) = 1 - P(a,x)  evaluated by its 
			Continued fraction representation  
				Q(a,x) = \Gamma(a,x)/\Gamma(a)  
		Limitations:
			1. Uses a continued fraction method - which does at most 100 iterations
			2. Uses the Lanzcros approximation for the Gamma function - see SFs::LogGamma.cpp
	*/
	const int ITMAX=100; // the maximum allowed iterations
	const double EPS=std::numeric_limits<double>::epsilon(); //is the relative accuracy
	const double DPMIN=std::numeric_limits<double>::min()/EPS; //is the smallest representable double
	int i;
	double an, b, c, d, del, h;

	log_gamma = LogGamma(a);
	b=x+1.0-a;
	c=1.0/DPMIN;
	d=1.0/b;
	h=d;
	for (i=1; i<=ITMAX; i++){
		an=-i*(i-a);
		b += 2.0;
		d = an*d+b;
		if (fabs(d) < DPMIN) d = DPMIN;
		c = b+an/c;
		if (fabs(c) < DPMIN) c = DPMIN;
		d = 1.0/d;
		del=d*c;
		h *= del;
		if(fabs(del-1.0) <= EPS) break;
	}
	if (i > ITMAX) UFs::MyError("Routine SFs::GammaQConFrac: a is too large, ITMAX is to small");
	gamma_q_cf=exp(-x+a*log(x)-log_gamma)*h;
	return 0; 
}
inline double SFs::Erf(const double x){
	/*	Dependences: GammaP.cpp  
		Info:
	 		This routine returns the Error function erf(x)   
	 			erf(x) = \sqrt{\frac{2}{pi}} \int_0^x \exp(-t^2) dt
			We dont calculate this by solving the integral: Instead we use the incomplete Gamma function P(x,a) 
				- See SFs::GammaP 
			erf(x) =P(0.5,x^2)  (x>=0)
		Limitations:
			1. Uses iterative methods - which does at most 100 iterations
			2. Uses the Lanzcros approximation for the Gamma function - see SFs::LogGamma
	*/
	return x < 0.0 ? -SFs::GammaP(0.5,x*x) : SFs::GammaP(0.5,x*x);
}
inline double SFs::Erfc(const double x){
	/*	Dependences: GammaP.cpp  
	 	Info:
	 		This routine returns the complement of the Error function erfc(x)   
				erfc(x) = \sqrt{\frac{2}{pi}} \int_x^\infty \exp(-t^2) dt
	 		We don't calculate this by solving the integral: Instead we use the incomplete Gamma function
	 			P(x,a) - See SFs::GammaP 
	 			erfc(x) =1-P(0.5,x^2)  (x>=0)
	 	Limitations:
	 		1. Uses iterative methods - which does at most 100 iterations
	 		2. Uses the Lanzcros approximation for the Gamma function - see LogGamma.cpp
	*/
	return x < 0.0 ? 1.0+SFs::GammaP(0.5,x*x) : 1.0-SFs::GammaP(0.5,x*x);
}
inline double SFs::ErfFast(const double x){
	/*	This routine returns the complement Error function erfc(x) with 
		fractional error everywhere less then 1.2x10^-7
	 */
	double t, z, ans;
	z =fabs(x);
	t = 1.0/(1.0 +0.5*z);
	ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*(0.17087277))))))))));
	return (x >= 0.0 ? 1.0-ans : ans-1.0);
}
inline double SFs::ErfcFast(const double x){
	/*	This routine returns the complement Error function erfc(x) with 
	 	fractional error everywhere less then 1.2x10^-7
	*/ 
	double t, z, ans;
	z =fabs(x);
	t = 1.0/(1.0 +0.5*z);
	ans = t*exp(-z*z-1.26551223+t*(1.00002368+t*(0.37409196+t*(0.09678418+t*(-0.18628806+t*(0.27886807+t*(-1.13520398+t*(1.48851587+t*(-0.82215223+t*(0.17087277))))))))));
	return (x >= 0.0 ? ans : 2.0-ans);
}
#endif /* SpecialFunctions_h */

