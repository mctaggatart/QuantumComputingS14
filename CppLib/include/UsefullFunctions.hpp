/*
Name:  UsefullFunctionss.hpp
Author: Jay Gambetta/Botan Khani

Dependences: UtilityFunctions.hpp  
Brief Discription: a namespace containing same usefull functions. They may be simply but they 
	do come in handy
Limitations: None

Version History
	v0: Januarary 3rd, 2009
	v1: September 12, 2009, added OutputMatrix and the load matrix file (cannot load imaginary numbers at the moment)
	v2: February 25, 2010: 	- Fixed load matrix file so it CAN load complex and imaginary numbers.  
								It's really to open and save matrices with this tool, just read the comments in the function. 
								Highly recommended.
							- Added sleep to make the CPU "sleep", useful for waiting for slow things like gnuplot 
							- Added LoadTable, which loads columns from a delimited file into a double pointer
							- Added BubbleSort, which sorts data1 in ascending order, can also force data2 to be sorted according to data1
							- Added FindMax
							- Added LoadVector, loads data from a delimited file into a Vector data type
							- Added LoadVectorCount, gives you the length of the data in a delimited file
							
*/
#ifndef UsefullFunctions_h
#define UsefullFunctions_h
#include "UtilityFunctions.hpp"
#include "Matrix.hpp"
#include <fstream>

namespace USs {//Utility functions 
	
	template<class T> 
	inline T Max(const T a, const T b){
		// returns the maximum of a and b for any type
		return b>a ? b : a;
	}
	template<class T>
	inline T Min(const T a, const T b){
		// returns the minimum of a and b for any type
		return b>a ? a : b;
	}
	template<class T> 
	inline void Swap(T &a, T &b){
		// sways the values of a and b for any type
		T temp =a; a=b; b=temp;
	}
	template<class T>
	inline int Sign(T&a){
		//returns the sign of a as an integer
		return (a>=0 ? 1 : -1);
	}
	std::vector<double> Polar2Cart(const std::vector<double>& polar);
	std::vector<double> Cart2Polar(const std::vector<double>& cart);
	template <class T>
	inline void OutputMatrix(matrix<T>& rho);
	inline double FindMax(double* y_values, int N);
	inline void FourierBubbleSort(double* &data1, int n);
	inline void FourierBubbleSort(double* &data1, double* &data2, int n);
	inline void Sleep(int seconds);
	inline vector<double> 	LoadVector		(const std::string filename, const std::string delimiter);
	inline int 				LoadVectorCount	(const std::string filename, const std::string delimiter);
	
		
	
	//File manipulation
	inline double** 						LoadTable		(const std::string filename, const std::string delimiter, const size_t max_col_size, int& N);
	inline matrix<std::complex<double> >** 	LoadMatrixFile	(const std::string filename, size_t& num_groups, size_t& group_size, size_t& rows, size_t& cols);
	
	void MakeGnuplotScript(size_t nlines, char* filename, char* xlabel, char* ylabel, char* datafile, bool islog);
	
	// pulse shapes -- depracated -- use Contol class instead
	inline double SquarePulse(const double t, const double ton, const double td, const double amp){
		// height is amp and duration is td .. thus area is amp*td
		//n pi/2 pulse area = n pi/2 -> amp  = n pi/2td
	 	 return (t>ton && t < td+ton) ? amp : 0.0;
	}
	inline double TanhPulse(const double t, const double ton, const double td,  const double amp, const double sigma){
		// height is amp and duration is td area is approximately amp*td
		//n pi/2 pulse area = n pi/2 -> amp  = n pi/2td
		double one_on_sigma=1.0/sigma;
		return 0.5*amp*(std::tanh(one_on_sigma*(t-ton)) + std::tanh(-one_on_sigma*(t-ton-td)) );
	}
	inline double Gaussian(const double t, const double ton, const double td, const double amp, const double sigma){
		//normalized so that the area is amp
		// npi/2 pulse is done by timeing by n pi/2
		return amp*std::exp(-0.5*pow((t-ton-td*0.5),2)/std::pow(sigma,2)) / (sigma*M_SQRT2*std::sqrt(M_PI));
	}	
	inline double ShiftedGaussian(const double t, const double ton, const double td, const double amp, const double sigma){
		//normalized so that the area is amp
		// npi/2 pulse is done by timeing by n pi/2
		//starts and ends at 0
		return M_PI*amp*(std::exp(-0.5*pow((t-ton-td*0.5)/sigma,2)) - std::exp(-0.5*pow((-ton-td*0.5)/sigma,2)) ) / (sigma*M_SQRT2*std::sqrt(M_PI))/erf(td*0.5*M_SQRT1_2/sigma);
	}
	
	inline double GaussianDerivative(const double t, const double ton, const double td, const double amp, const double sigma){
		//normalized so that the area is amp
		// npi/2 pulse is done by timeing by n pi/2
		return -2*0.5*(t-ton-td*0.5)/std::pow(sigma,2) * M_PI*amp*std::exp(-0.5*pow((t-ton-td*0.5)/sigma,2)) / (sigma*M_SQRT2*std::sqrt(M_PI))/erf(td*0.5*M_SQRT1_2/sigma);
	}
	
	inline double TruncatedGaussian(const double t, const double ton, const double td, const double amp, const double sigma){
		//normalized to approximately amp
		// npi/2 pulse is done by timeing by amp = n pi/2		if(t>ton && t<td+ton)
			return amp*std::exp(-0.5*std::pow((t-ton-td*0.5),2)/std::pow(sigma,2))/(sigma*M_SQRT2*std::sqrt(M_PI))/erf(td*0.5*M_SQRT1_2/sigma);	}
	
	//PhaseRamping
	
	inline double RandControl(const double min, const double max){
		srand(time(0));
		double ctrl = rand();
		return min + (ctrl/RAND_MAX)*(max - min);
	}
	
	//RandControlRamped
	
	
	
//	inline double StarkShift(amp, delta)
//	{
	
//	}
	
	
	// double TanhTanPulse(const double t, const double ton, const double t, const double amp, const double sigma){  
	// 	// the raw tanh-tan pulse function
	//     if (ton<=t && t<ton+10.0)
	// 		return amp*(std::tanh(fp_steepness*std::tan((t-ton-5.0)*M_PI/fp_edgetime))+1.0)/2.0;
	// 	else if (ton+10.0<=t && t<toff)
	// 		return amp;
	// 	if (t>td && t<=td)
	// 		return amp*(std::tanh(fp_steepness*std::tan((-t+ton+5.0)*M_PI/fp_edgetime))+1.0)/2.0;
	// 	else return 0.0;
	// }
	
}
//------------Member funtions------------------
inline std::vector<double> USs::Polar2Cart(const std::vector<double>& polar){
	//Converts a polar vector (amplitude,angle) to cartension vector (x,y)
	std::vector<double> cart(2); 
	cart[0] = polar[0]*std::cos(polar[1]);
	cart[1] = polar[0]*std::sin(polar[1]);
	return cart;
}
inline std::vector<double> USs::Cart2Polar(const std::vector<double>& cart){
	//Converns a cartension vector (x,y) to a polar vector (amplitude, angle)
	std::vector<double> polar(2);
	
	//Angle, first we have to get the correct quadrature (all stations to central)
	if (cart[0] > 0 && cart[1] >= 0) //quad 1
		polar[1] = std::atan(cart[1]/cart[0]);
	else if (cart[0] < 0 && cart[1] >= 0 ) // quad 2
		polar[1] = M_PI - std::atan(-cart[1]/cart[0]);
	else if (cart[0] < 0 && cart[1] <= 0) //quad 3
		polar[1] = M_PI + std::atan(cart[1]/cart[0]);
	else if (cart[0] > 0 && cart[1] <= 0) //quad 4
		polar[1] = 2*M_PI - std::atan(-cart[1]/cart[0]);
	else if (cart[0] == 0 && cart[1] > 0) // along x axis and above
		polar[1] = M_PI*0.5;
	else if (cart[0] == 0 && cart[1] < 0) // along x axis and below
		polar[1] = - M_PI*0.5;
	else 
		polar[1] = 0.0;// define the case when both are zero to be 0
	
	//amplitude 
	polar[0] = std::sqrt(std::pow(cart[0],2)+std::pow(cart[1],2));
	return polar;
}
//Output matrix with complex numbers (Matrix::SetOutputStyle can not do this)
template <class T>
inline void USs::OutputMatrix(matrix<T>& rho){
	//Initial fock state, rho=|n><n| 
	int num_row = rho.GetRows();
	int num_col = rho.GetColumns();
	for (int row=0; row < num_row; row++)
	{
		for (int col=0; col < num_col; col++)
		{
			std::cout << rho(row,col) << "\t";
		}
		std::cout << "\n";
	}
}
//Find the max of a pointer array of size N.
inline double FindMax(double* y_values, int N)
{
	double y_max	=y_values[0];
	//First find the maximum y_value
	for (int i=1; i < N; ++i)
		if ( ABS(y_values,i) > y_max) y_max=ABS(y_values,i);
		
	return y_max;
}
//Hey, sort data1 in ascending order!
inline void FourierBubbleSort(double* &data1, int n)
{
	int SubArrayEnd = n -1;
	double temp;

	while (SubArrayEnd > 0)
	{
		int nextEnd = 0;
		
		for (int j = 0; j < (n - 1); ++j)
		{
			if (data1[j] > data1[j+1])
			{
				temp		 = data1[j];
				data1[j] 	 = data1[j+1]; 
				data1[j+1]	 = temp;
				
				nextEnd = j;
			}
		}
	
		SubArrayEnd = nextEnd;
	}
}
//Hey, it sorts your two sets of data according to data1 in ascending order!
inline void FourierBubbleSort(double* &data1, double* &data2, int n)
{
	int SubArrayEnd = n -1;
	double temp;

	while (SubArrayEnd > 0)
	{
		int nextEnd = 0;
		
		for (int j = 0; j < (n - 1); ++j)
		{
			if (data1[j] > data1[j+1])
			{
				temp		 = data1[j];
				data1[j] 	 = data1[j+1]; 
				data1[j+1]	 = temp;
				
				temp		 = data2[2*j];
				data2[2*j] 	 = data2[2*j+2];
				data2[2*j+2] = temp;
				
				nextEnd = j;
			}
		}
	
		SubArrayEnd = nextEnd;
	}
}
//Make the CPU busy for a few seconds
inline void Sleep(int seconds)
{
	time_t start_time, cur_time;

    time(&start_time);
    do
    {
    	time(&cur_time);
	}
    while((cur_time - start_time) < seconds);
}
//Loads a vector from a simple delimited file (typically comma)
inline vector<double> LoadVector(const std::string filename, const std::string delimiter)
{
	//This is pretty straight forward to use, each element in the file is 
	//delimited by a delimiter.  You can choose any delimiter, it doesn't even
	//have to be a single character.  Currently only reads real numbers.
	std::string line, temp;
	std::ifstream myfile (filename.c_str());
	
	size_t commaLoc, currentLoc=0, prevLoc=0;
	size_t d_size=delimiter.size();
	size_t index=0;
	size_t max_size=1000, inc_size=1000;
	
	double* x;
	x=new double[max_size];
	//bool verbose=yes;
	while (! myfile.eof() )		//Loop through each row
	{
		getline (myfile,line);
		size_t string_size=line.size();
		
		while(currentLoc < string_size)		//Loop through each delimited element in the line
		{
			currentLoc=line.find(delimiter, prevLoc); 
			
			std::istringstream stuff(line.substr(prevLoc, currentLoc-prevLoc-d_size+1));	//Real number
			stuff>>x[index];
			
			index   = index + 1;

			//Reallocate the array if it's greater than the maximum size
			if (index > max_size)
			{
				double* temp;
				temp=new double [max_size];	
				for (size_t j = 0; j < max_size; j++)
				  temp[j] = x[j];
				delete [] x; 
				
				x=new double[max_size+inc_size];
				for (size_t j = 0; j < max_size; j++)
					x[j] = temp[j];
				max_size+=inc_size;
				delete [] temp;
			}

			prevLoc = currentLoc + d_size;
		}
	}
	
	index--;
	
	//Create the vector, now that we know the size
	vector<double> Vector(index);

	//Transfer the array to the vector
	for (size_t j=0; j<index; j++)
	{
		Vector[j]=x[j];
		if (verbose==yes) std::cout << "X" << j << ": " << Vector[j] << "\n";
	}
	
	return Vector;
	
	//Garbage collection
	delete [] x;
}
//Outputs the length of a vector from a simple delimited file (typically comma)
inline int LoadVectorCount(const std::string filename, const std::string delimiter)
{
	//This is pretty straight forward to use, each element in the file is 
	//delimited by a delimiter.  You can choose any delimiter, it doesn't even
	//have to be a single character.  Currently only reads real numbers.
	std::string line, temp;
	std::ifstream myfile (filename.c_str());
	
	size_t commaLoc, currentLoc=0, prevLoc=0;
	size_t d_size=delimiter.size();
	size_t index=0;
	size_t max_size=1000, inc_size=1000;
	
	double* x;
	x=new double[max_size];
	//bool verbose=yes;
	while (! myfile.eof() )		//Loop through each row
	{
		getline (myfile,line);
		size_t string_size=line.size();
		
		while(currentLoc < string_size)		//Loop through each delimited element in the line
		{
			currentLoc=line.find(delimiter, prevLoc); 
			
			std::istringstream stuff(line.substr(prevLoc, currentLoc-prevLoc-d_size+1));	//Real number
			stuff>>x[index];
			
			index   = index + 1;

			//Reallocate the array if it's greater than the maximum size
			if (index > max_size)
			{
				double* temp;
				temp=new double [max_size];	
				for (size_t j = 0; j < max_size; j++)
				  temp[j] = x[j];
				delete [] x; 
				
				x=new double[max_size+inc_size];
				for (size_t j = 0; j < max_size; j++)
					x[j] = temp[j];
				max_size+=inc_size;
				delete [] temp;
			}

			prevLoc = currentLoc + d_size;
		}
	}
	
	return index-1;
	
	//Garbage collection
	delete [] x;
}
//Loads a multiple columns containing vectors from a simple delimited file (typically comma)
inline double** LoadTable(const std::string filename, const std::string delimiter, const size_t max_col_size, int& N)
{
	//This is pretty straight forward to use, each element in the file is 
	//delimited by a delimiter.  You can choose any delimiter, it doesn't even
	//have to be a single character.  Currently only reads real numbers.
	std::string line, temp;
	std::ifstream myfile (filename.c_str());
	
	size_t commaLoc, currentLoc=0, prevLoc=0;
	size_t d_size=delimiter.size();
	size_t row_index=0, col_index=0;
	size_t max_size=1000, inc_size=1000;

	double** x;
	x=new double* [max_size];
	for (size_t j=0; j<max_size; ++j)
		x[j]= new double [max_col_size];
	bool verbose=0;
	
	//Check if the file exists
	if (myfile.is_open ()==false) 
	{	
		std::cout << "\nOops, file " << filename << " does not exists.\n";
		N=1;
		return x;
		exit(1);
	}
		
	while (! myfile.eof() )		//Loop through each row
	{
		getline (myfile,line); 
		if (verbose==1) cout << row_index << ":" << line << "\n";
		size_t string_size=line.size();
		col_index=0;
		prevLoc=0;
		currentLoc=0;

		while(currentLoc < string_size)		//Loop through each delimited element in the line
		{
			if (col_index >= max_col_size) break;

			currentLoc=line.find(delimiter, prevLoc); 
			
			std::istringstream stuff(line.substr(prevLoc, currentLoc-prevLoc-d_size+1));	//Real number
			stuff>>x[row_index][col_index]; 
			if (verbose==1) cout << x[row_index][col_index] << " || ";
			
			col_index   = col_index + 1;

			prevLoc = currentLoc + d_size;
		}
		if (verbose==1) cout << "\n";
		
		row_index=row_index+1;
		//Reallocate the array if it's greater than the maximum size
		if (row_index > max_size-1)
		{
			double** temp;
			temp=new double* [max_size];	
			for (size_t j = 0; j < max_size; j++)
			{
				temp[j] = new double [max_col_size];
				for (size_t k = 0; k < max_col_size; k++)
				{
					temp[j][k] = x[j][k];
				}
			}
			delete [] x; 

			x=new double*[max_size+inc_size];
			for (size_t j = 0; j < max_size; j++)
			{
				x[j] = new double [max_col_size];
				for (size_t k = 0; k < max_col_size; k++)
				{
					x[j][k] = temp[j][k];
				}
			}
			for (size_t j = max_size; j < max_size+inc_size; j++)
				x[j] = new double [max_col_size];
			max_size+=inc_size;
			delete [] temp;
		}
	}

	if (col_index >= max_col_size) std::cout << "You idiot!  There are more columns in this file than you specified\n";
	N=row_index-1;
	return x;
	
	//Garbage collection
	delete [] x;
}
//Load matrix data from a formatted file.  The formatted file can easily be created in any program such as Matlab.
inline matrix<std::complex<double> >** USs::LoadMatrixFile(const std::string filename, size_t& num_groups, size_t& group_size, size_t& rows, size_t& cols)
{
	/* Matrice file loader
		This function works by opening a file containing groups of matrices
		and placing their elements into an array of matrix classes
		
		The format of the file is specified as follows:
		-Use any delimiter as long as it is defined in the function parameter
		-The first line of the file must contain the following positive integer
		 parameters separated by spaces:
			[Number of Groups] [Number of Matrices per Group] [row dimension] [col dimension]
			ie. If I had a two control Hamiltonians for ten points to optimize, 
			    and my system had 20 levels, then the first line of the file would be:
				10 2 20 20
		-Subsequent lines contain matrix elements, one row is one line, columns
		 are delimited.  Each element that has an imaginary number must 
		 contain the real number, a comma, and the imaginary number with no spaces
		 in between.  Otherwise they may contain a real number without any commas
		 for the imaginary.  There are no spaces or lines between different matrices or groups,
		 it is assumed that they have the same dimensions.  Any empty line between the matrices
		 is considered the start of a new group.
	*/
	
	std::string line, temp;
	size_t commaLoc, currentLoc, prevLoc=0;
	matrix<std::complex<double> >** Matrix;
	size_t row=0, col, index1=0, index2=0;
	
	std::ifstream myfile (filename.c_str());	
	getline (myfile,line);

	//Get the number of groups
	currentLoc=line.find(" ", prevLoc);
	temp=line.substr(prevLoc, currentLoc);
	num_groups=UFs::stoi(temp);
	prevLoc=currentLoc;

	//Get the group size	
	currentLoc=line.find(" ", prevLoc+1);
	temp=line.substr(prevLoc, currentLoc);
	group_size=UFs::stoi(temp);
	prevLoc=currentLoc;
	
	//Get the number of rows
	currentLoc=line.find(" ", prevLoc+1);
	temp=line.substr(prevLoc, currentLoc);
	rows=UFs::stoi(temp);
	prevLoc=currentLoc;
	
	//Get the number of columns
	currentLoc=line.find(" ", prevLoc+1);
	temp=line.substr(prevLoc, currentLoc);
	cols=UFs::stoi(temp);
	prevLoc=currentLoc;
		
	//Initialize the matrix array
	Matrix=new  matrix<std::complex<double> > *[num_groups];
	
	for (size_t group=0; group < num_groups; ++group){
		Matrix[group] = new  matrix<std::complex<double> > [group_size];
				
		for(size_t k=0; k < group_size; ++k){
			Matrix[group][k].initialize(rows,cols);
		}
	}
	
	//Obtain the elements for the matrix array from the file
	while (! myfile.eof() )		//Loop through each row
	{
		getline (myfile,line);
		
		if (row >= rows)
		{
			index2++;
			if (index2 >= group_size)
			{
				index1++;	//Group index
				index2=0;	//Subgroup index
			}
			row=0;
		}
		
		prevLoc=0;
		size_t string_size=line.size();
		
		currentLoc=line.find(" ", prevLoc);

		if (currentLoc < string_size)
		{
			col=0;
			while (currentLoc < string_size)	//Loop through each column
			{	
				temp=line.substr(prevLoc, currentLoc-prevLoc);
				size_t temp_size=temp.size();
				
				commaLoc=temp.find(",", 0);
				
				std::istringstream stuff(temp.substr(0, commaLoc));					//Real number
				double x;
				stuff>>x;
				
				double y;
				if (commaLoc <= temp_size)	//Check if a comma was located
				{
					std::istringstream stuff2(temp.substr(commaLoc+1, temp_size));  //Imaginary number
					stuff2>>y;
				}
				else
				{
					y=0;
				}

				Matrix[index1][index2](row, col)=std::complex<double>(x,y);
				col++;
				
				prevLoc=currentLoc; 
				currentLoc=line.find(" ", prevLoc+1);
			}
			if (verbose==yes)
			{
				std::cout << "--------------------Data Loaded-----------------------------" << std::endl;
				USs::OutputMatrix(Matrix[index1][index2]);
				std::cout << "------------------------------------------------------------" << std::endl;
			}
			row++;
		}else{
			if (row>0){index1++;}
			index2=0;
			row=0;
		}
	} 
	myfile.close();

	return Matrix;
}
void USs::MakeGnuplotScript(size_t nlines, char* filename, char* xlabel, char* ylabel, char* datafile, bool islog)
{
	char tempfilename[80];
	ofstream dataout;	
	if(filename!=NULL)
			UFs::OpenFile(strcat(strcpy(tempfilename,filename), ".p"),dataout, 16);
	else return;

	dataout<<"set term postscript eps enhanced color blacktext 'Helvetica' 24\n";
	dataout<<"set terminal postscript\n"; 
	dataout<<"set output '| epstopdf --filter >" << filename << ".pdf'\n";
	dataout<<"set ylabel 'Controls (GHz)'\n";
	dataout<<"set xlabel 'Time during Gate (ns)'\n";
	
	if (islog)
		dataout<<"set log y\n";
	else dataout<<"unset log y\n";	
	
//	dataout<<"set xr [1:20]\n";
//	dataout<<"set yr [-5:5]\n";
	dataout<<"set nokey\n";
	dataout<<"unset title\n";



	dataout<<"plot ";
	
	for(size_t j=2; j<nlines+2; j++)
	{	dataout << "'" << datafile << "' u 1:" << j << "w l";
		if(j<nlines+1) dataout << ",";
	}
	
	dataout<<"\nset output\n";
	dataout<<"quit\n";
	
	dataout.close();


}

#endif /* UsefullFunctions_h */

