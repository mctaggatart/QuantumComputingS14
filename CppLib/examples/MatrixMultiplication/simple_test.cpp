#include "MatrixExponential.hpp"
#include "QuantumOperations.hpp"

using namespace std;


int main (int argc, char const *argv[]) {

	verbose=no;
	cout << "Running program " << argv[0] << endl;
	
	matrix<double> a(3,3), b(3,3), c(3,3);
	a.SetOutputStyle(Matrix);
	b.SetOutputStyle(Matrix);
	c.SetOutputStyle(Matrix);
	
	a(0,0) = 1.0;
	a(1,1) = 2.0;
	a(2,2) = 3.0;
	
	cout << a << endl;
	
	b(0,1) = 1.0;
	b(1,2) = 2.0;
	b(2,0) = 3.0;
	
	cout << b << endl;
	
	c = a*b;
	
	cout << c << endl;

	return 0;
}
