/*
Name: Grape.cpp
Author: Jay Gambetta

Dependences: GrapeUnitary.hpp  
Brief Discription: This demonstrates grape for unitary for PHI_4 on a four qubit system with hamlitonian 

\sum_i pi*\omega_i*Z_i + \sum_{i<j} pi*J_ij*(Z_iZ_j + X_iX_j + Y_iY_j)

# Nuclei: <frequency><tab><name><tab><type>
-2994.64	C1	C
-25511.36	C2	C
-21582.82	C3	C
-29442.59	C4	C
# Reference frequencies:
refFreq	C	C1

We pulse at -16000 to be in the middle of the spectrum so you have to add 16000 to all the chemical shift terms.

# Couplings:  <frequency><tab><name><tab><name> (must be in lex. order.)
  20.81	C1	C2
   0.73	C1	C3
   3.51	C1	C4
  34.83	C2	C3
   0.59	C2	C4
  36.08	C3	C4

The control Hamiltonians are
\omega_rf_x \sum_i X_i + \omega_rf_y \sum_i Y_i

and have a maximum value of 2*pi*16.7e3.


For this system my code can find single qubit 180 pulses is 0.5-1 minute with no robustness 
using 300 timesteps of 2us each.  I can find a CNOT gate between C3 and C4 with 600 timesteps of 
20us each in a couple minutes from a good guess but it may take several tries to get a good guess.  
This is on a 3GHz machine.

Version History
	v0: Feb  11th, 2009.
	
time ./GrapeUnitary3.e GrapeUnitary3.dat

real	0m59.024s
user	0m57.992s
sys	0m0.369s


*/
#include <OptimizeEvolution.hpp>
using namespace std;

int main (int argc, char const *argv[]){
	
	verbose=no;
	cout << "Running program " << argv[0] << endl;
	
	//Grape inputs
	size_t num_time=60, dim = 16,numDis=1, typeDis=2, num_controls =2;
	size_t max_iter=1000;
	double tolerance=std::numeric_limits<double>::min(), fidelity=0.99, base_a=2.0, epsilon=5000, dt=0.00002, tgate=dt*double(num_time);

 matrix<complex<double> >* dis;
	dis= new matrix<std::complex<double> >[numDis*typeDis];
	for(int k=0; k<typeDis*numDis; ++k){
	  dis[k].initialize(dim,dim); //NEED TO MAKE A A GLOBAL VARIABLE SET!!
}	
		MOs::Destroy(dis[0]);
	dis[1]=(MOs::Dagger(dis[0]))*(dis[0]);
 
	OptimizeEvolution sys(dim, num_controls, num_time, dt, dis, numDis, typeDis,"Unitary3");
	sys.SetNumericalParameters(fidelity=0.99, base_a, epsilon, tolerance, max_iter);
	//EDITTT
	sys.SetOppsDesired(dis);	
	matrix<complex<double> > sx(2,2), sy(2,2), sz(2,2), si(2,2), Hdrift(dim,dim), Hcontrol(dim,dim);
	MOs::Pauli(sx, sy, sz);
	MOs::Identity(si);
	
	matrix<complex<double> > Z1(MOs::TensorProduct(sz,MOs::TensorProduct(si,MOs::TensorProduct(si,si))));
	matrix<complex<double> > Z2(MOs::TensorProduct(si,MOs::TensorProduct(sz,MOs::TensorProduct(si,si))));
	matrix<complex<double> > Z3(MOs::TensorProduct(si,MOs::TensorProduct(si,MOs::TensorProduct(sz,si))));
	matrix<complex<double> > Z4(MOs::TensorProduct(si,MOs::TensorProduct(si,MOs::TensorProduct(si,sz))));
	matrix<complex<double> > X1(MOs::TensorProduct(sx,MOs::TensorProduct(si,MOs::TensorProduct(si,si))));
	matrix<complex<double> > X2(MOs::TensorProduct(si,MOs::TensorProduct(sx,MOs::TensorProduct(si,si))));
	matrix<complex<double> > X3(MOs::TensorProduct(si,MOs::TensorProduct(si,MOs::TensorProduct(sx,si))));
	matrix<complex<double> > X4(MOs::TensorProduct(si,MOs::TensorProduct(si,MOs::TensorProduct(si,sx))));
	matrix<complex<double> > Y1(MOs::TensorProduct(sy,MOs::TensorProduct(si,MOs::TensorProduct(si,si))));
	matrix<complex<double> > Y2(MOs::TensorProduct(si,MOs::TensorProduct(sy,MOs::TensorProduct(si,si))));
	matrix<complex<double> > Y3(MOs::TensorProduct(si,MOs::TensorProduct(si,MOs::TensorProduct(sy,si))));
	matrix<complex<double> > Y4(MOs::TensorProduct(si,MOs::TensorProduct(si,MOs::TensorProduct(si,sy))));
	
	double omega1=-2994.64+16000;
	double omega2=-25511.36+16000;
	double omega3=-21582.82+16000;
	double omega4=-29442.59+16000;
	double J12 = 20.81;
	double J13 = 0.73;
	double J14 = 3.51;
	double J23 = 34.83;
	double J24 = 0.59;
	double J34 = 36.08;
	
	
	//Drift
	Hdrift= M_PI*(omega1*Z1+omega2*Z2+omega3*Z3+omega4*Z4)+M_PI*(J12*(Z1*Z2 + X1*X2 + Y1*Y2) +J13*(Z1*Z3 + X1*X3 + Y1*Y3)+J14*(Z1*Z4 + X1*X4 + Y1*Y4) +J23*(Z2*Z3 + X2*X3 + Y2*Y3) + J24*(Z2*Z4 + X2*X4 + Y2*Y4) + J34*(Z3*Z4 + X3*X4 + Y3*Y4));
	//	cout <<Hdrift<<"hdrift"<<endl;
	sys.SetHdrift(Hdrift);
	
	Control u0(M_PI*(X1+X2+X3+X4), dt, num_time, 1, &sys, "u2_control0");	
	u0.RandomControl(-1000, 1000);
	
	Control u1(M_PI*(Y1+Y2+Y3+Y4), dt, num_time, 1, &sys, "u2_control1");
	u1.RandomControl(-1000, 1000);
	
	matrix<complex<double> > U_desired(dim,dim), CNOT12(dim,dim), cnot(4,4), goal(dim,dim);
	MOs::Null(cnot);
	//goal.Initalize();
	cnot(0,0)=complex<double>(1,0);
	cnot(1,1)=complex<double>(1,0);
	cnot(2,3)=complex<double>(1,0);
	cnot(3,2)=complex<double>(1,0);
	CNOT12=MOs::TensorProduct(si,MOs::TensorProduct(si,cnot));
	goal(15,15)=complex<double>(1,0);
	sys.SetOmega(5.0);
	sys.SetRhoInitial(goal);
	U_desired=CNOT12; 
	sys.SetUDesired(U_desired);
	//cout<<U_desired<<"U Des"<<endl;
	sys.SetTrueRhoDesired(U_desired);
	//run grape
	sys.UnitaryTransfer();
	
	return 0;
}
