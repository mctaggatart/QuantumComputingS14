set xlabel "h"
set logscale y
set logscale x
set key left top
set ylabel "error"
set samples 100

set yrange[1e-15:1]
	
# The 2.22045e-16 is the machine relative error

plot "TestEuler.dat" u 1:5 title 'Euler', "TestRombergEuler.dat" u 1:5 title 'Euler Romberg', 0.199526*x title '0.199526h', 0.251189*x**2 title '0.251189h^2', 6.30957e-18*x**-1 title '6.30957e-18/h', 2.22045e-16*x**-1 title 'numeric_limits<double>::epsilon()/h'
