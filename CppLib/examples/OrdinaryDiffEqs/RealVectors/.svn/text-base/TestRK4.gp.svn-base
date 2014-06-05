set xlabel "h"
set logscale y
set logscale x
set key left top
set ylabel "error"
set samples 100

set yrange[1e-17:1]
	
# The 2.22045e-16 is the machine relative error

plot "TestRK4.dat" u 1:5 title 'RK4', "TestRombergRK4.dat" u 1:5 title 'RK4 Romberg', 0.00316228*x**4 title '0.00316228h^4', 0.00316228*x**5 title '0.00316228h^5', 6.30957e-18*x**-1 title '6.30957e-18/h', 2.22045e-16*x**-1 title 'numeric_limits<double>::epsilon()/h'
