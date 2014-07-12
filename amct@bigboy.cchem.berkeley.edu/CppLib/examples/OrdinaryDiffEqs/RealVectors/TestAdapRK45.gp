set xlabel "eps_abs"
set logscale y
set logscale x
set key left top
set ylabel "error and h"
set samples 100

set yrange[1e-17:1]
	

plot "TestAdapRK45.dat" u 1:6 title 'Error', "TestAdapRK45.dat" u 1:2 title 'h'

