set xlabel "eps_abs"
set logscale y
set logscale x
set key left top
set ylabel "error and h"
set samples 100

set yrange[1e-15:1]
	

plot "TestComplexRK45.dat" u 1:8 title 'Error', "TestComplexRK45.dat" u 1:2 title 'h'

