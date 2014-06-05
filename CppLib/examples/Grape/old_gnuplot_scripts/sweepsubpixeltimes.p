set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'sweepsubpixeltimes.eps'
unset log y
set ylabel "total cpu time (ms)"
set xlabel "number of subpixels"
# set xr [0.0:0.022]
set nokey
unset title

plot 'times.dat' u 1:2 w l, 'times.dat' u 1:3 w l

set output
quit

