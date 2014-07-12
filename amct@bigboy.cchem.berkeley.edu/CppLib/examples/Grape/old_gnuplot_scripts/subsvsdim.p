set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'subsvsdim.eps'
unset log y
set ylabel "required subpixels"
set xlabel "dimension ()"
# set xr [0.0:0.022]
set nokey
unset title

plot 'subsvsdim.dat' u 1:2 w l, 'subsvsdim.dat' u 1:3 w l

set output
quit

