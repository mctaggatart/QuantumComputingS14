set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'sweepsubpixels.eps'
set log y
set ylabel "error"
set xlabel "number of subpixels"
# set xr [0.0:0.022]
set nokey
#set title "error vs subpixels"
unset title

plot 'com_evolrand_frq10_Hd5_Hc3.dat' u 1:2 w l, 'com_evolrand_frq10_Hd5_Hc3.dat' u 1:3 w l

set output
quit

