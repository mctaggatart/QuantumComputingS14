set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'ctrlDRAGA.eps'
unset log y
set ylabel 'Control Amplitude (2pi GHz)'
set xlabel 'Time during Pulse(ns)'
set nokey
plot 'ctrlDRAGA.dat' u 1:2 title 'Ex' linewidth 6 w l, 'ctrlDRAGA.dat' u 1:3 title 'Ey' linewidth 6 w l, 'ctrlDRAGA.dat' u 1:4 title '{/Symbol d}_1' linewidth 6 w l
set output
quit
