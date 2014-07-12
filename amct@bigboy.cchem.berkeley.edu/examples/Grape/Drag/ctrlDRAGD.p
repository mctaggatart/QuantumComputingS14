set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'ctrlDRAGD.eps'
set nokey
unset log y
set ylabel 'Control Amplitude (2pi GHz)'
set xlabel 'Time during Pulse(ns)'
plot 'ctrlDRAGD.dat' u 1:2 title 'Ex' linewidth 6 w l, 'ctrlDRAGD.dat' u 1:3 title 'Ey' linewidth 6 w l, 'ctrlDRAGD.dat' u 1:4 title '{/Symbol d}_1' linewidth 6 w l
set output
quit
