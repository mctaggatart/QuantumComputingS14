set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'ctrlDRAGE.eps'
set nokey
unset log y
set ylabel 'Control Amplitude (2pi GHz)'
set xlabel 'Time during Pulse(ns)'
plot 'ctrlDRAGE.dat' u 1:2 title 'Ex' linewidth 6 w l, 'ctrlDRAGE.dat' u 1:3 title 'Ey' linewidth 6 w l, 'ctrlDRAGE.dat' u 1:4 title '{/Symbol d}_1' linewidth 6 w l
set output
quit
