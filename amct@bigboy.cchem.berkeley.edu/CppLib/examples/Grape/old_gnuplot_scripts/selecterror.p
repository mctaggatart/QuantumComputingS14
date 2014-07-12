set term postscript eps color blacktext "Helvetica" 24
set output 'selecterror_grape.eps'
set log y
set nokey
set ylabel 'Gate Error'
set xlabel 'Gate Time(ns)'
set xr[0:40]

#set size 0.75,0.75

plot 'sel_gaussian.dat' u 1:2  title  'gauss' linewidth 6 w l,  'sel_drag2.dat' u 1:2 title 'drag' linewidth 6 w l,'sel_grape.dat' u 1:2  title  'grape 1try' linewidth 6 w l

set output
quit
