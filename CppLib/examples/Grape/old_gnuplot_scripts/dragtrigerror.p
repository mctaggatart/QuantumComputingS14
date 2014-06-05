set term postscript eps color blacktext "Helvetica" 24
set output 'dragtrigerror.eps'
set log y
set ylabel 'Gate Error'
set xlabel 'Gate Time(ns)'
set nokey
set xr[0:5]
set yr[0.00000000001:1]

#set size 0.75,0.75
plot 'gauss3lev.dat' u 1:2  title  'gauss' linewidth 6 w l,  'dragtrig3lev.dat' u 1:2 title 'drag trig (3 levels)' linewidth 6 w l,  'grape02lin3.dat' u 1:2 title 'grape 4controls' linewidth 6 w l

set output
quit
