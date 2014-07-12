set term postscript eps color blacktext "Helvetica" 24
set output 'dragtrigerrorladder.eps'
set log y
set ylabel 'Gate Error'
set xlabel 'Gate Time(ns)'
set nokey
set xr[1:5]
set yr[0.000001:1]

#set size 0.75,0.75
plot 'gauss3lev.dat' u 1:2  title  'gauss' linewidth 6 w l,  'drag5lev.dat' u 1:2 title 'drag (5 levels)' linewidth 6 w l,'dragtrig5lev.dat' u 1:2 title 'drag trig (5 levels)' linewidth 6 w l,  'grape02linladder.dat' u 1:2 title 'grape 4controls' linewidth 6 w l


set output
quit
