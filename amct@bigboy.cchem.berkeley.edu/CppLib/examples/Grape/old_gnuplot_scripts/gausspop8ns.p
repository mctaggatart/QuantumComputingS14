set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'popgauss8ns.eps'
unset log y
set ylabel 'Populations'
set xlabel 'Time during Pulse(ns)'
set nokey
set xr [0:8]
#plot 'gausspop8ns.dat' u 1:2 title '|1>' w l, 'gausspop8ns.dat' u 1:3 title 'in phase|0>' w l, 'gausspop8ns.dat' u 1:4 title 'off-phase |0>' w l, 'gausspop8ns.dat' u 1:5 title '|2>' w l
plot 'gausspop8ns.dat' u 1:2 title '|1>' linewidth 4 w l, 'gausspop8ns.dat' u 1:3 title 'in phase|0>' linewidth 4 w l,  'gausspop8ns.dat' u 1:5 title '|2>' linewidth 4 w l
set output
quit