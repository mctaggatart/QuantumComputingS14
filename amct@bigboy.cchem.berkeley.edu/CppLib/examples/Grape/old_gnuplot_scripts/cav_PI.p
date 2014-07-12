set term postscript eps color blacktext "Helvetica" 24
set output 'cav_PI.eps'
set log y
set xr [0.0:30]
set ylabel 'Gate Error'
set xlabel 'Gate Time (ns)'
#set title "Qubit3L+Cavity3L Delta=1.5GHz anharm=330MHz "
#plot 'cav_PI_tan0.9ns.dat' title 'Tan1ns + const d' w l, 'cav_PI_tanD5ns.dat' title 'Tan5ns +d(t) + Deriv' w l,'cav_PI_gauss1.5ns.dat' title 'Gauss1.5ns + d(t)' w l, 'cav_PI_tan4ns.dat' title 'Tan4ns + d(t)' w l

set size 0.7, 0.7

set nokey

plot 'cav_PI_tan0.9ns.dat' title 'Square' linewidth 6 w l, 'cav_PI_tanD5ns.dat' title 'Drag' linewidth 6 w l,'cav_PI_gauss1.5ns.dat' title 'Stark Shift' linewidth 6 w l
set output
quit