set term postscript eps color blacktext "Helvetica" 24
set output 'cav_q2_PI.eps'
set log y
# set xr [0.0:0.022]
set ylabel '1 - ||U*Unot||^2'
set xlabel 'Gate Time (ns)'
set title "Qubit2L+Cavity2L Delta=1.5GHz anharm=330MHz "
plot 'cav_PI_th_gauss1.5.dat' title 'Gauss 1.5*sigma + d(x^2)' w l, 'cav_PI_th2_gauss1.5.dat' title 'Gauss 1.5*sigma + d(x^4)' w l, 'cav_PI_th2b_gauss1.5.dat' title 'Gauss 1.5*sigma + d(x^2) +deriv' w l, 'cav_PI_th_square1ns.dat' title 'Square 1ns rise + const d' w l
set output
quit