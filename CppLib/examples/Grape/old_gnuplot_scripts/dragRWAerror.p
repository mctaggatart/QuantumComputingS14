set term postscript eps color blacktext "Helvetica" 24
set output 'dragRWAerror.eps'
set log y
set ylabel 'Error'
set xlabel 'Gate Time (ns)'
# set xr [0.0:0.022]
set title "Rotating wave approx error with \omega = 4GHz*2PI"
#plot 'rwa_gauss_24.dat' u 1:2 w l,  'rwa_drag_24.dat' u 1:2 w l, 'rwa_dragO3_24.dat' u 1:2 w l, 'rwa_stark_24.dat' u 1:2 w l, 'rwa_deriv_24.dat' u 1:2 w l, 'rwa_stark03_24.dat' u 1:2 w l, 'rwa_derivO3_24.dat' u 1:2 w l, 'rwa_avgstark_24.dat' u 1:2 w l
plot 'rwa_gauss_24.dat' u 1:2 title  'gauss' w l, 'rwa_dragO3_24.dat' u 1:2  title 'drag full order' w l, 'rwa_deriv_24.dat' u 1:2 title 'deriv only'  w l, 'rwa_avgstark_24.dat'u 1:2  title 'const det' w l
set output
quit
