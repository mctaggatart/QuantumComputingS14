set term postscript eps color blacktext "Helvetica" 24
set output 'Dr03_2.5_levels.eps'
set log y
set ylabel 'Error'
set xlabel 'Gate Time (ns)'
# set xr [0.0:0.022]
set title "Leakage error with anharm. \Delta=-2.5"
plot 'Dr03_3L_2.5.dat' u 1:2 title 'Drag 1 leakage above' w l, 'Dr03_4L_2.5.dat' u 1:2 title 'Drag 2 leakages above' w l, 'Dr03_5L_2.5.dat' u 1:2 title 'Drag 3 leakages above' w l, 'Ga2_3L_2.5.dat' u 1:2 title 'Gauss 1 leakage above' w l, 'Ga2_4L_2.5.dat' u 1:2 title 'Gauss 2 leakages above' w l
set output
quit
