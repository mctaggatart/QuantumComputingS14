set term postscript eps color blacktext "Helvetica" 24
set output '2chan_2.5.eps'
set log y
set ylabel 'Error'
set xlabel 'Gate Time (ns)'
# set xr [0.0:0.022]
set title "1-2 transition leakage error w anharm=-2.5"
plot 'Dr03_2chan_2.5.dat' u 1:2 title 'Drag O3' w l, 'Dr01_2chan_2.5.dat' u 1:2 title 'Drag 01' w l, 'Dr0_2chan_2.5.dat' u 1:2 title 'Gauss + Deriv' w l, 'Ga2_2chan_2.5.dat' u 1:2 title 'Gauss a=2' w l, 'Grap_2ch_XYZ_2.5.dat' u 1:2 title 'Grape XYZ' w l,'Grap_2ch_XY_2.5.dat' u 1:2 title 'Grap XY' w l, 'Grap_2ch_X_2.5.dat' u 1:2 title 'Grp X' w l
set output
quit
