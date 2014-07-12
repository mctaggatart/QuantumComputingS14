set term postscript eps color blacktext "Helvetica" 24
set output '2chcomparedelta_5ns.eps'
set log y
set ylabel 'Error'
set xlabel 'Delta_0 (GHz/2\pi)'
# set xr [0.0:0.022]
set title "Leakage error for 2 channels \Delta_2=-2.5"
plot 'delta5compareDrag.dat' w l,  'delta5compareDragL.dat' w l, 'delta5compareAdi.dat' w l, 'delta5compareGauss.dat' w l
set output
quit
