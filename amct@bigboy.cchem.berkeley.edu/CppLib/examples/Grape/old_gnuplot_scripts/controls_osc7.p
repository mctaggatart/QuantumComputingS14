set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'controls_osc7.eps'
unset log y
set ylabel "Control Amplitude (V)"
set xlabel "time during pulse (ns)"
# set xr [0.0:0.022]
set nokey
unset title

plot 'T7_spline.dat' w l, 'T7.dat' w l, 'T7_filter.dat' w l , 'T7_pix.dat' w l
#'squarepulse.dat' w l

set output
quit

