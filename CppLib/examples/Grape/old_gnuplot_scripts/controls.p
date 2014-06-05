set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'controls_gaussfilter_horns_buffer.eps'
unset log y
set ylabel "Frequency difference"
set xlabel "time during pulse (ns)"
# set xr [0.0:0.022]
set nokey
unset title

plot 'truecontrol_1_final.dat' w l, 'truecontrol_1_start.dat' w l

#'squarepulse.dat' w l

set output
quit

