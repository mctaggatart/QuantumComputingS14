set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'envelope.eps'
unset log y
set ylabel "amplitude (arb units)"
set xlabel "time (ns)"
# set xr [0.0:0.022]
set nokey
#set title "Leakage error for Harmonic 2 channels {/Symbol D}_2=-2.5, T_g=5ns"
unset title
plot 'sampled0_pixels.dat' w l, 'truecontrol0_gaussfilt.dat' w l, 'sampled0_cubinterpol.dat' w l
set output
quit

