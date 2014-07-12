set term postscript eps color blacktext "Helvetica" 24
set output 'Drg_2chan_2.5.eps'
set log y
set ylabel 'Error'
set xlabel 'Gate Time (ns)'
# set xr [0.0:0.022]
set title "1-2 transition leakage error w anharm=-2.5"
plot 'Ga2_2chan_2.5.dat' u 1:2 title 'Gauss a=2' w l, 'Adi_2ch_5l_2.5_c.dat' u 1:2 title 'G+Der1+Stark a=2' w l, 'Ga3_Adi_2ch_5l_2.5_03.dat' u 1:2 title 'G+Stark + O3 a=3' w l, 'Ga4_Adi_2ch_5l_2.5_03.dat' u 1:2 title 'G+Stark + O3 a=4' w l, 'Ga5_Adi_2ch_5l_2.5_03.dat' u 1:2 title 'G+Stark + O3 a=5' w l
set output
quit
