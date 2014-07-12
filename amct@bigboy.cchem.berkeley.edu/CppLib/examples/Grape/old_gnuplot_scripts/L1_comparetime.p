set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'L1_comparetime.eps'
set log y
unset title
set nokey
set ylabel 'Gate Error 1-F_g'
set xlabel 'Gate Time (ns)'
set yr [0.0000001:1]
set xr [2:9]
#set title "Single leakage channel along anharmonic ladder"
plot 'L1_gauss1.5.dat' u 1:2 title  'gauss {/Symbol s}=t_g/1.5' linewidth 6 w l, 'L1_DRAGA.dat' u 1:2  title 'adiabatic' linewidth 6 w l, 'L1_DRAGB.dat' u 1:2 title 'deriv only' linewidth 6  w l, 'L1_DRAGC.dat'u 1:2  title 'DRAG 3 orders' linewidth 6 w l, 'grapelinladder.dat'u 1:2  title 'GRAPE optimized' linewidth 6 w l

#'L1_DRAGF.dat'u 1:2  title 'GRAPE optimized' w l
#, 'L1_DRAGG.dat'u 1:2  title 'GRAPE optimized' w l

set output
quit

