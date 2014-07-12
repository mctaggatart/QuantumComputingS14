set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'numericdrag.eps'
set log y
unset title
set nokey
set ylabel 'Gate Error 1-F_g'
set xlabel 'Gate Time (ns)'
set yr [0.0000001:1]
set xr [1:8]
#set title "Single leakage channel along anharmonic ladder"
plot 'timesSSavg2_sweeptime.dat' w l title 'optimized (X+ avg Z)', 'timesSSavg_sweeptime.dat' w l title 'X+ avg Z', 'gauss_sweeptime.dat' w l title 'Gaussian', '1DragO3_sweeptime.dat' w l title 'Drag Order 3', 'timesSS_sweeptime.dat' w l title ' X + Stark Shift Corr', 'dragB_sweeptime.dat' w l title 'half drag', 'dragB2_sweeptime.dat' w l title 'Opt X+derY', 'dragC2_sweeptime.dat' w l title 'Opt X+ derY+ avg Z', '1Drag2Phot_sweeptime.dat' w l title 'Opt X +derY +avgZ + 0-2'

#'L1_DRAGF.dat'u 1:2  title 'GRAPE optimized' w l
#, 'L1_DRAGG.dat'u 1:2  title 'GRAPE optimized' w l

set output
quit

