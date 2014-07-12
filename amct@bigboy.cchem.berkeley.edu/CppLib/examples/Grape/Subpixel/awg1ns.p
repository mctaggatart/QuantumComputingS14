set term postscript eps enhanced color blacktext "Helvetica" 24
set output "| epstopdf --filter > awg1ns.pdf"
unset log
set ylabel "Control amplitude(V)"
set xlabel "Time during gate"
set nokey
set xr[0:7]
set yr[-0.1:0.2]
unset title
unset arrow

plot 'T6_filter.dat' u 1:2 w l, 'T6_pix.dat' u 1:2 w l, 'T6_spline.dat' u 1:2 w l, 'T6.dat' u 1:2 w l

set output
quit

