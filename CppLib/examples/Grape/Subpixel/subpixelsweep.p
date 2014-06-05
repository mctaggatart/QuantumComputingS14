set term postscript eps enhanced color blacktext "Helvetica" 24
set output "| epstopdf --filter > subpixelsweep_test.pdf"
set log y
set ylabel "Gate error(1-{/Symbol F)"
set xlabel "Number of subpixels"
set nokey
set xr[1:30]
set yr[0.00000000000005:1]
unset title
unset arrow

plot 'EnvelopeCompare_test_grainconfigs.dat' u 1:2 w l lt 5 lw 3, 'EnvelopeCompare_test_grainconfigs.dat' u 1:4 w l lt 1 lw 3, 'EnvelopeCompare_test_grainconfigs.dat' u 1:6 w l lt 2 lw 3, 'EnvelopeCompare_test_grainconfigs.dat' u 1:5 w l lt 4 lw 3

set output
quit

