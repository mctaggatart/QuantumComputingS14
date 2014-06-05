set term postscript eps enhanced color blacktext "Helvetica" 24
set output "| epstopdf --filter > pixelvalidity.pdf"
set log y
set ylabel "Gate error(1-{/Symbol F)"
set xlabel "Gate time (ns)"
set nokey
set xr[3:8]
set yr[0.00000000000005:1]
unset title
unset arrow
set arrow from 4,0.00000000000005 to 4, 1 nohead lt 7

plot 'EnvelopeCompare_test_sweeptimes.dat' u 1:3 w l lt 5 lw 3, 'EnvelopeCompare_test_sweeptimes.dat' u 1:5 w l lt 1 lw 3, 'EnvelopeCompare_test_sweeptimes.dat' u 1:6 w l lt 3 lw 3,  'EnvelopeCompare_test_sweeptimes.dat' u 1:7 w l lt 2 lw 3,  'EnvelopeCompare_test_sweeptimes.dat' u 1:8 w l lt 4 lw 3 

set output
quit

