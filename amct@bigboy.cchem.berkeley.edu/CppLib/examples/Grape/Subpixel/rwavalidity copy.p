set term postscript eps enhanced color blacktext "Helvetica" 24
set output 'rwavalidity.eps'
set log y
set ylabel "Gate Error"
set xlabel "RWA validity (Gate Time * Qubit Frequency)"
 set xr [4:14]
set nokey
unset title

plot 'RWAcompare_sweeptimes.dat' u 1:3 w l, 'RWAcompare_sweeptimes.dat' u 1:4 w l, 'RWAcompare_sweeptimes.dat' u 1:6 w l, 'RWAcompare_sweeptimes.dat' u 1:9 w l, 'RWAcompare_sweeptimes.dat' u 1:10 w l, 'RWAcompare_sweeptimes.dat' u 1:11 w l, 'RWAcompare_sweeptimes.dat' u 1:13 w l, 'RWAcompare_sweeptimes.dat' u 1:14 w l, 'RWAcompare_sweeptimes.dat' u 1:15 w l

set output
quit

