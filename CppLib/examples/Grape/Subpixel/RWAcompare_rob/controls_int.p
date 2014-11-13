set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >RWAcompare_rob/controls_int.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'RWAcompare_rob/controls_int.dat' u 1:2w l,'RWAcompare_rob/controls_int.dat' u 1:3w l
set output
quit
