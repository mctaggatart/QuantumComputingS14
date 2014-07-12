set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal png
set output 'controls_final.png'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'controls_final.dat' u 1:2 w l lw 4
set output
quit
