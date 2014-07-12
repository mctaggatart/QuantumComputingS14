set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal png
set output 'fidels.png'
set xlabel 'Time during Gate (ns)'
set log y
set nokey
unset title
plot 'fidels.dat' u 1:2 w l lw 4,'fidels.dat' u 1:3 w l lw 4
set output
quit
