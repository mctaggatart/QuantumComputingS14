set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >Unitary4/fidels.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
set log y
set nokey
unset title
plot 'Unitary4/fidels.dat' u 1:2w l,'Unitary4/fidels.dat' u 1:3w l
set output
quit
