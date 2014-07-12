set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >Unitary4/controls_int.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'Unitary4/controls_int.dat' u 1:2w l,'Unitary4/controls_int.dat' u 1:3w l
set output
quit