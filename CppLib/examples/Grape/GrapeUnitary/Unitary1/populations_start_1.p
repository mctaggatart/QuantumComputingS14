set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >Unitary1/populations_start_1.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'Unitary1/populations_start_1.dat' u 1:2w l,'Unitary1/populations_start_1.dat' u 1:3w l
set output
quit
