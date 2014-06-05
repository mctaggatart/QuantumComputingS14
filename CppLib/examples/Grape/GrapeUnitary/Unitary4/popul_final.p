set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >Unitary4/popul_final.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
set nokey
unset title
plot 'Unitary4/popul_final.dat' u 1:2w l,'Unitary4/popul_final.dat' u 1:3w l,'Unitary4/popul_final.dat' u 1:4w l
set output
quit
