set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >OffResCoupling/controls_int.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'OffResCoupling/controls_int.dat' u 1:2w l,'OffResCoupling/controls_int.dat' u 1:3w l,'OffResCoupling/controls_int.dat' u 1:4w l,'OffResCoupling/controls_int.dat' u 1:5w l,'OffResCoupling/controls_int.dat' u 1:6w l
set output
quit
