set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >OffResCoupling_delme_1GHz/controls_start.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'OffResCoupling_delme_1GHz/controls_start.dat' u 1:2w l,'OffResCoupling_delme_1GHz/controls_start.dat' u 1:3w l,'OffResCoupling_delme_1GHz/controls_start.dat' u 1:4w l,'OffResCoupling_delme_1GHz/controls_start.dat' u 1:5w l,'OffResCoupling_delme_1GHz/controls_start.dat' u 1:6w l,'OffResCoupling_delme_1GHz/controls_start.dat' u 1:7w l
set output
quit
