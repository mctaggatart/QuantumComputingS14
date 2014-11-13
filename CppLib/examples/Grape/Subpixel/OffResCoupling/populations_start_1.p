set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >OffResCoupling/populations_start_1.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'OffResCoupling/populations_start_1.dat' u 1:2w l,'OffResCoupling/populations_start_1.dat' u 1:3w l,'OffResCoupling/populations_start_1.dat' u 1:4w l,'OffResCoupling/populations_start_1.dat' u 1:5w l
set output
quit
