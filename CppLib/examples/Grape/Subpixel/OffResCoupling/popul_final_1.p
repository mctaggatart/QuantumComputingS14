set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >OffResCoupling/popul_final_1.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'OffResCoupling/popul_final_1.dat' u 1:2w l,'OffResCoupling/popul_final_1.dat' u 1:3w l,'OffResCoupling/popul_final_1.dat' u 1:4w l,'OffResCoupling/popul_final_1.dat' u 1:5w l
set output
quit
