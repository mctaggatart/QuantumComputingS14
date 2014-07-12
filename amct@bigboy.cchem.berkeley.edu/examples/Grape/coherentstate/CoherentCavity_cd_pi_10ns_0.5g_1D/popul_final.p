set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal png
set output 'popul_final.png'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'popul_final.dat' u 1:2 w l lw 4,'popul_final.dat' u 1:3 w l lw 4,'popul_final.dat' u 1:4 w l lw 4,'popul_final.dat' u 1:5 w l lw 4,'popul_final.dat' u 1:6 w l lw 4,'popul_final.dat' u 1:7 w l lw 4,'popul_final.dat' u 1:8 w l lw 4,'popul_final.dat' u 1:9 w l lw 4,'popul_final.dat' u 1:10 w l lw 4,'popul_final.dat' u 1:11 w l lw 4
set output
quit
