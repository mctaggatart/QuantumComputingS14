set term postscript eps enhanced color blacktext "Helvetica" 24
set terminal postscript 
set output "| epstopdf --filter >offrescontrols.pdf"
set ylabel "Controls (GHz)"
set xlabel "Time during Gate (ns)"
unset log y
set xr [1:20]
set yr [-3:3]
set nokey
unset title

plot 'splited0_final.dat' w l, 'splited1_final.dat' w l,'splited2_final.dat' w l,'splited3_final.dat' w l,'splited4_final.dat' w l


set output
quit

