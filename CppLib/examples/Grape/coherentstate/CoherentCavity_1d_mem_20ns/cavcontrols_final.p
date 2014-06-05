set term postscript eps enhanced color blacktext "Helvetica" 24
set terminal postscript 
set output "| epstopdf --filter >cavcontrols_final.pdf"
set ylabel "Controls (GHz)"
set xlabel "Time during Gate (ns)"
unset log y
set xr [1:20]
set yr [-5:5]
set nokey
unset title

plot 'cavcontrol0_final.dat' w l, 'cavcontrol1_final.dat' w l,'cavcontrol2_final.dat' w l


set output
quit

