set term postscript eps enhanced color blacktext "Helvetica" 24
set terminal postscript 
set output "| epstopdf --filter >offreserrorvssubpixels_delme.pdf"
set log y
set ylabel 'Gate error(1-{/Symbol F)'
set xlabel "Number of subpixels"
set xr [1:100]
set yr [0.000005:1]
set nokey
unset title
unset arrow

plot 'OffResCoupling_test_grainconfigs.dat' w l lw 5, 'OffResCoupling_05GHZ_grainconfigs.dat' w l lw 5, 'OffResCoupling_delme_grainconfigs.dat' w l lw 5

set output
quit

