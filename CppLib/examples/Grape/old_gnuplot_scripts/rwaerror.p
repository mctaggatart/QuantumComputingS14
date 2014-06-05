set term postscript eps color blacktext "Helvetica" 24
set output 'rwaerror.eps'
set log y

set nokey
set ylabel 'Error'
set xlabel 'Gate Time(number of qubit oscillations)'
set xr[0:20]

#plot 'rwa_gauss_24.dat' u 1:2 title  'gauss' w l, 'rwa_dragO3_24.dat' u 1:2  title 'drag full order' w l, 'rwa_deriv_24.dat' u 1:2 title 'deriv only'  w l, 'rwa_avgstark_24.dat' u 1:2  title 'const det' w l

#plot 'rwadet.dat' title 'const det' w l, 'rwa.dat' w l, 'rwaadb.dat' w l, 'rwa4alphdrag.dat' w l, 'rwa1.5drag.dat' w l, 'rwa1.5dragb.dat' w l



#set size 0.75,0.75

#plot 'rwa_gauss_24.dat' u ($1*4):2  title  'gauss' linewidth 6 w l,  'rwa1.5drag.dat' u ($1*4):2 title 'drag' linewidth 6 w l,  'rwaerror.dat' u ($1*4):5 title 'grape' linewidth 6 w l

plot 'rwa_gauss_24.dat' u ($1*4):2  title  'gauss'linewidth 6  w l, 'rwadet.dat' u ($1*4):2  title 'const det'linewidth 6  w l,  'rwa1.5drag.dat' u ($1*4):2 title 'drag full order'linewidth 6  w l, 'rwa1.5dragb.dat' u ($1*4):2 title 'deriv only'linewidth 6  w l,  'rwagrapeerr.dat' u ($1*4):2 title 'grape' linewidth 6 w l,  'rwagrape.dat' u ($1*4):2 title 'grape' linewidth 6 w l


set output
quit
