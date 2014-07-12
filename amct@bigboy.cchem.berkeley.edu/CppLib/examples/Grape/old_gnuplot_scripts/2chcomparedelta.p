set term postscript eps enhanced color blacktext "Helvetica" 24
set output '2chcomparedelta6ns.eps'
set log y
set ylabel 'Error'
set xlabel '{/Symbol D}_0 ({/Symbol D}_2)'
 set xr [-4:4]
#set title "Leakage error for Harmonic 2 channels {/Symbol D}_2=-2.5, T_g=6ns"
#plot 'delta8compareDrag.dat' w l,  'delta8compareDragL.dat' w l, 'delta8compareAdi.dat' w l, 'delta8compareGauss.dat' w l
#plot 'delta6nsG2.dat' w l title 'G2',  'delta6nsG2Ad3.dat' w l title 'Adi2O3', 'delta6nsG2D3.dat' w l title 'D2O3', 'delta6nsG2L1Dr3.dat' w l title 'L1Dr2O3', 'delta6nsG2L2Dr3.dat' w l title 'L2Dr2O3'
set nokey
set size 0.7,0.7
plot 'delta6nsG2.dat' linewidth 6 w l title 'G2', 'delta6nsG2L1Dr3.dat' linewidth 6  w l title 'L1Dr2O3', 'delta6nsG2L2Dr3.dat' linewidth 6  w l title 'L2Dr2O3'
set output
quit

apple_dore413@hotmail.com An Yen Geoyung

