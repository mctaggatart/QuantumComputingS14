set term postscript eps enhanced color blacktext "Helvetica" 24
set output '2chcomparetime.eps'
set log y
set nokey
# set xr [0.0:0.022]
#set title "Leakage error for 2 channels {/Symbol D}=0.4*2pi*GHz"
set ylabel 'Gate Error'
set xlabel 'Gate Time (ns)'
#plot 'drg1_2ch_5l_2.5.dat' w l, 'adi_2ch_5l_2.5.dat' w l,'drg2_2ch_5l_2.5.dat' w l
#plot '2ch_adiO3_g4.dat' w l,'2ch_adiO1_g3_2.5.dat' w l, '2ch_adiO3_g3.dat' w l, '2ch_gauss3_2.5.dat' w l, '2ch_gauss4_2.5.dat' w l, '2ch_drgO3b_g3_2.5.dat' w l, '2ch_drgO3b_g4_2.5.dat' w l, '2ch_drgO1c_g3_2.5.dat' w l

#plot '2ch_gauss3_2.5.dat'title 'Gauss G3 (t_{g}=3{/Symbol s})' w l, '2ch_drgO1c_g3_2.5.dat' title 'DRAG D3 (G3 + y=d_{t}x/4{/Symbol D}) '  w l, '2ch_drgO3b_g3_2.5.dat' title 'D3 + x correc' w l, '2ch_drgO3b_g4_2.5.dat' title 'D4 + x correc' w l, '2ch_adiO3_g3.dat' w l, '2ch_adiO3_g4.dat' title 'G4 + chirp + x corr' w l


#set size 0.7,0.7
#plot '2ch_gauss3_2.5.dat'title 'Gauss 3{/Symbol s}' linewidth 6 w l, '2ch_drgO3b_g4_2.5.dat' title 'Drag2' linewidth 6 w l

plot 'gauss2chan.dat' linewidth 6 w l, 'grape2chan.dat' linewidth 6 w l, 'halfdrag2chan.dat' linewidth 6 w l


#plot 'what.dat' w l, '2ch_adiO3_g3.dat' w l, '2ch_drgO3b_g3_2.5.dat' w l

set output
quit
