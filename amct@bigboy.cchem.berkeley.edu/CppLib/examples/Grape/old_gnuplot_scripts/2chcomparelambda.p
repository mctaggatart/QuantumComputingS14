set term postscript eps enhanced color blacktext "Helvetica" 24
set output '2chcomparelambda5ns_pps.eps'
set log y
set ylabel 'Error'
set xlabel '{/Symbol l}_0 ({/Symbol l}_2)'
# set xr [0.0:0.022]
set nokey
#set title "Leakage error for Harmonic 2 channels {/Symbol D}_2=-2.5, T_g=5ns"
unset title
#plot 'lamb5nsHarm2.5_G3.dat' title 'G3' w l, 'lamb5nsHarm2.5_drgC_G2.dat' title 'C' w l, 'lamb5nsHarm2.5_drgC2_G2.dat' title 'C2' w l, 'lamb5nsHarm2.5_drgB2_G3.dat' title 'B' w l, 'lamb5nsHarm2.5_drgA2_G3.dat'  title 'A' w l
plot 'lamb5nsHarm2.5_G3.dat' title 'G3' linewidth 6 w l, 'lamb5nsHarm2.5_drgC_G2.dat' title 'C' linewidth 6 w l,'lamb5nsHarm2.5_drgB2_G3.dat' title 'B' linewidth 6 w l
set output
quit
