set term postscript eps enhanced color blacktext "Helvetica" 24
set output "| epstopdf --filter >rwavalidity.pdf"
set log y
set ylabel "Gate error (1- {/Symbol F})"
set xlabel "Gate time (2{/Symbol p}/{/Symbol w}_1)"
set xr [2:6]
set yr [0.0000001:1]
set nokey
unset title
unset arrow

plot 'RWAcompare_sweepconfigs_old.dat' u ($1*2):3 w l lt 3, 'RWAcompare_sweepconfigs_old.dat' u ($1*2):4 w l lt 3, 'RWAcompare_sweepconfigs_old.dat' u ($1*2):5 w l lt 3, 'RWAcompare_sweepconfigs_old.dat' u ($1*2):6 w l lt 3, 'RWAcompare_sweepconfigs_old.dat' u ($1*2):7 w l lt 3, 'RWAcompare_sweepconfigs_old.dat' u ($1*2):8 w l lt 3, 'RWAcompare_sweepconfigs_old.dat' u ($1*2):9 w l lt 3, 'RWAcompare_sweepconfigs_old.dat' u ($1*2):10 w l lt 3, 'RWAcompare_sweepconfigs_old.dat' u ($1*2):11 w l lt 3, 'RWAcompare2_sweepconfigs_old.dat' u ($1*2):3 w l lt 2, 'RWAcompare2_sweepconfigs_old.dat' u ($1*2):4 w l lt 2, 'RWAcompare2_sweepconfigs_old.dat' u ($1*2):5 w l lt 2, 'RWAcompare2_sweepconfigs_old.dat' u ($1*2):6 w l lt 2, 'RWAcompare2_sweepconfigs_old.dat' u ($1*2):7 w l lt 2, 'RWAcompare2_sweepconfigs_old.dat' u ($1*2):8 w l lt 2, 'RWAcompare2_sweepconfigs_old.dat' u ($1*2):9 w l lt 2, 'RWAcompare2_sweepconfigs_old.dat' u ($1*2):10 w l lt 2, 'RWAcompare2_sweepconfigs_old.dat' u ($1*2):11 w l lt 2, 'RWAcompare_sweepconfigs_old.dat' u ($1*2):2 w l lw 3 lt 4


set output
quit



