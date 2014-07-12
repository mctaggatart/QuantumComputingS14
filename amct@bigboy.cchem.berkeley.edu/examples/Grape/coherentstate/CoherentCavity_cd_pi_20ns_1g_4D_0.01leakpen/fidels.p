set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >CoherentCavity_cd_pi_20ns_1g_4D_0.01leakpen/fidels.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
set log y
set nokey
unset title
plot 'CoherentCavity_cd_pi_20ns_1g_4D_0.01leakpen/fidels.dat' u 1:2w l,'CoherentCavity_cd_pi_20ns_1g_4D_0.01leakpen/fidels.dat' u 1:3w l
set output
quit
