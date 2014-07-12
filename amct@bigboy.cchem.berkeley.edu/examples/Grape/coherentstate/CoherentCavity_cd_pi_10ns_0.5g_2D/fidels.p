set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >CoherentCavity_cd_pi_10ns_0.5g_2D/fidels.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
set log y
set nokey
unset title
plot 'CoherentCavity_cd_pi_10ns_0.5g_2D/fidels.dat' u 1:2w l,'CoherentCavity_cd_pi_10ns_0.5g_2D/fidels.dat' u 1:3w l
set output
quit
