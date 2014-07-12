set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >CoherentCavity_cd_pi_10ns_0.5g_4D/controls_start.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'CoherentCavity_cd_pi_10ns_0.5g_4D/controls_start.dat' u 1:2w l,'CoherentCavity_cd_pi_10ns_0.5g_4D/controls_start.dat' u 1:3w l,'CoherentCavity_cd_pi_10ns_0.5g_4D/controls_start.dat' u 1:4w l
set output
quit
