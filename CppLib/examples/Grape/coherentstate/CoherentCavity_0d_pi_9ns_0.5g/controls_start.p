set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >CoherentCavity_0d_pi_9ns_0.5g/controls_start.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
set nokey
unset title
plot 'CoherentCavity_0d_pi_9ns_0.5g/controls_start.dat' u 1:2w l,'CoherentCavity_0d_pi_9ns_0.5g/controls_start.dat' u 1:3w l,'CoherentCavity_0d_pi_9ns_0.5g/controls_start.dat' u 1:4w l
set output
quit
