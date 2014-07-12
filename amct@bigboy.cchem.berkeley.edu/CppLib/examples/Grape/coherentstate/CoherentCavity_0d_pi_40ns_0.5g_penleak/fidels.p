set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >CoherentCavity_0d_pi_40ns_0.5g_penleak/fidels.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
set log y
set nokey
unset title
plot 'CoherentCavity_0d_pi_40ns_0.5g_penleak/fidels.dat' u 1:2w l
set output
quit
