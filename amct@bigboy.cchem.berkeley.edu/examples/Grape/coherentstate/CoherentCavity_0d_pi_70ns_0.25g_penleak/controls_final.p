set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >CoherentCavity_0d_pi_70ns_0.25g_penleak/controls_final.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'CoherentCavity_0d_pi_70ns_0.25g_penleak/controls_final.dat' u 1:2w l,'CoherentCavity_0d_pi_70ns_0.25g_penleak/controls_final.dat' u 1:3w l,'CoherentCavity_0d_pi_70ns_0.25g_penleak/controls_final.dat' u 1:4w l
set output
quit
