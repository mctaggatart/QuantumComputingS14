set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >CoherentCavity_cd_pi_10ns_1g_2D_0.01leakpen/popul_final_0.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'CoherentCavity_cd_pi_10ns_1g_2D_0.01leakpen/popul_final_0.dat' u 1:2w l,'CoherentCavity_cd_pi_10ns_1g_2D_0.01leakpen/popul_final_0.dat' u 1:3w l,'CoherentCavity_cd_pi_10ns_1g_2D_0.01leakpen/popul_final_0.dat' u 1:4w l,'CoherentCavity_cd_pi_10ns_1g_2D_0.01leakpen/popul_final_0.dat' u 1:5w l,'CoherentCavity_cd_pi_10ns_1g_2D_0.01leakpen/popul_final_0.dat' u 1:6w l,'CoherentCavity_cd_pi_10ns_1g_2D_0.01leakpen/popul_final_0.dat' u 1:7w l,'CoherentCavity_cd_pi_10ns_1g_2D_0.01leakpen/popul_final_0.dat' u 1:8w l,'CoherentCavity_cd_pi_10ns_1g_2D_0.01leakpen/popul_final_0.dat' u 1:9w l
set output
quit
