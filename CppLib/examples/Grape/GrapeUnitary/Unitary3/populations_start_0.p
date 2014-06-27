set term postscript eps enhanced color blacktext 'Helvetica' 24
set terminal postscript
set output '| epstopdf --filter >Unitary3/populations_start_0.pdf'
set ylabel 'Controls (GHz)'
set xlabel 'Time during Gate (ns)'
unset log y
set nokey
unset title
plot 'Unitary3/populations_start_0.dat' u 1:2w l,'Unitary3/populations_start_0.dat' u 1:3w l,'Unitary3/populations_start_0.dat' u 1:4w l,'Unitary3/populations_start_0.dat' u 1:5w l,'Unitary3/populations_start_0.dat' u 1:6w l,'Unitary3/populations_start_0.dat' u 1:7w l,'Unitary3/populations_start_0.dat' u 1:8w l,'Unitary3/populations_start_0.dat' u 1:9w l,'Unitary3/populations_start_0.dat' u 1:10w l,'Unitary3/populations_start_0.dat' u 1:11w l,'Unitary3/populations_start_0.dat' u 1:12w l,'Unitary3/populations_start_0.dat' u 1:13w l,'Unitary3/populations_start_0.dat' u 1:14w l,'Unitary3/populations_start_0.dat' u 1:15w l,'Unitary3/populations_start_0.dat' u 1:16w l,'Unitary3/populations_start_0.dat' u 1:17w l
set output
quit
