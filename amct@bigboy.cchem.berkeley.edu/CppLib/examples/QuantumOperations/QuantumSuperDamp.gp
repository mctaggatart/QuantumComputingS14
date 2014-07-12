set xlabel "Dimensions 2^n"
set logscale y
set key left top
set ylabel "Time [ms]"


plot "TestSuperDamp.dat" u 1:2 title 'Damp - forloop' w l, "TestSuperDamp.dat" u 1:3 title 'Damp - Blas  ' w l,  "TestSuperDamp.dat" u 1:4 title 'Damp - Implemented' w l, "TestSuperDamp.dat" u 1:5 title 'Damp - Blas multiplication' w l, "TestMatlabDamp.dat" u 2:1 title 'Damp - matlab' w l
