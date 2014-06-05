set xlabel "Dimensions 2^n"
set logscale y
set key left top
set ylabel "Time [ms]"


plot "TestExpectation.dat" u 1:2 title 'Exp CBLAS' w l, "TestExpectation.dat" u 1:3 title 'Exp ForLoop' w l, "TestExpectation.dat" u 1:4 title 'Exp Implemented' w l, "TestMatlabExp.dat" u 2:1 title 'Exp - matlab' w l
