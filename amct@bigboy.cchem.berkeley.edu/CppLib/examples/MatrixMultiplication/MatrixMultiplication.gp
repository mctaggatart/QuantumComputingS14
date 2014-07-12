set xlabel "Dimensions 2^n"
set logscale y
set key left top
set ylabel "Time [ms]"


plot "TestComplex.dat" u 1:2 title 'Complex-Method CBLAS' w l, "TestComplex.dat" u 1:3 title 'Complex-Method Implemented' w l, "TestComplex.dat" u 1:4 title 'Complex-Method ForLoop' w l, "TestReal.dat" u 1:2 title 'Real-Method CBLAS' w l, "TestReal.dat" u 1:3  title 'Real-Method Implemented' w l, "TestReal.dat" u 1:4 title 'Real-Method ForLoop' w l, "TestMatlab.dat" u 2:1 title 'matlab-real' w l, "TestMatlabComplex.dat" u 2:1 title 'matlab-complex' w l

#plot "TestComplex.dat" u 1:3 title 'Complex-Method CBLAS' w l, "TestReal.dat" u 1:3 title 'Real-Method CBLAS' w l, "TestMatlab.dat" u 2:(1*$1) title 'matlab-real' w l, "TestMatlabComplex.dat" u 2:(1*$1) title 'matlab-complex' w l


