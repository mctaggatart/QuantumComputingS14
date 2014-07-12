set xlabel "time, t"
set key left top
set ylabel "norm, purity, concurrence"
set samples 1000

J = 0.06*2*pi

plot "TestQubitXY.dat" u 1:2 w l title 'norm', "TestQubitXY.dat" u 1:3 w l title 'purity', "TestQubitXY.dat" u 1:4 w l title 'con', sin(2.0*x*J)**2

