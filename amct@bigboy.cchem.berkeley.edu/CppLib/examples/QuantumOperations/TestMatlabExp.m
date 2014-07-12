for i=1:12
    dim = 2^i;
    A=rand(dim,dim)+1i*rand(dim,dim);
    A=A+A';
    rho=rand(dim,dim)+1i*rand(dim,dim);
    rho=rho+rho';
    tic;
	c=trace(A*rho);
	e = 1000.0*toc;
    data(i)=e;
    dims(i)=i;
end
temp = [data;dims];
temp = temp'
save TestMatlabExp.dat -ascii -double temp