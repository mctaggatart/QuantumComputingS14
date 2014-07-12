for i=1:10
    dim = 2^i;
    A=rand(dim,dim)+1i*rand(dim,dim);
    A=A+A';
    rho=rand(dim,dim)+1i*rand(dim,dim);
    rho=rho+rho';
    tic;
	c=-1i*A*rho+1i*rho*A;
	e = 1000.0*toc;
    data(i)=e;
    dims(i)=i;
end
temp = [data;dims];
temp = temp'
save TestMatlabHami.dat -ascii -double temp