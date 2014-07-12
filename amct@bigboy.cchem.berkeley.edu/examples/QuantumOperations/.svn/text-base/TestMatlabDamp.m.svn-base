for i=1:8
    dim = 2^i;
    A=rand(dim,dim)+1i*rand(dim,dim);
    A=A+A';
    rho=rand(dim,dim)+1i*rand(dim,dim);
    rho=rho+rho';
    tic;
	nn=A'*rho;
	C=2.0*A*rho*A'-nn*rho-rho*nn;
    e = 1000.0*toc;
    data(i)=e;
    dims(i)=i;
end
temp = [data;dims];
temp = temp'
save TestMatlabDamp.dat -ascii -double temp