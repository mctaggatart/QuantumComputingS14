for i=1:11
    dim = 2^i;
    A=rand(dim,dim)+1i*rand(dim,dim);
    A=A+A';
    tic;
	c=expm(-A*1i);
	e = 1000.0*toc;
    data(i)=e;
    dims(i)=i;
end
temp = [data;dims];
temp = temp'
save TestMatlabExp.dat -ascii -double temp