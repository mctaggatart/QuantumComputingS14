for i=1:12
    dim = 2^i;
    A=rand(dim,dim);
    B=rand(dim,dim);
    tic;
    C=A*B; 
    e = 1000.0*toc;
    data(i)=e;
    dims(i)=i;
end
temp = [data;dims];
temp = temp'
save TestMatlab.dat -ascii -double temp