Ns=3:24;
n=length(Ns);
err=zeros(n,1);

dt=1.e-3;

figure
for i=1:n
    N=Ns(i);
    disp(['N=',num2str(N)]);
    [uf,err(i),ttot]=ns(N,N,10.,dt,'kov',0);
end

figure; semilogy(Ns,err);
