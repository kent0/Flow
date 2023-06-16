Ns=15:60;
n=length(Ns);
errn=zeros(n,1);

figure
for i=1:n
    N=Ns(i);
    disp(['N=',num2str(N)]);
    dt=4.e-4;
    [uf,errn(i),ttot,dt]=ns(N,N,.1,dt,'walsh',0);
end

figure; semilogy(Ns,errn);
pause

dts=4.e-4*1.2.^(0:20);
n=length(dts);
errt=zeros(n,1);

N=60;

figure
for i=1:n
    dt=dts(i);
    disp(['dt=',num2str(dt)]);
    [uf,errt(i),ttot,dt]=ns(N,N,.1,dt,'walsh',0);
    dts(i)=dt;
end

figure; loglog(dts,errt); hold on; loglog(dts,dts.^3);
