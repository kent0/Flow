function Fx=filt(z,wz,nf,df)

N=length(z)-1;

Lz=legendre(z,N);
Bh=sparse(diag(wz));

for i=1:(N+1); Lz(:,i)=Lz(:,i)*sqrt(1/(Lz(:,i)'*Bh*Lz(:,i))); end
Fx=Lz*sparse(diag([ones(N+1-nf,1);1-df./(nf:-1:1)']))*Lz'*Bh;
