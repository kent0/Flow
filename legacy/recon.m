function uf=recon(u,zetas)

nb=length(u)-1;
uf=zetas(:,1:nb+1)*u;
