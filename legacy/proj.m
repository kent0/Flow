function u=proj(uf,zetas,Op,nb);

u=zeros(nb,1);
ufm0=Op(uf-zetas(:,1));

u(1)=1;
for i=2:nb+1
   u(i)=dot(ufm0,zetas(:,i));
end
