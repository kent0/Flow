function u=proj2(uf,zetas,Op,nb);

u=zeros(nb,1);
vf=Op(uf);

for i=1:nb
   u(i)=dot(vf,zetas(:,i));
end
