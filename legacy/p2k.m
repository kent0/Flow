function uk=p2k(snaps,z,Op)

nsnap=size(snaps,2);
nb=size(z,2)-1;

uk=zeros(nb,nsnap);

for i=1:nsnap
    Opui=Op(snaps(:,i));
    for j=2:nb+1
        uk(j-1,i)=dot(Opui,z(:,j));
    end
end
