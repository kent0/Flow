clear ps;

n=round(size(csnap,1)*.5);
uu=csnap(1:n,:);

[rho,ps(1)]=max(abs(uu(1:n,1)));
U=uu(:,1);
I=speye(n);
P=I(:,ps(1));

for l=2:nb+1
    c=(P'*U(1:n,:))\(P'*uu(1:n,l));
    r=uu(1:n,l)-U(1:n,:)*c;
    [rho,ps(l)]=max(abs(r));
    U=[U uu(:,l)];
    P=[P I(:,ps(l))];
end

PUi=inv(P'*U);
