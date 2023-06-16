clear ps;

n=round(size(csnap,1)*.5);
uu=csnap;

[rho,ps(1)]=max(uu(1:n,1).^2+uu(1+n:end,1).^2);
U=uu(:,1);
I=speye(n);
P=I(:,ps(1));

for l=2:nb
    c=(P'*U(1:n,:))\(P'*uu(1:n,l));
    r=uu(:,l)-U*c;
    [rho,ps(l)]=max(r(1:n).^2+r(n+1:end).^2);
    U=[U uu(:,l)];
    P=[P I(:,ps(l))];
end

PUi=inv(P'*U);
