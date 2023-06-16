function [PUi,ps,U]=deim(uu,nc);

n=size(uu,1);

ps=[];

[rho,ps(1)]=max(abs(uu(:,1)));
U=uu(:,1);
I=speye(n);
P=I(:,ps(1));

for l=2:nc
    c=(P'*U)\(P'*uu(:,l));
    r=uu(:,l)-U*c;
    [rho,ps(l)]=max(abs(r));
    U=[U uu(:,l)];
    P=[P I(:,ps(l))];
end

PUi=inv(P'*U);
