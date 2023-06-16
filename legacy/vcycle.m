function z=vcycle(r0,H,Hc,A,Ac,J,sigma,nsmooth,ndim,emin,emax,stype,if1)

%stype=2;

r=r0;
u=r*0;
D=diag(H);
Di=1./D;

n=round((length(r0)/ndim)^(1/ndim));
m=n/2;

%z=zwgll(n+1);
%z=z(2:end-1);

%[X,Y]=ndgrid(z,z);
%
%size(X)
%size(r0)
%rmag=sqrt(r(1:n*n).^2+r(n*n+1:end).^2);
%umag=sqrt(u(1:n*n).^2+u(n*n+1:end).^2);
%mesh(X,Y,reshape(umag,n,n))
%norm(rmag)
%pause

% pre-smoothing
sigma=0.5;

if stype==1 % jacobi
  for ks=1:nsmooth  
      s=sigma*(Di.*r);
      u=u+s;
      r=r-H*s;
  end
elseif stype==2 % chebyshev
  if stype == 2; edel=.5*(emax-emin); end
  for ks=1:nsmooth  
      thk=(pi/nsmooth)*(ks-.5);
      lamk=emin+edel*(cos(thk)+1);
      sigma=1/lamk;
      s=sigma*Di*r;
      u=u+s;
      r=r-H*s;
  end
end

if size(Hc,1) == 2 || if1
    nn=length(r)/2;
    rc=[J'*r(1:nn);J'*r(nn+1:end)];
    s=Hc\rc;
else
    l=m/2;
    J1=[];
    J1=interpn(m);

    if ndim == 1
        Jc=J1;
    elseif ndim == 2
        Jc=kron(J1,J1);
    elseif ndim == 3
        Jc=kron(kron(J1,J1),J1);
    end

    Acc=zeros(l*l*2);
    Hcc=zeros(l*l*2);

    Acc(1:l*l,1:l*l)         = Jc'*Ac(1:m*m,1:m*m)*Jc;
    Acc(l*l+1:end,1:l*l)     = Jc'*Ac(m*m+1:end,1:m*m)*Jc;
    Acc(1:l*l,l*l+1:end)     = Jc'*Ac(1:m*m,m*m+1:end)*Jc;
    Acc(l*l+1:end,l*l+1:end) = Jc'*Ac(m*m+1:end,m*m+1:end)*Jc;

    Hcc(1:l*l,1:l*l)         = Jc'*Hc(1:m*m,1:m*m)*Jc;
    Hcc(l*l+1:end,1:l*l)     = Jc'*Hc(m*m+1:end,1:m*m)*Jc;
    Hcc(1:l*l,l*l+1:end)     = Jc'*Hc(1:m*m,m*m+1:end)*Jc;
    Hcc(l*l+1:end,l*l+1:end) = Jc'*Hc(m*m+1:end,m*m+1:end)*Jc;

    rc=[J'*r(1:n*n);J'*r(n*n+1:end)];

%   if1=true;
    s=vcycle(rc,Hc,Hcc,Ac,Acc,Jc,sigma,nsmooth,ndim,emin,emax,stype,if1);
end

    u=u+[J*s(1:m*m); J*s(m*m+1:end)];
    umag=sqrt(u(1:n*n).^2+u(n*n+1:end).^2);
    mesh(X,Y,reshape(umag,n,n));
r=r0-H*u;

% post-smoothing

if stype == 1 % jacobi
  for ks=1:nsmooth  
      s=sigma*(Di.*r);
      u=u+s;
      r=r-H*s;
  end
elseif stype == 2 % chebyshev
  for ks=1:nsmooth  
      thk=(pi/nsmooth)*(ks-.5);
      lamk=emin+edel*(cos(thk)+1);
      sigma=1/lamk;
      s=sigma*Di*r;
      u=u+s;
      r=r-A*s;
  end
end
z=u;
