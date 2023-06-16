function z=vcycles(r0,H,Hc,A,Ac,Di,J,sigma,nsmooth,ndim,emin,emax,stype,if1)

%stype=2;

r=r0;
u=r*0;
%D=diag(H);
%Di=1./D
%pause;

n=round(length(r0)^(1/ndim));
m=n/2;

% pre-smoothing
sigma=0.5;
sigma=4/3;
%sigma=1/3;

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
      s=sigma*(Di.*r);
      u=u+s;
      r=r-H*s;
  end
end

%if size(Hc,1) == 1 || if1
if true;
    rc=J'*r;
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

    Acc=Jc'*Ac*Jc;
    Hcc=Jc'*Hc*Jc;

    Dci=diag(1/diag(Hcc));
    Dci=Dci/max(eig(Dci*Hcc));

    rc=J'*r;

%   if1=true;
    s=vcycles(rc,Hc,Hcc,Ac,Acc,Dci,Jc,sigma,nsmooth,ndim,emin,emax,stype,if1);
end

u=u+J*s;
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
      s=sigma*(Di.*r);
      u=u+s;
      r=r-H*s;
  end
end
z=u;
