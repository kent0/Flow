n=8;
m=n/2;
maxit=500;
tol=1.e-8;

xs=rand(n*n*2,1);
J0=interpn(n);
J=kron(J0,J0);
H=acmat;
H=amat;
Hc=zeros(m*m*2,m*m*2);
A=amat;
Ac=zeros(m*m*2,m*m*2);

Di=1/max(eig(A));

xs=xs+(Di*Di*Di)*(A*(A*(A*xs)));
b=A*xs;

Hc(1:m*m,1:m*m)         = J'*H(1:n*n,1:n*n)*J;
Hc(m*m+1:end,1:m*m)     = J'*H(n*n+1:end,1:n*n)*J;
Hc(1:m*m,m*m+1:end)     = J'*H(1:n*n,n*n+1:end)*J;
Hc(m*m+1:end,m*m+1:end) = J'*H(n*n+1:end,n*n+1:end)*J;

Ac(1:m*m,1:m*m)         = J'*A(1:n*n,1:n*n)*J;
Ac(m*m+1:end,1:m*m)     = J'*A(n*n+1:end,1:n*n)*J;
Ac(1:m*m,m*m+1:end)     = J'*A(1:n*n,n*n+1:end)*J;
Ac(m*m+1:end,m*m+1:end) = J'*A(n*n+1:end,n*n+1:end)*J;

sigma=0.66;
nsmooth=3;
ndim=2;
emin=0;
emax=1;
stype=1;
if1=false;

x=H\b;
%if1=true;
z=b*0;
z=z+vcycle(b-H*z,H,Hc,A,Ac,J,sigma,nsmooth,ndim,emin,emax,stype,if1);
z=z+vcycle(b-H*z,H,Hc,A,Ac,J,sigma,nsmooth,ndim,emin,emax,stype,if1);
z=z+vcycle(b-H*z,H,Hc,A,Ac,J,sigma,nsmooth,ndim,emin,emax,stype,if1);
z=z+vcycle(b-H*z,H,Hc,A,Ac,J,sigma,nsmooth,ndim,emin,emax,stype,if1);
z=z+vcycle(b-H*z,H,Hc,A,Ac,J,sigma,nsmooth,ndim,emin,emax,stype,if1);

nn=length(b)/2;

rmag=sqrt(b(1:nn).^2+b(nn+1:end).^2);
xmag=sqrt(x(1:nn).^2+x(nn+1:end).^2);
xsmag=sqrt(xs(1:nn).^2+xs(nn+1:end).^2);
zmag=sqrt(z(1:nn).^2+z(nn+1:end).^2);

xc=linspace(0,1,sqrt(length(b)/2)+2);
xc=linspace(0,1,sqrt(length(b)/2));

[X,Y]=ndgrid(xc,xc);

rmag=reshape(rmag,sqrt(nn),sqrt(nn));
xmag=reshape(xmag,sqrt(nn),sqrt(nn));
zmag=reshape(zmag,sqrt(nn),sqrt(nn));
xsmag=reshape(xsmag,sqrt(nn),sqrt(nn));

mesh(X,Y,zmag);
title('zmag');
figure
mesh(X,Y,xmag);
title('xmag');
figure
mesh(X,Y,xsmag);
title('xsmag');
figure
mesh(X,Y,rmag);
title('rmag');
