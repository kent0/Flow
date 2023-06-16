function [x,rs,it]=pcg(b,a,minv,itmax,tol,ue,nxyz)

n=length(b);
r=b;
z=minv(r);
p=z;
x=0*b;
rzold=r'*z;

xx=[];
yy=[];

if length(nxyz)==2
    xx=1:nxyz(1);
    yy=1:nxyz(2);
end

rs=zeros(itmax+1,1);
it=0;
rs(1)=1.;

for k=1:itmax
    w=a(p);
    alpha=(r'*z)/(p'*w);
    x=x+alpha*p;
    r=r-alpha*w;
    e=ue-x;
    rs(k+1)=norm(e)/norm(ue);
    if length(nxyz)==2
        mesh(xx,yy,reshape(e,nxyz(1),nxyz(2))); pause
    end
%   rs(k)=sqrt(r'*r);
    it=it+1;
    if rs(k) < tol
        break;
    end
    z=minv(r);
    rznew=r'*z;
    beta=rznew/rzold;
    rzold=rznew;
    p=z+p*beta;
end
if length(nxyz)==2
    e=ue-x;
    rs(k+1)=norm(e)/norm(ue);
    mesh(xx,yy,reshape(e,nxyz(1),nxyz(2)));
    pause
end
rs=rs(1:it);
