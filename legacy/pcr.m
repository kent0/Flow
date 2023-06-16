function [x,rs,it]=pcr(b,a,minv,itmax,tol,nxyz)

n=length(b);

r=minv(b);
p=r;
x=0*b;

xx=[];
yy=[];

if length(nxyz)==2
end

ar=a(r);
w=ar;
b0=r'*ar;
bold=b0;

rs=zeros(itmax+1,1);
it=0;
rs(1)=1.;

for k=1:itmax
    v=minv(w);

    alpha=bold/(w'*v);

    x=x+alpha*p;
    r=r-alpha*v;
    ar=a(r);

    bnew=r'*ar;
    beta=bnew/bold;
    bold=bnew;

    p=r+beta*p;
    w=ar+beta*w;

    rs(k+1)=sqrt(bnew/bold);
    if length(nxyz)==2
        mesh(xx,yy,reshape(r,nxyz(1),nxyz(2))); pause
    end

    it=it+1;
    if rs(k) < tol; break; end
end

if length(nxyz)==2
    e=ue-x;
    rs(k+1)=norm(e)/norm(ue);
    mesh(xx,yy,reshape(e,nxyz(1),nxyz(2)));
    pause
end
rs=rs(1:it);
