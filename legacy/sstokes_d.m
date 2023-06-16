close all
nu=1;
N=511;

Nx=N;
Ny=N;

[Ax,Ay,Abx,Aby,Bx,By,Bxi,Byi,Bbx,Bby,Bbxm,Bbym,Bbpxi,Bbpyi,Dbx,Dby, ...
    Dbxm,Dbym,Ibx,Iby,Ibxm,Ibym,Jxm,Jym,Jxl,Jyl,Jpx,Jpy,Jpxl,Jpyl,Rx,Ry, ...
    Xn,Yn,Xm,Ym,Xl,Yl,Xpn,Ypn]=setup(Nx,Ny);

hinv=@(b,bdti,ifupd) h2dinv(b,Ax,Ay,Bx,By,Rx,Ry,nu,bdti,ifupd);
hop=@(ub,bdti) h2d(ub,Abx,Aby,Bbx,Bby,Rx,Ry,nu,bdti);
mass=@(u) a2u(By,Bx,u);
ainv=@(b) hinv(b,0.,false);

uf=[reshape((.5*(Yn+1)).^(Ny),[],1);reshape(Yn*0,[],1)];
u=a2u(Ry,Rx,uf);
ub=uf-a2u(Ry',Rx',a2u(Ry,Rx,uf));
u=rand(size(u,1),1);
p=zeros((Nx-1)*(Ny-1),1);

n2=round(length(uf)/2);
umag0=sqrt(uf(1:n2).^2+uf(n2+1:end).^2);

DT=@(p) dt2d(p,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,true);
DbT=@(p) dt2d(p,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,true);
D=@(uv) d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,true);
Db=@(uv) d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,false);

g=Db(ub);

t=0;
dt=1.e15;

%us=hinv(-hop(ub,0),0.,true);
p=g*0;
c=u*0;
[us,ps,bdti1]=bdfext_step3(u,ub,hinv,hop,mass,DT,c,p,t,dt,1);
uf=a2u(Ry',Rx',us)+ub;
res=-(g+D(us));
%M=@(x) x;
M=@(x) a2u(Bbpyi,Bbpxi,x);
p=ssolve(D,ainv,DT,res,M);
uf=a2u(Ry',Rx',DT(p));
u=us+hinv(DT(p),0,false);
uf=a2u(Ry',Rx',u)+ub;

umag=sqrt(uf(1:n2).^2+uf(n2+1:end).^2);

figure; mplot_ref(umag,Xn,Yn,Jxl,Jyl,Xl,Yl);
%figure; mplot_ref(uf(1:n2),Xn,Yn,Jxl,Jyl,Xl,Yl);
%figure; mplot_ref(uf(n2+1:end),Xn,Yn,Jxl,Jyl,Xl,Yl);
%figure; mplot_ref(p,Xpn,Ypn,Jpxl,Jpyl,Xl,Yl);
