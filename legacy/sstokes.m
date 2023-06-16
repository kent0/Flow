close all
clear all
nu=1;
N=15;

Nx=N;
Ny=N;

exact=@(x,y,t) [1-reshape(y.*y,[],1);reshape(y,[],1)*0];
exact=@(x,y,t) [reshape(y.^Ny,[],1);reshape(y,[],1)*0];
exact=@(x,y,t) [reshape(y,[],1)*0;reshape(y,[],1)*0];
ax=0; bx=1; ay=-1; by=1;

[Ax,Ay,Abx,Aby,Bx,By,Bxi,Byi,Bbx,Bby,Bbxm,Bbym,Bbpxi,Bbpyi, ...
          Dbx,Dby,Dbxm,Dbym,Dbpx,Dbpy,Ibx,Iby,Ibxm,Ibym,Ibpx,Ibpy, ...
          Jxm,Jym,Jxl,Jyl,Jpx,Jpy,Jpxl,Jpyl,Rx,Ry, ...
          Xn,Yn,Xm,Ym,Xl,Yl,Xpn,Ypn]=setup(Nx,Ny,ax,bx,ay,by);

ue=exact(Xn,Yn,0);

ainv=@(b) a2dinv(b,Ax,Ay,Bx,By,Rx,Ry);
aop=@(ub) a2u(Ry,Rx,a2d(ub,Abx,Aby,Bbx,Bby));

uf=exact(Xn,Yn,0); n2=round(length(uf)/2);
u=a2u(Ry,Rx,uf);
ub=uf-a2u(Ry',Rx',a2u(Ry,Rx,uf));
p=reshape(Xpn*0,[],1);
f=uf*0; f=[f(1:n2)+2;f(n2+1:end)];

DT=@(p) dt2d(p,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,true);
Db=@(uv) d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,false);
D=@(uv) d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,true);

S=@(p) D(ainv(DT(p)));
M=@(x) a2u(Bbpyi,Bbpxi,x);

us=ainv(a2u(Ry*Bby,Rx*Bbx,f)-aop(ub));
res=-Db(ub+a2u(Ry',Rx',us));
p=pcg(S,res,1.e-10,1000,M);
uf=a2u(Ry',Rx',us+ainv(DT(p)))+ub;

umag=sqrt(uf(1:n2).^2+uf(n2+1:end).^2);

figure; mplot_ref(umag,Xn,Yn,Jxl,Jyl,Xl,Yl); title('||u||');
figure; mplot_ref(uf(1:n2),Xn,Yn,Jxl,Jyl,Xl,Yl); title('u');
figure; mplot_ref(uf(n2+1:end),Xn,Yn,Jxl,Jyl,Xl,Yl); title('v');
figure; mplot_ref(p,Xpn,Ypn,Jpxl,Jpyl,Xl,Yl); title('p');
figure; mplot_ref(Db(uf),Xpn,Ypn,Jpxl,Jpyl,Xl,Yl); title('D(u,v)');
