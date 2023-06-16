close all
nu=1;
N=63;
%N=5;
Nx=N;
Ny=N;

[Ax,Ay,Abx,Aby,Bx,By,Bxi,Byi,Bbx,Bby,Bbxm,Bbym,Bbpxi,Bbpyi, ...
          Dbx,Dby,Dbxm,Dbym,Dbpx,Dbpy,Ibx,Iby,Ibxm,Ibym,Ibpx,Ibpy, ...
          Jxm,Jym,Jxl,Jyl,Jpx,Jpy,Jpxl,Jpyl,Rx,Ry, ...
          Xn,Yn,Xm,Ym,Xl,Yl,Xpn,Ypn]=setup(Nx,Ny);

hinv=@(b,bdti,ifupd) h2dinv(b,Ax,Ay,Bx,By,Rx,Ry,nu,bdti,ifupd);
hop=@(ub,bdti) h2d(ub,Abx,Aby,Bbx,Bby,Rx,Ry,nu,bdti);
mass=@(u) a2u(By,Bx,u);
ainv=@(b) hinv(b,0.,false);

uf=[reshape((.5*(Yn+1)).^(Ny),[],1);reshape(Yn*0,[],1)];
u=a2u(Ry,Rx,uf);
ub=uf-a2u(Ry',Rx',a2u(Ry,Rx,uf));
%sstokes;
%u=rand(size(u,1),1);
p=zeros((Nx-1)*(Ny-1),1);

n2=round(length(uf)/2);
umag0=sqrt(uf(1:n2).^2+uf(n2+1:end).^2);

DT=@(p) dt2d(p,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,true);
DbT=@(p) dt2d(p,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,false);
D=@(uv) d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,true);
Db=@(uv) d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,false);

g=Db(ub);

T=.1;
nsteps=200;
iostep=10;
dt=T/nsteps;
t=0.;


Bbpx=inv(Bbpxi);
Bbpy=inv(Bbpyi);

Abpx=Dbpx'*Bbpxi*Dbpx; Abpx=.5*(Abpx+Abpx');
Abpy=Dbpy'*Bbpyi*Dbpy; Abpy=.5*(Abpy+Abpy');

[Sx,lamx]=eig(Abpx,full(Bbpxi));
[Sy,lamy]=eig(Abpy,full(Bbpyi));

Abpx=Dbpx'*Bbpx*Dbpx; Abpx=.5*(Abpx+Abpx');
Abpy=Dbpy'*Bbpy*Dbpy; Abpy=.5*(Abpy+Abpy');

[Sx,lamx]=eig(Abpx,full(Bbpx));
[Sy,lamy]=eig(Abpy,full(Bbpy));

Dia=diag(kron(sparse(lamy),Ibpx)+kron(Ibpy,sparse(lamx)));
Dinv=1./Dia; Dinv(1)=0.;
%disp(sort(Dinv));

M=@(x) a2u(inv(Bbpyi),inv(Bbpxi),x);
M=@(x) a2u(Bbpyi,Bbpxi,x);
M=@(x) a2u(Sy,Sx,Dinv.*a2u(Sy',Sx',x));
%M=@(x) x;

for istep=1:nsteps
    c=u*0;
    [us,ps,bdti1]=bdfext_step3(u,ub,hinv,hop,mass,DT,c,p,t,dt,istep);
    res=-(g+D(us));
    if (istep<4); hinv2=@(x) hinv(x,bdti1,false); end
    dp=dp*0;
    dp=ssolve(D,hinv2,DT,res,M,dp);
    p=dp+ps;
    u=us+hinv(DT(dp),bdti1,false);
    uf=a2u(Ry',Rx',u)+ub;
    t=t+dt;
    if mod(istep,iostep)==0
        umag=sqrt(uf(1:n2).^2+uf(n2+1:end).^2);
        mplot_ref(umag,Xn,Yn,Jxl,Jyl,Xl,Yl);
%       mplot_ref(uf(1:n2),Xn,Yn,Jxl,Jyl,Xl,Yl);
%       mplot_ref(uf(n2+1:end),Xn,Yn,Jxl,Jyl,Xl,Yl);
%       mplot_ref(p,Xpn,Ypn,Jpxl,Jpyl,Xl,Yl);
        pause(.1);
    end
end
