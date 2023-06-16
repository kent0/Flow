nu=.01;
Nx=21;
Ny=21;

[Ax,Ay,Abx,Aby,Bx,By,Bxi,Byi,Bbx,Bby,Bbxm,Bbym,Dbx,Dby,Dbxm,Dbym, ...
   Ibx,Iby,Ibxm,Ibym,Jx,Jy,Jpx,Jpy,Rx,Ry,Xn,Yn,Xm,Ym]=setup(Nx,Ny);

hinv=@(b,bdti,ifupd) h2dinv(b,Ax,Ay,Bx,By,Rx,Ry,nu,bdti,ifupd);
hop=@(ub,bdti) h2d(ub,Abx,Aby,Bbx,Bby,Rx,Ry,nu,bdti);
mass=@(u) a2u(By,Bx,u);

uf=zeros(2*(Nx+1)*(Ny+1),1);
n2=round(length(uf)/2);
uf=[reshape((.5*(Yn+1)).^Ny,[],1);reshape(Yn*0,[],1)];
uf=uf+rand(size(uf,1),1);
%uf=uf(1:n2);
u=a2u(Ry,Rx,uf);
ub=uf-a2u(Ry',Rx',a2u(Ry,Rx,uf));

DT=@(p) dt2d(p,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,true);
DbT=@(p) dt2d(p,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,true);
D=@(uv) d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,true);
Db=@(uv) d2d(uv,Bbx,Bby,Dbx,Dby,Ibx,Iby,Jpx,Jpy,Rx,Ry,false);

g=-a2u(Ry,Rx,a2d(ub,Abx,Aby,Bbx,Bby));

mplot_ref(uf(1:n2),Xn,Yn,Jx,Jy,Xm,Ym);

T=10;
t=0.;
nsteps=100;
dt=T/nsteps;
iostep=10;

for istep=1:nsteps
    c=u*0;
    p=u*0;
    [u,bdti1]=bdfext_step2(u,ub,hinv,hop,mass,c,p,t,dt,istep);
    if (mod(istep,iostep)==0)
        uf=a2u(Ry',Rx',u)+ub;
%       disp(size(uf))
%       pause
%       mplot_ref(uf(1:n2),Xn,Yn,Jx,Jy,Xm,Ym);
        mplot_ref(uf(n2+1:end),Xn,Yn,Jx,Jy,Xm,Ym);
        drawnow
        pause(.1);
    end
    t=t+dt;
end
