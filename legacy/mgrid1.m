
ifconv=true;

Nx=N;
Ny=N;
nu=1;

cname='ldc';
[exact,f,dom,ifexact]=problem(cname,nu,Ny);

bc='dddd';

nf=round(Nx*.1);
nf=0;
df=.05;

ifplot=true;
%ifplot=false;

[exact,f,dom,ifexact]=problem(cname,nu,Ny);

[Abx,Aby,Bbx,Bby,Bbxm,Bbym,Bbxn2,Bbyn2,Bbpx,Bbpy,Dbx,Dby,Dbxm,Dbym,Fx,Fy,Ibx,Iby,Ibxm,Ibym,Jxnm,Jynm,Jxmn,Jymn,Jxnn2,Jynn2,Jxn2m,Jyn2m,Jxn2l,Jyn2l,Jxml,Jyml,Jxnl,Jynl,Jxpn,Jypn,Jxpl,Jypl,Rx,Ry,Xn,Yn,Xn2,Yn2,Xm,Ym,Xl,Yl,Xpn,Ypn]=setops(Nx,Ny,dom,bc,nf,df);

x=reshape(Xn,[],1); y=reshape(Yn,[],1);
xp=reshape(Xpn,[],1); yp=reshape(Ypn,[],1);

[Ax,Ay,Bx,By,Dx,Dy,Ix,Iy]=rops(Abx,Aby,Bbx,Bby,Dbx,Dby,Rx,Ry);

[Bxi,Byi,Bbpxi,Bbpyi]=invops(Bx,By,Bbpx,Bbpy);

[Ab,B,Bb,Binv,Hb,Hinv,DbT,DT,Db,D,Einv,C,R,RT, ...
   vort,sf,maskc]=setfuns(Abx,Aby,Ax,Ay, ...
    Bbx,Bby,Bx,By,Bxi,Byi,Bbxm,Bbym,Dbx,Dby,Dbxm,Dbym, ...
    Ibx,Iby,Ix,Iy,Ibxm,Ibym,Jxnm,Jynm,Jxpn,Jypn,Fx,Fy,Rx,Ry,nu);

t=0.1;
u0=R(exact(x,y,t));
n=length(u0)/2;
nn=sqrt(n);

%[Sx,lamx]=eig(full(Ax),full(Bx));
%[Sy,lamy]=eig(full(Ay),full(By));
%D=diag(kron(Iy,sparse(lamx))+kron(sparse(lamy),Ix));
%emax=max(abs(D));

aop=@(s) R(Hb(RT(s),0));
H=opmat(aop,n);
H=kron(Ax,Iy)+kron(Ix,Ay);
H=kron(Ax,By)+kron(Bx,Ay);
A=1*H;
Di=1./diag(H);
Dif=diag(Di);

Di=Di/max(eig(Dif*H));
J1=interpn(nn);
J=kron(J1,J1);
Hc=J'*H*J;
Ac=J'*A*J;
x=rand(n,1);
b=H*x;
tol=1.e-8;
z=0*x;

%emax=max(eig(Ax))+max(eig(Ay));
%emax=2*max(eig(Ax))*max(eig(Bx));
%emin=emax/4;

emax=1;
emin=emax*efac;

itmax=1000;
res=zeros(itmax+1,1);
res(1)=1.;
it=0;

for i=1:itmax
    z=z+vcycles(b-H*z,H,Hc,A,Ac,Di,J,sigma,nsmooth,ndim,emin,emax,stype,if1);
    res(i+1)=norm(x-z)/norm(x);
    it=it+1;
    if res(i+1) < tol; break; end
end

res=res(1:it);
