function [u,ttot,iits,ress,pmat,bmat,amat,c1mat,c2mat]=sns(Nx,Ny,tol,cname,nu,sig0,ptype,u0);

figure;

sigrat=0.5;
itmax=1000;
iits=zeros(itmax,1);
ress=zeros(itmax,1);

clear advance Einv DT;

%nu=1./7500;
%nu=1./200;

ifconv=true;

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

if ~ifconv; C=@(vf,uf) R(uf)*0; end

pdim=2;

[p1,p2,pm,pv,psf,pp,pdiv,plot_m1,plot_m2,plot_mm,plot_m3]=setplots(Xn,Yn,Xn2,Yn2,Xm,Ym,Xpn,Ypn,Jxnl,Jynl,Jxn2l,Jyn2l,Jxml,Jyml,Jxpl,Jypl,Xl,Yl,pdim,vort,sf,Db);

t=0;
ue=exact(x,y,t);
if strcmp(cname,'kov'); pe=kov_p(xp,yp,nu); end
if strcmp(cname,'walsh'); pe=walsh_p(xp,yp,0,nu); end
ub=maskc(ue);
p=xp*0; dp=p;

if length(u0) < 2;
    uf=ue;
else
    if size(u0,2) == 1
        N_old=round(sqrt(size(u0,1)*.5)-1);
        uf=interp_uv(reshape(u0,[],1),N_old,N_old,Nx,Ny);
    else
        Nx_old=size(u0,1)-1;
        Ny_old=size(u0,2)-1;
        uf=interp_uv(reshape(u0,[],1),Nx_old,Ny_old,Nx,Ny);
    end
end

uf=uf-maskc(uf)+ub;
uf=incomp2(uf,Binv,Einv,Db,DT,RT);

nsteps=1000;
iostep=1;
t=0.;

uf_old=uf;

ttot=0.;
count=0;

ubsf=uf*0;
ee=1;
S=[];
s=[];

dt=1;
ue_next=exact(x,y,t+dt); ub=maskc(ue_next);
sigma=0;
it=0;

c1mat=[];
c2mat=[];
amat=[];
bmat=[];
pmat=[];

istep=0;
while ee > tol
    it=it+1;
    tic;

%   [ee,S]=sadvance(uf,ub,S,Hb,Hinv,Bb,Binv,Einv,Db,DT,R,RT,C,sigma,pm);
    [ee,s,S,iit,pmat,bmat,amat,c1mat,c2mat]=sadvance2(uf,S,Hb,Hinv,Bb,Binv,Einv,Db,DT,R,RT,C,pm,ptype);

    ttot=ttot+toc;
    delta=norm(s)/norm(uf);
%   sigma=sigma*sigrat+exp(-1.5*delta)*(1-sigrat)
    sigma=sigma*sigrat+exp(-sig0*delta)*(1-sigrat)
%   sigma=1.;

    uf=uf+sigma*RT(s);

    [b1,b2]=evalf(uf,C,Hb,Binv,Einv,Db,DT,R,RT);
    iits(it)=iit;
    ress(it)=norm(R(incomp2(RT(b1+b2),Binv,Einv,Db,DT,RT)));

    t=t+dt;istep=istep+1;

    disp(ee)

    ue=ue_next;
    ud=uf_old-uf;
    uf_old=uf;
%   if it > 10; break; end

    [diff,err]=post(uf,ue,ud,p,ubsf,t,dt,Xn,Yn,istep,iostep,nsteps,ifexact,ifplot, ...
        Bb,p1,p2,pm,pv,psf,pp,pdiv);
end

iits=iits(1:it);
ress=ress(1:it);

u=reshape(uf,[Nx+1 Ny+1 2]);

disp(['Solve time=',num2str(ttot)]);

aop=@(s) R(Hb(RT(s),0));
c1op=@(s) R(C(RT(s),uf));
c2op=@(s) R(C(uf,RT(s)));
bop=Binv;
acop=@(s) Binv(R(Hb(RT(s),0)+C(RT(s),uf)+C(uf,RT(s))));
pop=@(s)R(incomp2(RT(s),Binv,Einv,Db,DT,RT));
jop=@(s) pop(acop(s));

acmat=opmat(acop,length(s));
jmat=opmat(jop,length(s));
amat=opmat(aop,length(s));
c1mat=opmat(c1op,length(s));
c2mat=opmat(c2op,length(s));
bmat=opmat(bop,length(s));
pmat=opmat(pop,length(s));

