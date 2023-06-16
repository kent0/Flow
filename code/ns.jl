#function ns(Nx,Ny,T,dt,cname,nu,u0)
using Plots

#[u,err,ttot,dt]=
#N=33;
#cname="ldc";Nx=N;Ny=N;dt=1e-3;nu=1/1000;u0=[];

ifconv=true;

bc="dddd";

nf=Int(round(Nx*.1));
nf=0;
df=.05;

ifplot=true;
#ifplot=false;

(exact,f,dom,ifexact)=problem(cname,nu);

Abx,Aby,Bbx,Bby,Bbxm,Bbym,Bbxn2,Bbyn2,Bbpx,Bbpy,Dbx,Dby,Dbxm,Dbym,Fx,Fy,Ibx,Iby,Ibxm,Ibym,Jxnm,Jynm,Jxmn,Jymn,Jxnn2,Jynn2,Jxn2m,Jyn2m,Jxn2l,Jyn2l,Jxml,Jyml,Jxnl,Jynl,Jxpn,Jypn,Jxpl,Jypl,Rx,Ry,Xn,Yn,Xn2,Yn2,Xm,Ym,Xl,Yl,Xpn,Ypn=setops(Nx,Ny,dom,bc,nf,df);

x=reshape(Xn,:); y=reshape(Yn,:);
xp=reshape(Xpn,:); yp=reshape(Ypn,:);

(Ax,Ay,Bx,By,Dx,Dy,Ix,Iy)=rops(Abx,Aby,Bbx,Bby,Dbx,Dby,Rx,Ry);

(Bxi,Byi,Bbpxi,Bbpyi)=invops(Bx,By,Bbpx,Bbpy);

(Ab,B,Bb,Binv,Hb,Hinv,DbT,DT,Db,D,Einv,C,R,RT,
   vort,sf,maskc)=setfuns(Abx,Aby,Ax,Ay,
    Bbx,Bby,Bx,By,Bxi,Byi,Bbxm,Bbym,Dbx,Dby,Dbxm,Dbym,
    Ibx,Iby,Ix,Iy,Ibxm,Ibym,Jxnm,Jynm,Jxpn,Jypn,Fx,Fy,Rx,Ry,nu);

if ~ifconv; C=(vf,uf) -> R(uf)*0; end

pdim=2;

(p1,p2,pm,pv,psf,pp,pdiv,plot_m1,plot_m2,plot_mm,plot_m3)=setplots(Xn,Yn,Xn2,Yn2,Xm,Ym,Xpn,Ypn,Jxnl,Jynl,Jxn2l,Jyn2l,Jxml,Jyml,Jxpl,Jypl,Xl,Yl,pdim,vort,sf,Db);

#global t=0.0;
ue=exact(x,y,t);
if cname == "kov"; pe=kov_p(xp,yp,nu); end
if cname == "walsh"; pe=walsh_p(xp,yp,0,nu); end
ub=maskc(ue);
p=xp*0; dp=p;

if length(u0) < 2;
    uf=ue;
else
    if size(u0,2) == 1
        N_old=Int(round(sqrt(size(u0,1)*.5)-1));
        uf=interp_uv(reshape(u0,:),N_old,N_old,Nx,Ny);
    else
        Nx_old=size(u0,1)-1;
        Ny_old=size(u0,2)-1;
        uf=interp_uv(reshape(u0,:),Nx_old,Ny_old,Nx,Ny);
    end
end

nsteps=Int(round(T/dt));
dt=T/nsteps;
iostep=Int(round(5/dt));
#iostep=Int(round(0.1/dt));

uf_old=uf;

ttot=0.;
count=0;

ubsf=uf*0;

ulag=[];
clag=[];
plag=[];
tlag=[];

println("sum uf $(sum(uf))")

for istep=1:nsteps
    global uf,t,p,dt,uf_old,ub,ue,err,ulag,clag,plag,tlag;
#   println("time $t")
    ue_next=exact(x,y,t+dt); ub=maskc(ue_next);
    if bc=="pppp"; ub=ub*0; end

    usf,ps,bdti1=advance(
#   usf,ps,bdti1=advance(ulag,clag,plag,tlag,
        uf,ub,p,t,dt,istep,Hb,Hinv,Bb,Binv,DbT,R,RT,C);

    (uf,p)=incomp(usf,ps,bdti1,Binv,Einv,Db,DT,RT);

#   dudt=evalf(uf,C,Hb,Binv,Einv,Db,DT,R,RT);
#   ndudt=norm(dudt);
#   if ndudt < 1.e-4; break; end

    if (ifexact && istep<5); uf=ue_next; end

    ue=ue_next;
    ud=uf_old-uf;
    uf_old=uf;
    t=t+dt;

    (diff,err)=post(uf,ue,ud,p,ubsf,t,dt,Xn,Yn,istep,iostep,nsteps,
       ifexact,ifplot,Bb,p1,p2,pm,pv,psf,pp,pdiv);
end

u=reshape(uf,Nx+1,Ny+1,2);

println("Done!");

#end

