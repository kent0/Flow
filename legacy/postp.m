function postp(Nx,Ny,T,dt,cname,u0);

clear advance Einv;

nu=1./40;
nu=1./15000;

ifconv=true;

bc='dddd';
nf=round(Nx*.1);
nf=0;
df=.05;

nsteps=round(T/dt);
dt=T/nsteps;
iostep=round(.2/dt);
t=0.;

ifplot=true;

[exact,f,dom,ifexact]=problem(cname,nu,Ny);

[Abx,Aby,Bbx,Bby,Bbxm,Bbym,Bbpx,Bbpy, ...
    Dbx,Dby,Dbxm,Dbym,Fx,Fy,Ibx,Iby,Ibxm,Ibym, ...
    Jxm,Jym,Jxml,Jyml,Jxl,Jyl,Jpx,Jpy,Jpxl,Jpyl,Rx,Ry, ...
    Xn,Yn,Xm,Ym,Xl,Yl,Xpn,Ypn]=setops(Nx,Ny,dom,bc,nf,df);

x=reshape(Xn,[],1); y=reshape(Yn,[],1);

[Ax,Ay,Bx,By,Dx,Dy,Ix,Iy]=rops(Abx,Aby,Bbx,Bby,Dbx,Dby,Rx,Ry);

[Bxi,Byi,Bbpxi,Bbpyi]=invops(Bx,By,Bbpx,Bbpy);

ubsf=y*0;

[B,Bb,Binv,Hb,Hinv,DbT,DT,Db,D,Einv,C,R,RT, ...
   vort,sf,maskc]=setfuns(Abx,Aby,Ax,Ay, ...
    Bbx,Bby,Bx,By,Bxi,Byi,Bbxm,Bbym,Dbx,Dby,Dbxm,Dbym, ...
    Ibx,Iby,Ibxm,Ibym,Jxm,Jym,Jpx,Jpy,Fx,Fy,Rx,Ry,nu,ubsf);

if ~ifconv; C=@(uf) R(uf)*0; end

pdim=2;

[p1,p2,pm,pv,psf,pp,pdiv,plot_m1,plot_m2,plot_mm]=setplots(Xn,Yn,Xm,Ym,Xpn,Ypn,Jxl,Jyl,Jxml,Jyml,Jpxl,Jpyl,Xl,Yl,pdim,vort,sf,Db);

t=0;
ue=exact(x,y,t);
ub=maskc(ue);
p=reshape(Xpn,[],1)*0; dp=p;

if length(u0) < 2;
    uf=ue;
else
    Nx_old=size(u0,1)-1;
    Ny_old=size(u0,2)-1;
    uf=interp_uv(reshape(u0,[],1),Nx_old,Ny_old,Nx,Ny);
end

ndump=1500
%for idump=1:ndump
for idump=870:871
    dstr=sprintf('%06d',idump);
    dfname=[dstr,'.nx100.ny100.dat'];
%   ux=dlmread(dfname);
%   ux=reshape(ux,Nx+1,Ny+1,4);
%   uf=reshape(ux(:,:,1:2),[],1);
    uf=dlmread(dfname);
    time=idump*0.2;
    pm(uf); title(['Velocity Magnitude  TIME=' num2str(time,'%8.3e')],'FontName','Operator Mono');
    xlabel('x'); ylabel('y');
    set(gcf,'PaperUnits','inches','PaperPosition',[0 0 12.8 7.2])
    saveas(gcf,['umag.',dstr,'.png']);
end
