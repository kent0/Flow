N=63;
Nx=N;
Ny=N;
cname='ldc_reg';
bc='dddd';
nf=0;
df=0;
ifconv=true;

[exact,f,dom,ifexact]=problem(cname,nu,Ny);

[Abx,Aby,Bbx,Bby,Bbxm,Bbym,Bbxn2,Bbyn2,Bbpx,Bbpy,Dbx,Dby,Dbxm,Dbym,Fx,Fy,Ibx,Iby,Ibxm,Ibym,Jxnm,Jynm,Jxmn,Jymn,Jxnn2,Jynn2,Jxn2m,Jyn2m,Jxn2l,Jyn2l,Jxml,Jyml,Jxnl,Jynl,Jxpn,Jypn,Jxpl,Jypl,Rx,Ry,Xn,Yn,Xn2,Yn2,Xm,Ym,Xl,Yl,Xpn,Ypn]=setops(Nx,Ny,dom,bc,nf,df);

x=reshape(Xn,[],1); y=reshape(Yn,[],1);

[Ax,Ay,Bx,By,Dx,Dy,Ix,Iy]=rops(Abx,Aby,Bbx,Bby,Dbx,Dby,Rx,Ry);

[Bxi,Byi,Bbpxi,Bbpyi]=invops(Bx,By,Bbpx,Bbpy);

ubsf=y*0;

[Ab,B,Bb,Binv,Hb,Hinv,DbT,DT,Db,D,Einv,C,R,RT, ...
   vort,sf,maskc]=setfuns(Abx,Aby,Ax,Ay, ...
   Bbx,Bby,Bx,By,Bxi,Byi,Bbxm,Bbym,Dbx,Dby,Dbxm,Dbym, ...
   Ibx,Iby,Ix,Iy,Ibxm,Ibym,Jxnm,Jynm,Jxpn,Jypn,Fx,Fy,Rx,Ry,nu);

if ~ifconv; C=@(uf) R(uf)*0; end

pdim=2;

[p1,p2,pm,pv,psf,pp,pdiv,plot_m1,plot_m2,plot_mm,plot_m3]=setplots(Xn,Yn,Xn2,Yn2,Xm,Ym,Xpn,Ypn,Jxnl,Jynl,Jxn2l,Jyn2l,Jxml,Jyml,Jxpl,Jypl,Xl,Yl,pdim,vort,sf,Db);

nb=20;
ndump=2000;
us=load_snaps(Nx,Ny,ndump,'r7500/%06d.nx63.ny63.dat');
uf=dlmread('r7500/002000.nx63.ny63.dat');
nu=1./7500;

pod;
podc;

genops;
