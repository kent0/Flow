function setops(Nx,Ny,dom,bc,nf,df);

k=2;

ax=dom[1];
bx=dom[2];
ay=dom[3];
by=dom[4];

lx=bx-ax;
ly=by-ay;

Npx=Nx-k;
Npy=Ny-k;

Mx=Int(round((3*Nx+1)*.5+.1));
My=Int(round((3*Ny+1)*.5+.1));

Lx=Nx*4;
Ly=Nx*4;

zxn,wxn=zwgll(Nx);
zyn,wyn=zwgll(Ny);

wx=wxn;
wy=wyn;

zxn2,wxn2=zwgll(Nx*2);
zyn2,wyn2=zwgll(Ny*2);

zxm,wxm=zwgll(Mx);
zym,wym=zwgll(My);

zxl,wxl=zwgll(Lx);
zyl,wyl=zwgll(Ly);

zzx,_=zwgll(Nx-nf);
zzy,_=zwgll(Ny-nf);

xf=ax.+lx*(zzx)*.5;
yf=ay.+ly*(zzy)*.5;

zxpn,wxpn=zwgll(Npx);
zypn,wypn=zwgll(Npy);

#   [zxpn,wxpn]=zwgl(Npx+1);
#   [zypn,wypn]=zwgl(Npy+1);

wxn=lx*wxn*.5; wyn=ly*wyn*.5;
wxn2=lx*wxn2*.5; wyn2=ly*wyn2*.5;
wxm=lx*wxm*.5; wym=ly*wym*.5;
wxl=lx*wxl*.5; wyl=ly*wyl*.5;
wxpn=lx*wxpn*.5; wypn=ly*wypn*.5;

xn=ax.+lx*(zxn.+1)*.5; yn=ay.+ly*(zyn.+1)*.5;
xn2=ax.+lx*(zxn2.+1)*.5; yn2=ay.+ly*(zyn2.+1)*.5;
xm=ax.+lx*(zxm.+1)*.5; ym=ay.+ly*(zym.+1)*.5;
xl=ax.+lx*(zxl.+1)*.5; yl=ay.+ly*(zyl.+1)*.5;
xpn=ax.+lx*(zxpn.+1)*.5; ypn=ay.+ly*(zypn.+1)*.5;

Xn,Yn=ndgrid(xn,yn);
Xn2,Yn2=ndgrid(xn2,yn2);
Xm,Ym=ndgrid(xm,ym);
Xl,Yl=ndgrid(xl,yl);

Xpn,Ypn=ndgrid(xpn,ypn);

if bc[1] == 'p' && bc[2] == 'p'; i1=1; i2=1; end
if bc[3] == 'p' && bc[4] == 'p'; i3=1; i4=1; end

if bc[1] == 'd'; i1=2; elseif bc[1] == 'n'; i1=1; end
if bc[2] == 'd'; i2=1; elseif bc[2] == 'n'; i2=0; end
if bc[3] == 'd'; i3=2; elseif bc[3] == 'n'; i3=1; end
if bc[4] == 'd'; i4=1; elseif bc[4] == 'n'; i4=0; end

Rx=speye(Nx+1); Rx=Rx[i1:end-i2,:];
Ry=speye(Ny+1); Ry=Ry[i3:end-i4,:];

if bc[1] == 'p' && bc[2] == 'p'; Rx[1,end]=1.; end
if bc[3] == 'p' && bc[4] == 'p'; Ry[1,end]=1.; end

Bbx=spdiagm(wxn);
Bby=spdiagm(wyn);

Bbxn2=spdiagm(wxn2);
Bbyn2=spdiagm(wyn2);

Bbxm=spdiagm(wxm);
Bbym=spdiagm(wym);

Bbpx=spdiagm(wxpn);
Bbpy=spdiagm(wypn);

Dbx=deriv_mat(xn);
Dby=deriv_mat(yn);

Dbxm=deriv_mat(xm);
Dbym=deriv_mat(ym);

Abx=Dbx'*Bbx*Dbx; Abx=.5*(Abx+Abx');
Aby=Dby'*Bby*Dby; Aby=.5*(Aby+Aby');

Ibx=speye(Nx+1);
Iby=speye(Ny+1);

Ibxm=speye(Mx+1);
Ibym=speye(My+1);

Jxnm=interp_mat(xm,xn);
Jynm=interp_mat(ym,yn);

Jxmn=interp_mat(xn,xm);
Jymn=interp_mat(yn,ym);

Jxnn2=interp_mat(xn2,xn);
Jynn2=interp_mat(yn2,yn);

Jxn2m=interp_mat(xm,xn2);
Jyn2m=interp_mat(ym,yn2);

Jxn2l=interp_mat(xl,xn2);
Jyn2l=interp_mat(yl,yn2);

Jxml=interp_mat(xl,xm);
Jyml=interp_mat(yl,ym);

Jxnl=interp_mat(xl,xn);
Jynl=interp_mat(yl,yn);

Jxpn=interp_mat(xn,xpn);
Jypn=interp_mat(yn,ypn);

Jxpl=interp_mat(xl,xpn);
Jypl=interp_mat(yl,ypn);

Jxf=interp_mat(xf,xn);
Jyf=interp_mat(yf,yn);

Fx=filt(zxn,wx,nf,df);
Fy=filt(zyn,wy,nf,df);

return Abx,Aby,Bbx,Bby,Bbxm,Bbym,Bbxn2,Bbyn2,Bbpx,Bbpy,
          Dbx,Dby,Dbxm,Dbym,Fx,Fy,Ibx,Iby,Ibxm,Ibym,
          Jxnm,Jynm,Jxmn,Jymn,Jxnn2,Jynn2,Jxn2m,Jyn2m,
          Jxn2l,Jyn2l,Jxml,Jyml,Jxnl,Jynl,Jxpn,Jypn,Jxpl,Jypl,
          Rx,Ry,Xn,Yn,Xn2,Yn2,Xm,Ym,Xl,Yl,Xpn,Ypn
end

